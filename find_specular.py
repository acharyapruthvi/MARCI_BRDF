#!/usr/bin/env python3
"""
find_specular.py - scan every pixel of one or more ISIS cubes for
SPECULAR-REFLECTION viewing geometry and write the matches to a CSV.

This is stage 2 of the two-stage specular pipeline (stage 1 is
make_specular_geom.sh, which generates one multi-band phocube geometry cube
<basename>.geom.cub per source cube).

The specular condition, per pixel:
  1. localemission == localincidence      (within --tol-angle degrees)
  2. relative azimuth between the spacecraft and the sun == 180 degrees
     (within --tol-az degrees), i.e. the spacecraft sits on the OPPOSITE
     side of the local surface normal from the sun - the mirror direction
     in the principal plane.
Additionally, pixels with localemission or localincidence below --min-angle
are excluded: as both angles approach zero the azimuths become numerically
ill-defined and condition (2) is meaningless, so near-nadir/near-noon pixels
would otherwise flood the output with false "specular" hits.

The relative azimuth is computed as the absolute difference of the two
phocube azimuth backplanes, wrapped into [0, 180]:
  d      = |spacecraftazimuth - sunazimuth|  (mod 360)
  relaz  = min(d, 360 - d)
so the specular test is |relaz - 180| <= --tol-az.

Bands of the geom cube are identified BY GDAL BAND DESCRIPTION
(GetDescription(), the "Description = ..." string gdalinfo prints), never by
position, because phocube emits bands in its own canonical order - not the
order requested on its command line. A geom cube missing any required
backplane fails LOUDLY for that cube instead of silently mislabeling.

DN comes from the ORIGINAL cube at the same (line, sample); pixels whose DN
is an ISIS special pixel (Null/Lis/Lrs/His/Hrs) or NaN are excluded, as are
pixels where any required backplane is special/NaN.

NOTE ON OUTPUT VOLUME: specular hits come in contiguous PATCHES of adjacent
pixels, not isolated points - every pixel inside the tolerance window is a
row. Use --max-per-cube to keep only the N most-specular pixels per cube
(ranked by a normalized deviation score; see spec_score column), or
post-filter/cluster the CSV afterwards.

Output CSV columns (fixed, in this order):
  filename, ls, line, sample, dn, lat, lon, phase,
  localemission, localincidence, ei_diff,
  spacecraftazimuth, sunazimuth, relaz, offnadir, spec_score

  ei_diff    = |localemission - localincidence|            (deg)
  relaz      = relative spacecraft-sun azimuth in [0,180]  (deg)
  spec_score = hypot(ei_diff/tol_angle, (180-relaz)/tol_az); 0 = perfectly
               specular, 1 = on the edge of the tolerance window. Rows are
               sorted by this within each cube.
  line, sample are 0-BASED array indices into the (possibly camtrim-nulled)
  cube. Ls is parsed from the cube filename (3rd underscore field, last
  digit = decimal place); use --no-ls if filenames don't follow that scheme.

Usage:
  python find_specular.py --geom-dir /data/specular_geom -o specular.csv /data/MARCI_IMAGES/DONE/*/*.cub
  python find_specular.py --tol-angle 0.5 --tol-az 2 -o specular.csv myframe.l3.cub
  python find_specular.py --jobs 36 --max-per-cube 50 -o specular.csv ...

By default the geom cube is looked for alongside each original cube; use
--geom-dir if make_specular_geom.sh wrote them elsewhere (-d).
"""

import argparse
import csv
import glob
import math
import os
import re
import sys

try:
    from osgeo import gdal
    gdal.UseExceptions()
except ImportError:
    gdal = None

import numpy as np

# ISIS special-pixel DN values, matched defensively by value.
ISIS_SPECIAL_8BIT = (0, 255)
ISIS_SPECIAL_16BIT = (-32768, -32767, -32766, -32765, -32764)

BACKPLANES = ["lat", "lon", "phase", "localemission", "localincidence",
              "spacecraftazimuth", "sunazimuth", "offnadir"]

# Map a band's NORMALIZED GDAL description to its canonical backplane name.
# phocube does NOT emit bands in command-line flag order, so this description
# match - not band position - is what tells us which band is which.
DESC_TO_NAME = {
    "latitude": "lat",
    "longitude": "lon",
    "phaseangle": "phase",
    "localemissionangle": "localemission",
    "localincidenceangle": "localincidence",
    "sunazimuth": "sunazimuth",
    "spacecraftazimuth": "spacecraftazimuth",
    "offnadirangle": "offnadir",
}
for _n in BACKPLANES:
    DESC_TO_NAME.setdefault(_n, _n)

GEOM_SUFFIX = "geom"

COLUMNS = ["filename", "ls", "line", "sample", "dn", "lat", "lon", "phase",
           "localemission", "localincidence", "ei_diff",
           "spacecraftazimuth", "sunazimuth", "relaz", "offnadir",
           "spec_score"]


def err(*a):
    print("find_specular.py:", *a, file=sys.stderr)


def _norm_desc(s):
    """Normalize a band description for matching: lowercase, drop non-alnum."""
    return re.sub(r"[^a-z0-9]", "", (s or "").lower())


def special_mask(arr, dtype_name):
    """Boolean array: True where a pixel is an ISIS special pixel or NaN."""
    m = np.isnan(arr)
    if "Float" in dtype_name:
        # ISIS float special pixels are huge-magnitude sentinels.
        m |= (arr <= -3.4e38) | (arr >= 3.4e38)
    elif "Int16" in dtype_name:
        m |= np.isin(arr, ISIS_SPECIAL_16BIT)
    elif "Byte" in dtype_name:
        m |= np.isin(arr, ISIS_SPECIAL_8BIT)
    return m


def read_band_array(path, band=1):
    """(numpy array in NATIVE dtype, dtype_name) for one band of a cube.
    No float64 promotion - float32 halves memory and is ample for angles."""
    ds = gdal.Open(path)
    if ds is None:
        raise FileNotFoundError(f"GDAL could not open: {path}")
    bnd = ds.GetRasterBand(band)
    dtype_name = gdal.GetDataTypeName(bnd.DataType)
    arr = bnd.ReadAsArray()
    ds = None
    return arr, dtype_name


def geom_path(cube, geom_dir=None):
    base = cube[:-4] if cube.endswith(".cub") else cube
    bn = os.path.basename(base)
    d = geom_dir if geom_dir else os.path.dirname(cube)
    return os.path.join(d, f"{bn}.{GEOM_SUFFIX}.cub") if d else f"{bn}.{GEOM_SUFFIX}.cub"


def parse_ls(cube):
    """Ls from the cube basename: 3rd underscore field, LAST digit is the
    decimal place (B11_013726_2951_... -> 295.1). Raises ValueError if the
    filename does not follow the scheme (hard-fail by design)."""
    bn = os.path.basename(cube)
    parts = bn.split("_")
    if len(parts) < 3 or not re.fullmatch(r"[0-9]{2,}", parts[2]):
        raise ValueError(
            f"cannot parse Ls from filename '{bn}' (expected a numeric 3rd "
            f"underscore field, e.g. B11_013726_2951_...); use --no-ls to "
            f"skip Ls parsing")
    f = parts[2]
    return float(f[:-1] + "." + f[-1])


def scan_cube(cube, geom_dir, tol_angle, tol_az, min_angle, max_per_cube,
              want_ls, lat_range=None, ei_only=False):
    """Scan ONE cube for matching pixels.

    Two modes:
      specular (default): |localemission - localincidence| <= tol_angle AND
        the spacecraft-sun relative azimuth within tol_az of 180 deg.
      ei_only: ONLY |localemission - localincidence| <= tol_angle - the
        equal-angle locus regardless of azimuth. The azimuth bands are then
        not read during the search at all (less I/O); they are still sampled
        at hit pixels so the CSV columns stay identical, and spec_score
        becomes just ei_diff/tol_angle.

    Two-stage lazy read: only the bands the condition needs (plus latitude
    if lat_range is given) are read up front; DN and the remaining output
    bands are read ONLY if the condition matched somewhere. Most cubes have
    zero hits, so this skips more than half the I/O.

    lat_range: optional (min_lat, max_lat) - pixels outside it are excluded
    (replaces the old camtrim crop, using phocube's own latitude backplane).

    Returns (rows, reason) - rows is a list of row dicts (possibly empty),
    reason is a skip/diagnostic string or None. Raises ValueError on an
    unparseable Ls filename when want_ls (hard-fail: aborts the whole run)."""
    ls_val = parse_ls(cube) if want_ls else ""

    geom = geom_path(cube, geom_dir)
    if not os.path.isfile(geom):
        return [], f"missing geom cube: {os.path.basename(geom)}"

    ds = gdal.Open(geom)
    if ds is None:
        return [], f"GDAL could not open geom cube: {os.path.basename(geom)}"

    # Map backplane name -> band index BY DESCRIPTION (never position;
    # phocube emits its own canonical order). No pixel I/O yet.
    idx_of, dtype_of, descs = {}, {}, []
    for b in range(1, ds.RasterCount + 1):
        bnd = ds.GetRasterBand(b)
        desc = bnd.GetDescription()
        descs.append(desc)
        name = DESC_TO_NAME.get(_norm_desc(desc))
        if name is not None and name not in idx_of:
            idx_of[name] = b
            dtype_of[name] = gdal.GetDataTypeName(bnd.DataType)

    missing = [n for n in BACKPLANES if n not in idx_of]
    if missing:
        return [], (f"geom cube missing backplane(s) {missing}; "
                    f"band descriptions found: {descs}")

    def band(name):
        """Read one backplane in its NATIVE dtype (float32 for phocube)."""
        return ds.GetRasterBand(idx_of[name]).ReadAsArray()

    # --- STAGE 1: read ONLY what the condition needs -------------------------
    emi = band("localemission")
    inc = band("localincidence")

    bad = special_mask(emi, dtype_of["localemission"])
    bad |= special_mask(inc, dtype_of["localincidence"])

    scaz = suaz = None
    if not ei_only:
        scaz = band("spacecraftazimuth")
        suaz = band("sunazimuth")
        bad |= special_mask(scaz, dtype_of["spacecraftazimuth"])
        bad |= special_mask(suaz, dtype_of["sunazimuth"])

    lat = None
    if lat_range is not None:
        lat = band("lat")
        bad |= special_mask(lat, dtype_of["lat"])

    with np.errstate(invalid="ignore"):
        ei_diff = np.abs(emi - inc)
        hit = (~bad) \
            & (ei_diff <= tol_angle) \
            & (emi >= min_angle) \
            & (inc >= min_angle)
        if not ei_only:
            d = np.abs(scaz - suaz) % 360.0
            relaz = np.minimum(d, 360.0 - d)      # in [0, 180]
            az_dev = 180.0 - relaz                # 0 = perfectly opposed
            hit &= (az_dev <= tol_az)
        if lat_range is not None:
            hit &= (lat >= lat_range[0]) & (lat <= lat_range[1])

    lines, samples = np.nonzero(hit)
    if lines.size == 0:
        ds = None
        return [], None

    # --- STAGE 2: hits exist - now (and only now) read DN and the rest ------
    dn_arr, dn_dtype = read_band_array(cube)
    if dn_arr.shape != emi.shape:
        ds = None
        return [], (f"shape mismatch: DN {dn_arr.shape} vs geom {emi.shape} "
                    f"(was the geom cube made from a different cube?)")

    # drop hits whose DN is special/NaN
    keep = ~special_mask(dn_arr[lines, samples], dn_dtype)
    lines, samples = lines[keep], samples[keep]
    if lines.size == 0:
        ds = None
        return [], None

    # In ei_only mode the azimuths weren't read in stage 1 - read them now
    # (output columns only). relaz at hit pixels; NaN where special.
    if ei_only:
        scaz = band("spacecraftazimuth")
        suaz = band("sunazimuth")
        with np.errstate(invalid="ignore"):
            d = np.abs(scaz - suaz) % 360.0
            relaz = np.minimum(d, 360.0 - d)
        az_bad = special_mask(scaz, dtype_of["spacecraftazimuth"]) \
            | special_mask(suaz, dtype_of["sunazimuth"])
        relaz = relaz.astype("float64")
        relaz[az_bad] = math.nan

    # normalized deviation score: 0 = perfect geometry
    if ei_only:
        score = ei_diff[lines, samples].astype("float64") / tol_angle
    else:
        score = np.hypot(ei_diff[lines, samples].astype("float64") / tol_angle,
                         az_dev[lines, samples].astype("float64") / tol_az)
    order = np.argsort(score, kind="stable")
    if max_per_cube is not None and order.size > max_per_cube:
        order = order[:max_per_cube]
    lines, samples, score = lines[order], samples[order], score[order]

    # values for the output columns, sampled at hit indices only; special
    # pixels become NaN. Stage-1 bands are reused, the rest read here.
    stage1 = {"localemission": emi, "localincidence": inc,
              "spacecraftazimuth": scaz, "sunazimuth": suaz}
    if lat is not None:
        stage1["lat"] = lat
    col_vals = {}
    for name in BACKPLANES:
        arr = stage1.get(name)
        if arr is None:
            arr = band(name)
        v = arr[lines, samples].astype("float64")
        v[special_mask(v, dtype_of[name])] = math.nan
        col_vals[name] = v
    ds = None

    bn = os.path.basename(cube)
    rows = []
    for i, (ln, sp, sc) in enumerate(zip(lines, samples, score)):
        row = {"filename": bn, "ls": ls_val,
               "line": int(ln), "sample": int(sp),
               "dn": float(dn_arr[ln, sp]),
               "ei_diff": float(ei_diff[ln, sp]),
               "relaz": float(relaz[ln, sp]),
               "spec_score": float(sc)}
        for name in BACKPLANES:
            row[name] = float(col_vals[name][i])
        rows.append(row)
    return rows, None


def _worker(task):
    """Multiprocessing worker: returns (cube, rows, reason, error_str)."""
    cube, cfg = task
    try:
        rows, reason = scan_cube(cube, cfg["geom_dir"], cfg["tol_angle"],
                                 cfg["tol_az"], cfg["min_angle"],
                                 cfg["max_per_cube"], cfg["want_ls"],
                                 cfg.get("lat_range"),
                                 cfg.get("ei_only", False))
        return cube, rows, reason, None
    except ValueError as e:
        # Ls parse failure: hard-fail the whole run, per pipeline convention.
        return cube, [], None, f"FATAL:{e}"
    except Exception as e:
        return cube, [], None, str(e)


def main():
    ap = argparse.ArgumentParser(
        description="Scan every pixel of each cube for specular viewing "
                    "geometry (localemission == localincidence, spacecraft "
                    "opposite the sun in azimuth) and write matches to a CSV.")
    ap.add_argument("--ei-only", action="store_true",
                    help="match on |localemission - localincidence| <= "
                         "--tol-angle ALONE, ignoring the spacecraft-sun "
                         "azimuth condition entirely (--tol-az is then "
                         "unused). Finds the whole equal-angle locus, not "
                         "just mirror geometry; typically used with a wider "
                         "--tol-angle, e.g. --ei-only --tol-angle 5. "
                         "spec_score becomes ei_diff/tol_angle.")
    ap.add_argument("--tol-angle", type=float, default=1.0,
                    help="max |localemission - localincidence| in degrees "
                         "(default: 1.0)")
    ap.add_argument("--tol-az", type=float, default=5.0,
                    help="max deviation of the spacecraft-sun relative "
                         "azimuth from 180 degrees (default: 5.0)")
    ap.add_argument("--min-angle", type=float, default=1.0,
                    help="exclude pixels whose localemission or "
                         "localincidence is below this many degrees, where "
                         "azimuths are ill-defined and the specular test is "
                         "trivially/spuriously satisfied (default: 1.0)")
    ap.add_argument("--max-per-cube", type=int, default=None,
                    help="keep only the N most-specular pixels per cube "
                         "(lowest spec_score). Default: keep all matches.")
    ap.add_argument("--geom-dir", "--backplane-dir", default=None,
                    dest="geom_dir",
                    help="directory containing the .geom.cub cubes "
                         "(default: same directory as each input cube)")
    ap.add_argument("--jobs", "-j", type=int, default=1,
                    help="scan N cubes in parallel with multiprocessing "
                         "(default: 1)")
    ap.add_argument("--band", "-b", type=int, default=None, choices=[1, 2, 4],
                    help="only process cubes for this MARCI band, i.e. files "
                         "whose basename contains '.bandNNNN.' for the given "
                         "band (both even and odd framelets). Default: no "
                         "band filtering.")
    ap.add_argument("--lat-range", default=None, metavar="MIN,MAX",
                    help="only keep pixels whose latitude backplane is in "
                         "[MIN, MAX] degrees (southern latitudes NEGATIVE, "
                         "e.g. -90,-80). Replaces camtrim cropping - exact "
                         "same camera-model latitudes, zero ISIS cost. "
                         "NOTE: negative values need the = form, "
                         "--lat-range=-90,-80 (a space-separated value "
                         "starting with '-' is rejected by the parser).")
    ap.add_argument("--ls", default=None, metavar="MIN,MAX",
                    help="only process cubes whose Ls (parsed from the "
                         "filename) is in [MIN, MAX] inclusive - use the same "
                         "range you gave make_specular_geom.sh -L. "
                         "Incompatible with --no-ls.")
    ap.add_argument("--no-ls", action="store_true",
                    help="do not parse Ls from filenames (ls column left "
                         "empty). Without this flag an unparseable filename "
                         "aborts the whole run.")
    ap.add_argument("--no-header", action="store_true",
                    help="do not write the CSV header line")
    ap.add_argument("-o", "--output", help="write results to this CSV file "
                                           "(default: stdout)")
    ap.add_argument("cubes", nargs="*",
                    help="ORIGINAL cube files (not the .geom.cub cubes). "
                         "Default: every *.cub in the current directory, "
                         "excluding geom cubes.")
    args = ap.parse_args()

    if gdal is None:
        err("osgeo.gdal is not importable - install GDAL Python bindings.")
        sys.exit(1)
    if args.tol_angle <= 0 or args.tol_az <= 0:
        err("--tol-angle and --tol-az must be > 0")
        sys.exit(1)
    if args.min_angle < 0:
        err("--min-angle must be >= 0")
        sys.exit(1)
    if args.max_per_cube is not None and args.max_per_cube < 1:
        err("--max-per-cube must be >= 1")
        sys.exit(1)

    if args.cubes:
        cubes = args.cubes
    else:
        geom_suffix = f".{GEOM_SUFFIX}.cub"
        cubes = sorted(c for c in glob.glob("*.cub")
                       if not c.endswith(geom_suffix))
        if not cubes:
            err("no cubes given and no original *.cub found here")
            sys.exit(1)

    lat_range = None
    if args.lat_range is not None:
        try:
            lo, hi = (float(x) for x in args.lat_range.split(","))
        except ValueError:
            err(f"--lat-range must be MIN,MAX (got '{args.lat_range}')")
            sys.exit(1)
        if lo > hi:
            err(f"--lat-range min ({lo}) must be <= max ({hi}); remember "
                f"southern latitudes are negative (south pole: -90,-80)")
            sys.exit(1)
        lat_range = (lo, hi)

    ls_range = None
    if args.ls is not None:
        if args.no_ls:
            err("--ls and --no-ls are incompatible")
            sys.exit(1)
        try:
            lo, hi = (float(x) for x in args.ls.split(","))
        except ValueError:
            err(f"--ls must be MIN,MAX (got '{args.ls}')")
            sys.exit(1)
        if lo > hi:
            err(f"--ls min ({lo}) must be <= max ({hi})")
            sys.exit(1)
        ls_range = (lo, hi)

    if args.band is not None:
        tag = f".band{args.band:04d}."
        before = len(cubes)
        cubes = [c for c in cubes if tag in os.path.basename(c)]
        if not cubes:
            err(f"--band {args.band}: no input cube matched '{tag}' "
                f"({before} cube(s) before filtering)")
            sys.exit(1)

    missing = [c for c in cubes if not os.path.isfile(c)]
    if missing:
        for c in missing:
            err(f"cube not found: {c}")
        sys.exit(1)

    # Hard-fail EARLY on unparseable Ls filenames, before any heavy work.
    if not args.no_ls:
        for c in cubes:
            try:
                parse_ls(c)
            except ValueError as e:
                err(f"error: {e}")
                sys.exit(1)

    if ls_range is not None:
        before = len(cubes)
        cubes = [c for c in cubes
                 if ls_range[0] <= parse_ls(c) <= ls_range[1]]
        if not cubes:
            err(f"--ls {args.ls}: no input cube in range "
                f"({before} cube(s) before filtering)")
            sys.exit(1)

    cfg = {"geom_dir": args.geom_dir, "tol_angle": args.tol_angle,
           "tol_az": args.tol_az, "min_angle": args.min_angle,
           "max_per_cube": args.max_per_cube, "want_ls": not args.no_ls,
           "lat_range": lat_range, "ei_only": args.ei_only}
    tasks = [(c, cfg) for c in cubes]

    out = open(args.output, "w", newline="") if args.output else sys.stdout
    n_rows = 0
    n_hit_cubes = 0
    try:
        w = csv.DictWriter(out, fieldnames=COLUMNS, extrasaction="ignore")
        if not args.no_header:
            w.writeheader()

        def handle(result):
            nonlocal n_rows, n_hit_cubes
            cube, rows, reason, error = result
            bn = os.path.basename(cube)
            if error is not None:
                if error.startswith("FATAL:"):
                    err(f"error: {error[6:]}")
                    sys.exit(1)
                err(f"skip (error: {error}): {bn}")
                return
            if reason:
                err(f"skip ({reason}): {bn}")
                return
            if not rows:
                err(f"no specular pixels: {bn}")
                return
            n_hit_cubes += 1
            for r in rows:
                rounded = {k: (round(v, 3) if isinstance(v, float)
                               and not math.isnan(v) else v)
                           for k, v in r.items()}
                w.writerow(rounded)
                n_rows += 1

        if args.jobs > 1:
            import multiprocessing as mp
            with mp.Pool(args.jobs) as pool:
                for result in pool.imap_unordered(_worker, tasks):
                    handle(result)
        else:
            for t in tasks:
                handle(_worker(t))
    finally:
        if args.output:
            out.close()

    err(f"{n_rows} specular pixel(s) from {n_hit_cubes}/{len(cubes)} cube(s).")
    if n_rows == 0:
        sys.exit(2)


if __name__ == "__main__":
    main()