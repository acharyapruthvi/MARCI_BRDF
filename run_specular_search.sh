#!/usr/bin/env bash
#
# run_specular_search.sh - full specular-reflection search pipeline in one
#   command: generate phocube geometry cubes, scan them for specular pixels
#   with find_specular.py, append the hits to a CSV, and DELETE the geometry
#   cubes - all in bounded BATCHES so disk space is never a limiting factor.
#
#   (This supersedes make_specular_geom.sh; the geom cubes are now temporary
#   intermediates, not products.)
#
# PIPELINE, per batch of cubes:
#   1. generate : ONE worker per cube runs spice check -> spiceinit (if
#                 needed) -> phocube; all workers share one `xargs -P` pool
#                 (no phase barriers within a batch). There is NO camtrim
#                 stage: latitude cropping is done in find_specular.py from
#                 the latitude backplane (same camera-model values, one whole
#                 ISIS pass per cube cheaper, and no camtrim OOM risk).
#   2. scan     : ONE find_specular.py call scans the whole batch with its
#                 own multiprocessing pool (--jobs = -j) and appends rows
#                 (no header) to the output CSV.
#   3. clean    : the batch's geom cubes are deleted (unless -k).
#   Peak temp disk = one batch of geom cubes (~ -B cubes x 8 float32 bands at
#   full cube dimensions). Cleanup also fires on error/interrupt via trap.
#
# Cube selection (pick ONE):
#   -t PATH      test mode: every *.cub in the single directory PATH
#   (default)    bulk mode: every *.cub in every immediate subdirectory of
#                -r ROOT, i.e. ROOT/*/*.cub
#   cube.cub ... explicit cube paths
#
# Band selection:
#   -b N         only process cubes for MARCI band N (1, 2, or 4), i.e. files
#                whose basename contains ".band000N." - both the .even. and
#                .odd. framelet halves are kept. REQUIRED unless explicit
#                cube paths are given (then it further filters those paths).
#
# Ls (solar longitude) filtering:
#   -L MIN,MAX   only process cubes whose Ls (parsed from the filename: 3rd
#                underscore field, LAST digit = decimal place, e.g.
#                B04_011323_1800_... -> 180.0) is in [MIN, MAX] inclusive.
#                An unparseable filename aborts the WHOLE RUN.
#
# Usage:
#   run_specular_search.sh -b BAND -o FILE [options] [cube.cub ...]
#
# Options:
#   -b N        MARCI band to process: 1, 2, or 4 (see above)
#   -o FILE     output CSV (default: stdout). Header written once up front;
#               each batch appends rows.
#   -r ROOT     bulk-mode root dir (default: /data/MARCI_IMAGES/DONE)
#   -t PATH     test mode: process just this one directory's *.cub files
#   -L MIN,MAX  Ls range filter (see above)
#   -d DIR      scratch dir for geom cubes (default: an auto-created
#               .specular_geom_$$ next to this script, removed at exit)
#   -j N        max parallel workers, for BOTH cube generation and the
#               find_specular.py scan (default: all cores via `nproc`).
#               NOTE: this is also the memory throttle for the scan stage
#               (each scan worker holds several full-frame float32 arrays).
#   -B N        cubes per batch (default: 4 x -j). Larger = better core
#               saturation at the generate->scan boundary, more peak disk.
#   -k          keep the geom cubes instead of deleting them
#   -f          force: regenerate geom cubes even if they already exist
#   --crop-lat MIN,MAX  only keep hits whose latitude is in this band
#               (applied in find_specular.py from the latitude backplane; no
#               camtrim is run). REMEMBER: southern latitudes are NEGATIVE
#               (south pole: --crop-lat -90,-80).
#   --tol-angle DEG    passed to find_specular.py (default: its default, 1.0)
#   --tol-az DEG       passed to find_specular.py (default: its default, 5.0)
#   --min-angle DEG    passed to find_specular.py (default: its default, 1.0)
#   --max-per-cube N   passed to find_specular.py (default: keep all hits)
#   -w          run spiceinit via the SPICE Web Service if a cube needs it
#   -S          skip the spiceinit check entirely
#   --keep-prt  do NOT suppress ISIS session logs (print.prt)
#   --py PATH      python executable for find_specular.py (default: python3)
#   --script PATH  path to find_specular.py (default: alongside this script)
#   -h          show this help and exit
#
# IMPORTANT: place all -OPTION flags BEFORE any explicit cube paths.
#
# Examples:
#   run_specular_search.sh -b 1 -L 220,222 --crop-lat -90,-80 -j 24 -o specular_band1.csv
#   run_specular_search.sh -b 4 -t /data/MARCI_IMAGES/DONE/B04_011323_1800_MA_00N177W -o test.csv
#   run_specular_search.sh -b 2 --max-per-cube 50 --tol-az 2 -j 24 -o specular_band2.csv

set -eu

prog=${0##*/}
die()  { printf '%s: error: %s\n' "$prog" "$*" >&2; exit 1; }
note() { printf '%s [%s] %s\n' "$prog" "$(date +%H:%M:%S)" "$*" >&2; }
usage() { sed -n '3,86p' "$0" | sed 's/^#\{0,1\} \{0,1\}//'; }

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Must match COLUMNS in find_specular.py (the scan appends with --no-header).
header="filename,ls,line,sample,dn,lat,lon,phase,localemission,localincidence,ei_diff,spacecraftazimuth,sunazimuth,relaz,offnadir,spec_score"

# ---------------------------------------------------------------- options ----
band=""; outfile=""; geom_dir=""; web=false; skip_spice=false; force=false
keep=false; keep_prt=false
root="/data/MARCI_IMAGES/DONE"; test_dir=""; ls_range=""; njobs=""; batch_cubes_n=""
crop_minlat=""; crop_maxlat=""
tol_angle=""; tol_az=""; min_angle=""; max_per_cube=""
pybin="python3"; pyscript="$script_dir/find_specular.py"

# Old per-band backplane suffixes; recognized only so leftovers are never
# treated as source cubes.
names=(lat lon phase localemission localincidence spacecraftazimuth sunazimuth offnadir)

args=()
while (( $# )); do
  case "$1" in
    -b) band=$2; shift 2 ;;
    -o) outfile=$2; shift 2 ;;
    -r) root=$2; shift 2 ;;
    -t) test_dir=$2; shift 2 ;;
    -L) ls_range=$2; shift 2 ;;
    -d) geom_dir=$2; shift 2 ;;
    -j) njobs=$2; shift 2 ;;
    -B) batch_cubes_n=$2; shift 2 ;;
    -k) keep=true; shift ;;
    -f) force=true; shift ;;
    --crop-lat) IFS=',' read -r crop_minlat crop_maxlat _ <<<"$2"; shift 2 ;;
    --tol-angle) tol_angle=$2; shift 2 ;;
    --tol-az) tol_az=$2; shift 2 ;;
    --min-angle) min_angle=$2; shift 2 ;;
    --max-per-cube) max_per_cube=$2; shift 2 ;;
    -w) web=true; shift ;;
    -S) skip_spice=true; shift ;;
    --keep-prt) keep_prt=true; shift ;;
    --py) pybin=$2; shift 2 ;;
    --script) pyscript=$2; shift 2 ;;
    -h) usage; exit 0 ;;
    --) shift; args+=("$@"); break ;;
    -*) die "unknown option $1 (use -h for help)" ;;
    *)  args+=("$1"); shift ;;
  esac
done
set -- "${args[@]}"

is_num() { [[ $1 =~ ^[+-]?([0-9]+\.?[0-9]*)|(\.[0-9]+)([eE][+-]?[0-9]+)?$ ]]; }

# --- band: required (unless explicit cubes given), must be 1, 2, or 4 -------
band_tag=""
if [[ -n $band ]]; then
  case "$band" in
    1|2|4) band_tag=$(printf '.band%04d.' "$band") ;;
    *) die "-b must be 1, 2, or 4 (got '$band')" ;;
  esac
elif (( $# == 0 )); then
  die "-b BAND is required (1, 2, or 4) unless explicit cube paths are given"
fi
[[ -n $band_tag ]] && note "band filter: basenames containing '$band_tag' (even + odd framelets)"

do_crop=false
if [[ -n $crop_minlat || -n $crop_maxlat ]]; then
  [[ -n $crop_minlat && -n $crop_maxlat ]] || die "--crop-lat needs MIN,MAX"
  is_num "$crop_minlat" || die "--crop-lat min not numeric: '$crop_minlat'"
  is_num "$crop_maxlat" || die "--crop-lat max not numeric: '$crop_maxlat'"
  awk -v a="$crop_minlat" -v b="$crop_maxlat" 'BEGIN{ exit !(a<=b) }' \
    || die "--crop-lat min ($crop_minlat) must be <= max ($crop_maxlat)"
  do_crop=true
fi

for _v in tol_angle tol_az min_angle; do
  [[ -n ${!_v} ]] && { is_num "${!_v}" || die "--${_v//_/-} must be numeric: '${!_v}'"; }
done
if [[ -n $max_per_cube ]]; then
  [[ $max_per_cube =~ ^[0-9]+$ && $max_per_cube -ge 1 ]] || die "--max-per-cube must be a positive integer"
fi

command -v phocube >/dev/null 2>&1 || \
  die "phocube not on PATH - activate your ISIS environment first (e.g. 'conda activate isis')"
command -v xargs >/dev/null 2>&1 || die "xargs not on PATH (needed for parallelism)"
[[ -f $pyscript ]] || die "find_specular.py not found at: $pyscript (use --script to point at it)"
command -v "$pybin" >/dev/null 2>&1 || die "python executable not found: $pybin"

if [[ -n $njobs ]]; then
  [[ $njobs =~ ^[0-9]+$ && $njobs -ge 1 ]] || die "-j must be a positive integer"
else
  njobs=$(nproc 2>/dev/null || echo 1)
fi
note "parallelism: up to $njobs worker(s) (generation and scan)"

if [[ -n $batch_cubes_n ]]; then
  [[ $batch_cubes_n =~ ^[0-9]+$ && $batch_cubes_n -ge 1 ]] || die "-B must be a positive integer"
else
  batch_cubes_n=$(( njobs * 4 ))
  (( batch_cubes_n < 1 )) && batch_cubes_n=1
fi
note "batch size: $batch_cubes_n cube(s) per generate->scan->delete cycle"

ls_min=""; ls_max=""; ls_min_t=""; ls_max_t=""
if [[ -n $ls_range ]]; then
  read -r ls_min ls_max _ <<<"${ls_range//,/ }"
  [[ -n ${ls_min:-} && -n ${ls_max:-} ]] || die "could not parse Ls range '$ls_range' (need MIN,MAX)"
  is_num "$ls_min" || die "Ls range min not numeric: '$ls_min'"
  is_num "$ls_max" || die "Ls range max not numeric: '$ls_max'"
  awk -v a="$ls_min" -v b="$ls_max" 'BEGIN{ exit !(a<=b) }' \
    || die "Ls range min ($ls_min) must be <= max ($ls_max)"
  # Integer tenths for fork-free per-cube comparison (Ls in a filename is
  # already tenths: field "2205" -> 220.5 deg -> 2205).
  ls_min_t=$(awk -v x="$ls_min" 'BEGIN{ printf "%d", int(x*10 + 0.5) }')
  ls_max_t=$(awk -v x="$ls_max" 'BEGIN{ printf "%d", int(x*10 + 0.5) }')
fi

# try_parse_ls_tenths NAME : set global ls_tenths from NAME's basename (3rd
# underscore field, last digit = decimal place) in integer tenths. Returns 1
# (without dying) if the basename doesn't follow the scheme. Pure bash - no
# forks - because this runs once per file/dir in the archive and fork cost
# dominates discovery time, especially on WSL.
ls_tenths=""
try_parse_ls_tenths() {
  local f=${1##*/} rest field
  rest=${f#*_};    [[ $rest != "$f" ]] || return 1
  rest=${rest#*_}; field=${rest%%_*}
  [[ $field =~ ^[0-9]{2,}$ ]] || return 1
  ls_tenths=$((10#$field))
}

# parse_ls_tenths FILE : same, but a malformed CUBE filename is a hard error
# (dies in the MAIN shell, not a subshell, so the abort actually happens).
parse_ls_tenths() {
  try_parse_ls_tenths "$1" || \
    die "cannot parse Ls from filename '${1##*/}' (expected a numeric 3rd underscore field, e.g. B04_011323_1800_...)"
}

is_backplane() {
  local f=${1##*/} n
  [[ $f == *.geom.cub || $f == *.trim.cub ]] && return 0
  for n in "${names[@]}"; do
    [[ $f == *.$n.cub ]] && return 0
  done
  return 1
}

in_ls_range() {
  [[ -n $ls_range ]] || return 0
  parse_ls_tenths "$1"
  (( ls_tenths >= ls_min_t && ls_tenths <= ls_max_t ))
}

in_band() {
  [[ -n $band_tag ]] || return 0
  [[ ${1##*/} == *"$band_tag"* ]]
}

# ------------------------------------------------------ discover the cubes --
declare -a all_cubes=()
mode=""
n_seen=0

# progress_tick: called once per file examined during discovery.
progress_tick() {
  n_seen=$((n_seen + 1))
  if (( n_seen % 1000 == 0 )); then
    note "  discovery: $n_seen file(s) examined, ${#all_cubes[@]} selected so far..."
  fi
}

if (( $# )); then
  mode="explicit"
  for a in "$@"; do
    [[ $a == -* ]] && die "'$a' looks like an option but came after a cube path; put options first (see -h)"
    [[ -f $a ]] || die "cube not found: $a"
    progress_tick
    is_backplane "$a" && continue
    in_band "$a" || continue
    in_ls_range "$a" || continue
    all_cubes+=("$a")
  done
elif [[ -n $test_dir ]]; then
  mode="test"
  [[ -d $test_dir ]] || die "test directory not found: $test_dir"
  note "discovering cubes in $test_dir ..."
  shopt -s nullglob
  for c in "$test_dir"/*.cub; do
    progress_tick
    is_backplane "$c" && continue
    in_band "$c" || continue
    in_ls_range "$c" || continue
    all_cubes+=("$c")
  done
  shopt -u nullglob
else
  mode="bulk"
  [[ -d $root ]] || die "root directory not found: $root"
  note "discovering cubes under $root/*/ ..."
  shopt -s nullglob

  # LEVEL 1 - directories. Observation dir names carry the SAME Ls field as
  # their cubes (B04_011323_1800_MA_00N177W), so with -L we can reject a
  # whole directory from its NAME ALONE, never listing its contents. For a
  # narrow Ls window this skips ~99% of directories, which is where nearly
  # all the discovery time goes on large archives. A dir whose name doesn't
  # parse is NOT skipped - we descend and let the per-file checks (which DO
  # hard-fail on malformed cube names) decide.
  keep_dirs=()
  n_dirs=0
  for d in "$root"/*/; do
    d=${d%/}
    n_dirs=$((n_dirs + 1))
    if [[ -n $ls_range ]] && try_parse_ls_tenths "$d"; then
      (( ls_tenths >= ls_min_t && ls_tenths <= ls_max_t )) || continue
    fi
    keep_dirs+=("$d")
  done
  note "  directory pre-filter: ${#keep_dirs[@]} of $n_dirs dir(s) match the Ls range (others never listed)"

  # LEVEL 2 - files, only inside surviving dirs, with the band tag baked
  # into the glob so non-matching files never reach the loop. Per-file
  # checks still run as the authoritative filter (cheap at this point).
  band_glob="*.cub"
  [[ -n $band_tag ]] && band_glob="*${band_tag}*.cub"
  for d in "${keep_dirs[@]}"; do
    for c in "$d"/$band_glob; do
      progress_tick
      is_backplane "$c" && continue
      in_band "$c" || continue
      in_ls_range "$c" || continue
      all_cubes+=("$c")
    done
  done
  shopt -u nullglob
fi
(( ${#all_cubes[@]} )) || die "no cubes selected ($mode mode; examined $n_seen file(s); check -b/-L filters and -r/-t path)"
note "$mode mode: ${#all_cubes[@]} cube(s) selected (of $n_seen examined)"

# --- scratch dir for the temporary geom cubes ---------------------------------
if [[ -n $geom_dir ]]; then
  geom_dir_auto=false
else
  geom_dir="$script_dir/.specular_geom_$$"
  geom_dir_auto=true
fi
mkdir -p "$geom_dir"

# --- optional IsisPreferences that silences per-process session logs --------
pref_file=""
if ! $keep_prt; then
  pref_file=$(mktemp "${TMPDIR:-/tmp}/isisprefs.XXXXXX")
  printf 'Group = SessionLog\n  TerminalOutput = Off\n  FileOutput      = Off\nEndGroup\n' > "$pref_file"
fi

# --- cleanup: current batch's geom cubes + scratch dir + prefs, on any exit --
current_geoms=()
cleanup_batch() {
  $keep && { current_geoms=(); return 0; }
  (( ${#current_geoms[@]} == 0 )) && return 0
  rm -f "${current_geoms[@]}"
  current_geoms=()
}
cleanup_all() {
  cleanup_batch
  [[ -n $pref_file ]] && rm -f "$pref_file" 2>/dev/null || true
  if ! $keep && [[ ${geom_dir_auto:-false} == true && -d $geom_dir ]]; then
    rmdir "$geom_dir" 2>/dev/null || true
  fi
}
trap cleanup_all EXIT INT TERM

# ------------------------------------------------------------- cube worker --
# ONE worker does the ENTIRE generation chain for ONE cube:
#   getkey spice check -> spiceinit (if needed) -> phocube.
# Config from exported RSS_* variables. Arg: the source cube path.
process_one() {
  local cub=$1 bn out rc
  bn=${cub##*/}; bn=${bn%.cub}
  out="$RSS_GEOM_DIR/$bn.geom.cub"

  local pref=()
  [[ -n $RSS_PREF ]] && pref=(-preference="$RSS_PREF")

  if [[ -f $out ]]; then
    if [[ $RSS_FORCE == true ]]; then rm -f "$out"
    else printf 'SKIP_EXISTS\t%s\n' "${cub##*/}" >&2; return 0; fi
  fi

  # --- spice ---
  if [[ $RSS_SKIP_SPICE != true ]]; then
    if ! getkey "${pref[@]}" from="$cub" grpname=Kernels keyword=InstrumentPointing 2>/dev/null | grep -q .; then
      if [[ $RSS_WEB == true ]]; then
        spiceinit "${pref[@]}" from="$cub" web=true >/dev/null 2>&1 || rc=$?
      else
        spiceinit "${pref[@]}" from="$cub" >/dev/null 2>&1 || rc=$?
      fi
      if [[ -n ${rc:-} ]]; then
        printf 'SPICE_FAIL\t%s\n' "$cub" >&2
        return 1
      fi
    fi
  fi

  # --- phocube: ONE multi-band geom cube; bands identified BY DESCRIPTION
  # downstream (phocube emits its own canonical order, not this flag order).
  # No camtrim: latitude cropping happens in find_specular.py via the
  # latitude backplane - same camera-model values, zero extra ISIS cost. ---
  if ! phocube "${pref[@]}" from="$cub" to="$out" \
        LATITUDE=true LONGITUDE=true PHASE=true \
        LOCALEMISSION=true LOCALINCIDENCE=true \
        SPACECRAFTAZIMUTH=true SUNAZIMUTH=true OFFNADIRANGLE=true \
        EMISSION=false INCIDENCE=false DN=false \
        >/dev/null 2>&1; then
    printf 'PHOCUBE_FAIL\t%s\n' "$cub" >&2
    rm -f "$out"
    return 1
  fi

  printf 'GEOM_OK\t%s\n' "${cub##*/}" >&2
  return 0
}
export -f process_one
export RSS_GEOM_DIR="$geom_dir" RSS_FORCE="$force" RSS_SKIP_SPICE="$skip_spice"
export RSS_WEB="$web" RSS_PREF="$pref_file"

# ----------------------------------------------------------------- header ----
if [[ -n $outfile ]]; then
  mkdir -p "$(dirname "$outfile")" 2>/dev/null || true
  echo "$header" > "$outfile"
  note "output -> $outfile"
else
  echo "$header"
  note "output -> stdout"
fi

# find_specular.py flags shared by every batch. Band/Ls filtering already
# happened in bash, so neither --band nor --ls is passed; Ls PARSING stays on
# so the ls column is populated (bash validated the filenames when -L was
# given; MARCI names are uniform so mid-run parse failures abort loudly,
# which is the intended hard-fail behavior).
py_flags=(--geom-dir "$geom_dir" --jobs "$njobs" --no-header)
[[ -n $tol_angle ]]    && py_flags+=(--tol-angle "$tol_angle")
[[ -n $tol_az ]]       && py_flags+=(--tol-az "$tol_az")
[[ -n $min_angle ]]    && py_flags+=(--min-angle "$min_angle")
[[ -n $max_per_cube ]] && py_flags+=(--max-per-cube "$max_per_cube")
[[ -n $crop_minlat ]]  && py_flags+=("--lat-range=$crop_minlat,$crop_maxlat")

# =================================================================== run =====
ncubes=${#all_cubes[@]}
nbatches=$(( (ncubes + batch_cubes_n - 1) / batch_cubes_n ))
total_rows=0; total_gen_fail=0; batch_i=0
start=0
while (( start < ncubes )); do
  end=$(( start + batch_cubes_n ))
  (( end > ncubes )) && end=$ncubes
  batch=( "${all_cubes[@]:start:end-start}" )
  start=$end
  batch_i=$((batch_i + 1))
  note "batch $batch_i/$nbatches: ${#batch[@]} cube(s)"

  # --- 1. generate: one chained worker per cube, single pool ---
  note "  generating ${#batch[@]} geom cube(s) (up to $njobs in parallel; GEOM_OK lines stream as cubes finish)"
  printf '%s\n' "${batch[@]}" \
    | xargs -d '\n' -P "$njobs" -I{} bash -c 'process_one "$@"' _ {} || \
    note "  some cube(s) failed generation (see *_FAIL lines above)"

  # register the batch's geom cubes for cleanup; scan only the ones that exist
  scan_cubes=()
  current_geoms=()
  for cub in "${batch[@]}"; do
    bn=${cub##*/}; bn=${bn%.cub}
    g="$geom_dir/$bn.geom.cub"
    if [[ -f $g ]]; then
      scan_cubes+=("$cub")
      current_geoms+=("$g")
    else
      total_gen_fail=$((total_gen_fail + 1))
    fi
  done

  if (( ${#scan_cubes[@]} == 0 )); then
    note "  no geom cubes generated in this batch - nothing to scan"
    cleanup_batch
    continue
  fi

  # --- 2. scan: one find_specular.py call for the whole batch ---
  # exit 0 = hits written, exit 2 = zero hits (fine); anything else = warn.
  note "  scanning ${#scan_cubes[@]} cube(s) with find_specular.py (--jobs $njobs)"
  rows_before=0
  [[ -n $outfile ]] && rows_before=$(wc -l < "$outfile")
  scan_rc=0
  if [[ -n $outfile ]]; then
    "$pybin" "$pyscript" "${py_flags[@]}" "${scan_cubes[@]}" >> "$outfile" || scan_rc=$?
  else
    "$pybin" "$pyscript" "${py_flags[@]}" "${scan_cubes[@]}" || scan_rc=$?
  fi
  if (( scan_rc != 0 && scan_rc != 2 )); then
    note "  find_specular.py exited with status $scan_rc for this batch"
  fi
  if [[ -n $outfile ]]; then
    rows_after=$(wc -l < "$outfile")
    total_rows=$(( total_rows + rows_after - rows_before ))
  fi

  # --- 3. clean: delete this batch's geom cubes ---
  cleanup_batch
done

if [[ -n $outfile ]]; then
  note "done: $ncubes cube(s) selected, $total_gen_fail failed generation, $total_rows specular row(s) -> $outfile"
  (( total_rows > 0 )) || exit 2
else
  note "done: $ncubes cube(s) selected, $total_gen_fail failed generation (rows went to stdout)"
fi