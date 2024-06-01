#!bin/bash
################################
# This bash script processes the cub files using the method described in Robbins (2022)
################################
read -p "Enter the folder to search for .IMG files: " img_folder

# Check if the directory exists
image_file=$(find "$img_folder" -type f -name "*.IMG" -print -quit)

if [[ -n $image_file ]]; then
    image=$(basename "$image_file")
    image="${image%.*}"
else
    echo "No image file found in $directory"
fi
#Initial calibration and seperatig the files into different filters
marci2isis from=$img_folder/$image.IMG to=$img_folder/$image.cub flip=no 1> /dev/null
spiceinit from=$img_folder/$image.odd.cub cknadir=true 1> /dev/null
spiceinit from=$img_folder/$image.even.cub cknadir=true 1> /dev/null
marcical from=$img_folder/$image.odd.cub to=$img_folder/$image.odd.l1.cub 1> /dev/null
marcical from=$img_folder/$image.even.cub to=$img_folder/$image.even.l1.cub 1> /dev/null
explode from=$img_folder/$image.odd.l1.cub to=$img_folder/$image.odd.l1 1> /dev/null
explode from=$img_folder/$image.even.l1.cub to=$img_folder/$image.even.l1 1> /dev/null

## Processing for each filter
# Photometeric correction
for ((i=3; i<=3; i++))
do
  photomet from=$img_folder/$image.even.l1.band000$i.cub to=$img_folder/$image.even.l3.band000$i.cub frompvl=/home/pruthvi/anaconda3/envs/isis4/appdata/templates/photometry/Hapkehen_model.pvl 1> /dev/null
  photomet from=$img_folder/$image.odd.l1.band000$i.cub to=$img_folder/$image.odd.l3.band000$i.cub frompvl=/home/pruthvi/anaconda3/envs/isis4/appdata/templates/photometry/Hapkehen_model.pvl 1> /dev/null
done
#Figuring out the image size
echo $img_folder/$image.odd.l3.band0003.cub | python image_size.py > $img_folder/image_size.text
height=$(grep -oE 'Height: [0-9]+' $img_folder/image_size.text | awk '{print $2}')
width=$(grep -oE 'Width: [0-9]+' $img_folder/image_size.text | awk '{print $2}')
# New image size to crop the image and reduce the processing time. 
NSMAPLE=$(awk "BEGIN{printf(\"%.0f\", $width * 0.5)}")
SAMPLE=$(awk "BEGIN{printf(\"%.0f\", $width * 0.25)}")
# Cropping each filter
for ((i=3; i<=3; i++))
do
  crop from=$img_folder/$image.odd.l3.band000$i.cub to=$img_folder/$image.odd.band000$i.l3.crop.cub sample=$SAMPLE nsamples=$NSMAPLE 1> /dev/null
  crop from=$img_folder/$image.even.l3.band000$i.cub to=$img_folder/$image.even.band000$i.l3.crop.cub sample=$SAMPLE nsamples=$NSMAPLE 1> /dev/null
done
# Creating a projection map centers at 0E -90E up to -60N with a PPD of 24
maptemplate map=$img_folder/AA_equ_24ppd.map projection=polarstereographic clon=0 clat=-90 targopt=user targetname=Mars rngopt=user minlat=-90 maxlat=-60 minlon=-180 maxlon=180 maxlon=180 resopt=ppd resolution=24 1> /dev/null
# Polar stereographic projection
for ((i=3; i<=3; i++))
do
  #old code
  cam2map from=$img_folder/$image.odd.band000$i.l3.crop.cub to=$img_folder/$image.odd.band000$i.l2.cub map=$img_folder/AA_equ_24ppd.map pixres=map defaultrange=map 1> /dev/null
  cam2map from=$img_folder/$image.even.band000$i.l3.crop.cub to=$img_folder/$image.even.band000$i.l2.cub map=$img_folder/AA_equ_24ppd.map pixres=map defaultrange=map 1> /dev/null
  
  #Even files
  phocube from=$img_folder/$image.even.band000$i.l2.cub to=DONE/$image/$image.phase.even.band000$i.cub specialpixels=no emission=no incidence=no latitude=no longitude=no
  phocube from=$img_folder/$image.even.band000$i.l2.cub to=DONE/$image/$image.emission.even.band000$i.cub specialpixels=no phase=no emission=no incidence=no localemission=yes latitude=no longitude=no
  phocube from=$img_folder/$image.even.band000$i.l2.cub to=DONE/$image/$image.incidence.even.band000$i.cub specialpixels=no phase=no emission=no incidence=no localincidence=yes latitude=no longitude=no
  phocube from=$img_folder/$image.even.band000$i.l2.cub to=DONE/$image/$image.SA.even.band000$i.cub specialpixels=no phase=no emission=no incidence=no latitude=no longitude=no sunazimuth=yes
  phocube from=$img_folder/$image.even.band000$i.l2.cub to=DONE/$image/$image.SCA.even.band000$i.cub specialpixels=no phase=no emission=no incidence=no latitude=no longitude=no spacecraftazimuth=yes

done

for ((i=3; i<=3; i++))
do
  ls $img_folder/$image*.band000$i.l2.cub > $img_folder/merge.lis
  automos from=$img_folder/merge.lis mosaic=$img_folder/$image.band000$i.l2.cub priority=average 1> /dev/null
done
# Combining the difference frames into one final image
for ((i=3; i<=3; i++))
do
  explode from=$img_folder/$image.band000$i.l2.cub to=$img_folder/$image.band000$i.FINAL
  rm $img_folder/$image.band000$i.FINAL.band0002.cub
  mv $img_folder/$image.band000$i.FINAL.band0001.cub $img_folder/$image.band000$i.FINAL.cub
done

