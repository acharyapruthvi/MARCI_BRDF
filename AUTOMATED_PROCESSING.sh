#!/bin/bash
## This bash script is a wrapper around the ISIS_AUTOMATED.sh script which processes MARCI images from specific Ls date. 
eval "$(conda shell.bash hook)"
conda activate isis4

read -p "Enter the folder to search for .IMG files: " img_folder

total_images=$(find "$img_folder" -mindepth 1 -type d | wc -l) # Finds total number of images in an image folder
processed_images=0
last_completed_image=""

process_directory() {
    dir="$1"
    image_name=$(basename "$dir" .IMG)

    # If the file is already processed, then it ignores that cub file
    if [ -f "DONE/$image_name/$image_name.band0001.FINAL.cub" ]; then
        echo "Skipping $dir as $image_name.band0001.FINAL.cub already exists in DONE/$image_name/"
        return
    fi

    echo "Processing $dir..."

    # Processes the cub files and moves final processes PNG files to a directory
    mkdir -p DONE/$image_name
    echo $dir | bash -e ISIS_AUTOMATED.sh
    last_completed_image="$dir"
    mv $dir/*FINAL.cub DONE/$image_name/
    echo "Deleting things from $dir"
    find "$dir" -type f ! -name "*.IMG" -delete
    ls $dir > /dev/null
    echo "Finished processing $dir..."
    ((processed_images++))
}

export -f process_directory

find "$img_folder"/*/ -type d -print0 | xargs -0 -P "$(nproc)" -I {} bash -c 'process_directory "$@"' _ {}

echo "All images processed."
