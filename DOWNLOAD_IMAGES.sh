#!/bin/bash
################################
# This script automatically downloads MARCI cub files from the PDS website. 
# This script only downloads up to Ls 170
################################

# Specify the URL to search for text files
url="https://pds-imaging.jpl.nasa.gov/data/mro/mars_reconnaissance_orbiter/marci/"

# Create a directory to store the processed images
mkdir -p NEED_TO_PROCESS

# Prompt the user to enter the orbit number range
read -p "Enter the starting orbit number (e.g., 23): " orbit_number_start
read -p "Enter the ending orbit number (e.g., 30): " orbit_number_end

# Function to reset lower_limit
reset_lower_limit() {
    lower_limit=1700
    echo "RESET LOWER_LIMIT"
}

# Download a single file
download_file() {
    local orbit_number=$1
    local file="mrom_$(printf "%04d" $orbit_number)_md5.txt"
    echo "Downloading $file..."
    if wget -q "$url$file" -P PDS_TXT/; then
        echo "Download successful"
    else
        echo "Failed to download $file. Try another orbit number."
        exit 1
    fi
}

# Process a single file
process_file() {
    local file_path=$1
    local lower_limit=$2

    echo "Processing $file_path..."
    echo "------------------"
    while IFS= read -r line; do
        if [[ $line == *.IMG ]] && ( [[ $line == *"_MA_"* ]] || [[ $line == *"_MC_"* ]] ); then
            filename=$(echo "$line" | awk '{print $2}')
            image_name="${filename#*mrom}"
            date_part="${image_name%_*}"
            value="${date_part: -7:4}"
            if [[ $lower_limit -gt 1700 ]]; then
                lower_limit=1300
            fi

            echo "CHECKING $value IS WITHIN THE $lower_limit and 1700 BOUNDS"
            if [[ $lower_limit -lt 1700 ]] && [[ $value -ge $lower_limit ]] && [[ $value -le 1700 ]]; then
                img_name=$(basename "$filename" .IMG)
                done_folder_path="DONE/$img_name"
                need_to_process_folder_path="NEED_TO_PROCESS/$img_name"

                # Check if the folder already exists in DONE directory
                if [ -d "$done_folder_path" ]; then
                    echo "Folder $done_folder_path already exists. Skipping image."
                else
                    echo "CURRENT LOWER_LIMIT IS: $lower_limit"
                    echo "CURRENTLY DOWNLOADING $image_name"
                    echo "DOWNLOADING #$count IMAGES"
                    lower_limit=$((lower_limit + 5))
                    echo $lower_limit

                    # Check if the folder already exists in NEED_TO_PROCESS directory
                    if [ -d "$need_to_process_folder_path" ]; then
                        echo "Folder $need_to_process_folder_path already exists. Skipping download."
                    else
                        mkdir -p "$need_to_process_folder_path"
                        wget -q "https://pds-imaging.jpl.nasa.gov/data/mro/mars_reconnaissance_orbiter/marci/mrom$image_name" -P "$need_to_process_folder_path/"
                    fi

                    ((count++))
                    if ((count == 150)); then
                        echo "DONE DOWNLOADING $count IMAGES"
                        exit 0
                    fi
                fi
            else
                echo "RESETTING"
                lower_limit=1300
            fi
        fi
    done < "$file_path"
}

# Create a directory to store the downloaded text files
mkdir -p PDS_TXT

# Download the text files between mrom_0023 and mrom_0030
for ((i=orbit_number_start; i<=orbit_number_end; i++)); do
    download_file "$i" &
done

# Wait for all downloads to finish
wait

# Create an array of the downloaded text files
file_list=(PDS_TXT/*)

lower_limit=1300
count=0 # counter to track how many times a number has been echoed

# Process each file in parallel
for file_path in "${file_list[@]}"; do
    process_file "$file_path" "$lower_limit" &
done

# Wait for all processing to finish
wait

rm -rf PDS_TXT/*
