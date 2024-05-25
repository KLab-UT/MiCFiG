#!/bin/bash

destination="hold"

mkdir -p "$destination"

while IFS= read -r file; do
    if [ -f "$file" ]; then
        mv "$file" "$destination"
        echo "Moved $file to $destination"
    else
        echo "File $file does not exist"
    fi
done < keep_files.txt
