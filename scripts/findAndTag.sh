#!/bin/bash

rawDataDir=/Users/andrewhsu/Projects/PREP-NURA/NURA/data/OptimizedCyclic-WSSL
rawPNGDir=/Users/andrewhsu/Projects/PREP-NURA/NURA/plots/OptimizedCyclic-WSSL

dataDir=$(echo "$rawDataDir" | tr -d '\r' | xargs)
pngDir=$(echo "$rawPNGDir" | tr -d '\r' | xargs)

tagName="Red"
searchValue="1"
colNum="5"

while IFS= read -r -d '' -u 3 file; do
    echo "Scanning File: $(basename "$file")" 
    
    if awk -F ',' -v col="$colNum" -v val="$searchValue" \
        'BEGIN {found=0}  {if($col == val) {found=1; exit}} END {if (found) exit 0; else exit 1}' "$file"; then
        echo "Tagging CSV: $file"
        tag --add "$tagName" "$file"
        
        fileName=$(basename "$file" .csv)
        parentDir=$(basename "$(dirname "$file")")
        pngPath="${pngDir}/${parentDir}/${fileName}.png"
        
        if ls "$pngPath" >/dev/null 2>&1; then
            echo "Tagging PNG: $pngPath"
            tag --add "$tagName" "$pngPath"
        else
            echo "$pngPath does not exist."
        fi
    fi
done 3< <(find "$dataDir" -name "*.csv" -print0)