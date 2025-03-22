#!/bin/bash

# Find all files tracked by Git LFS (e.g., .dat files)
files=$(git lfs ls-files | awk '{print $1}')

# Loop through each file and push it individually
for file in $files; do
    # Add the file to git
    git add "$file"
    
    # Commit the file
    git commit -m "Add $file to LFS"
    
    # Push the file
    git push origin main
    
    # Optional: Print status after each push
    echo "Pushed $file"
done
