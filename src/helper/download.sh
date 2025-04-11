#!/bin/bash

root_dir=$1
url_file=$2

# mkdir $root_dir
while IFS= read -r line
do
    filepath=$(echo $line|cut -d" " -f1);
    dir=$(dirname $filepath)
    mkdir -p $root_dir/$dir;
    link=$(printf "$line"|awk -F'\t' '{print $2}')
    wget -c --wait=30 --no-verbose -t 2 --http-user=anonymous -O $root_dir/$filepath "$link";
    # -a $root_dir/pride_archive_rewrite_download_wget_log to keep log
    
done < $url_file
