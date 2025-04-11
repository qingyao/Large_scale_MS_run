#! /usr/bin/bash
year_month_folder=$1
taxID_name_file=../$year_month_folder/taxID_name.tsv
string_ver=$2
paxdb_ver=$3

while IFS=$'\t' read -r taxID _; 
do
    src=../rsc/$paxdb_ver/fasta/fasta.$string_ver.$taxID.fa.gz
    dest=../rsc/$paxdb_ver/fasta/fasta.$string_ver.$taxID.fa
    if [ ! -e $dest ]
    then
	wget -c -O $src "https://stringdb-downloads.org/download/protein.sequences.$string_ver/$taxID.protein.sequences.$string_ver.fa.gz"
	gunzip $src
	sleep 0.5
    fi
        
    src=../rsc/$paxdb_ver/links/$taxID.protein.links.$string_ver.txt.gz
    dest=../rsc/$paxdb_ver/links/$taxID.protein.links.$string_ver.txt
    filtered=../rsc/$paxdb_ver/links/$taxID.network_${string_ver}_900.txt
    if [ ! -e $dest ]    
    then
    wget -c -O $src "https://stringdb-downloads.org/download/protein.links.$string_ver/$taxID.protein.links.$string_ver.txt.gz"
    gunzip $src
    sleep 5
    fi

    if [ ! -e $filtered ]
    then
    awk -F' ' '($3>900){print}' $dest > $filtered
    fi
        
done < $taxID_name_file