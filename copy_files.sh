#!/usr/bin/env bash

region="CG DLPFC NAC OFC SUB"
prefix="5-Scripts-DrugFindR"

for r in $region; do
    files_10=$(fd 10% "$prefix/$r")
    files_5=$(fd 5% "$prefix/$r")
    for f10 in $files_10; do
        # Add an a after the first number in filename
        newname=$(echo ${f10/10%/50%} | sed -E 's/([0-9]+)-drug/\1a-drug/')
        echo "Copying ${f10} to ${newname}"
        cp -v $f10 $newname
        ./update_names.sh $newname
    done
    for f5 in $files_5; do
        newname=$(echo ${f5/5%/25%} | sed -E 's/([0-9]+)-drug/\1a-drug/')
        echo "Copying ${f5} to ${newname}"
        cp -v $f5 $newname
        ./update_names.sh $newname
    done
done