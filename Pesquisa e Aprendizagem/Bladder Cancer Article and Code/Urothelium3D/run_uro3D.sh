#!/bin/bash

a=$1"/"
mkdir -p $a

if [ ! -f "./"$a"tags_"$6".dat" ]; then

	./cpm_uro3D.out "$a" $2 $3 $4 $5 $6 # dir seed N1 J11 tumor id
else

	echo "WARNING: id $6 for directory $a already exists and has finished! In order to prevent previous data being overwritten, simulation will not run."
fi
