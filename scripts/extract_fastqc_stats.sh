#! /bin/bash

# this script is intended to be run in a folder of `...factqc.zip` reports and extract the `fastqc_data.txt` file that holds the stats.

allzips=(*fastqc.zip)

numzips=`expr $(echo ${#allzips[*]}) - 1`

for i in `seq 0 ${numzips}`
  do 
	zipname=`basename ${allzips[${i}]} .zip`
	unzip -p ${i} ${zipname}/fastqc_data.txt > ${zipname}_data.txt
  done

