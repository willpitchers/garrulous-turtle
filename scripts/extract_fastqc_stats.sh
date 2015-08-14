#! /bin/bash

# this script is intended to be run in a folder of `...factqc.zip` reports and extract the `fastqc_data.txt` file that holds the stats.
# added: bundle up the reulting txt files into an archive for ease of backup, then remove them to keep directory clean 

allzips=(*fastqc.zip)

numzips=`expr $(echo ${#allzips[*]}) - 1`

for i in `seq 0 ${numzips}`
  do 
	zipname=`basename ${allzips[${i}]} .zip`
	unzip -p ${zipname}.zip ${zipname}/fastqc_data.txt > ${zipname}_data.txt
  done

dirname=`basename $(pwd)`

zip ${dirname}_qc.zip *_data.txt

rm *_data.txt
