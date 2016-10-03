#! /bin/bash

rm missing_files.txt
readarray Fastqs < fastq.list
readarray Readfiles < readfiles.list
readarray Fishfiles < fishfiles.list

for a in ${Fastqs[@]}
  do
  if [ ! -f ${a} ]
    then echo ${a}  >> missing_files.txt
  fi
done

for a in ${Readfiles[@]}
  do
  if [ ! -f ${a}.trimmed.fq ]
    then echo ${a}.trimmed.fq >> missing_files.txt
  fi
done

for a in ${Fishfiles[@]}
  do
  if [ ! -f ${a}.trimmed.aligned.sam ]
    then echo ${a}.trimmed.aligned.sam >> missing_files.txt
  fi
done

for a in ${Fishfiles[@]}
  do
  if [ ! -f ${a}.trimmed.aligned.sorted.sam ]
    then echo ${a}.trimmed.aligned.sorted.sam >> missing_files.txt
  fi
done

for a in ${Fishfiles[@]}
  do
  if [ ! -f ${a}.trimmed.aligned.dedup.bam ]
    then echo ${a}.trimmed.aligned.dedup.bam >> missing_files.txt
  fi
done

for a in ${Fishfiles[@]}
  do
  if [ ! -f ${a}.trimmed.aligned.dedup.realigned.bam ]
    then echo ${a}.trimmed.aligned.dedup.realigned.bam >> missing_files.txt
  fi
done

for a in ${Fishfiles[@]}
  do
  if [ ! -f ${a}.trimmed.aligned.dedup.realigned.recalibrated.bam ]
    then echo ${a}.trimmed.aligned.dedup.realigned.recalibrated.bam >> missing_files.txt
  fi
done

#for a in ${Readfiles[@]}
#  do
#  if [ ! -f ${a}.recalibrated.boot.bam ]
#    then echo ${a}.recalibrated.boot.bam >> missing_files.txt
#  fi
#done

#for a in ${Readfiles}
#  do
#  if [ ! -f ${a}.raw_variants.vcf ]
#    then echo ${a}.raw_variants.vcf >> missing_files.txt
#  fi
#done

#for a in ${Readfiles}
#  do
#  if [ ! -f ${a}.raw_variants.vcf.gz ]
#    then echo ${a}.raw_variants.vcf.gz >> missing_files.txt
#  fi
#done

#for a in ${Readfiles}
#	do
#	if [ ! -f ${a}.recalibrated.boot.bam ]
#	 then echo ${a}.recalibrated.boot.bam >> missing_files.txt
#	fi
#done

#for a in ${Readfiles}
#    do
#    if [ ! -f ${a}.raw_variants.boot.vcf ]
#     then echo ${a}.raw_variants.boot.vcf >> missing_files.txt
#    fi
#done

