#! /bin/bash

rm missing_files.txt

for a in *fastq.gz
  do
  if [ ! -f `basename ${a} _001.fastq.gz`_pe.trimmed.fq ]
    then echo `basename ${a} _001.fastq.gz`_pe.trimmed.fq >> missing_files.txt
  fi
  if [ ! -f `basename ${a} _001.fastq.gz`_se.trimmed.fq ]
    then echo `basename ${a} _001.fastq.gz`_se.trimmed.fq >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.dedup.bam ]
    then echo `basename ${a} .aligned.sam`.dedup.bam >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.realigned.bam ]
    then echo `basename ${a} .aligned.sam`.realigned.bam >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.recalibrated.bam ]
    then echo `basename ${a} .aligned.sam`.recalibrated.bam >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.recalibrated.boot.bam ]
    then echo `basename ${a} .aligned.sam`.recalibrated.boot.bam >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.raw_variants.vcf ]
    then echo `basename ${a} .aligned.sam`.raw_variants.vcf >> missing_files.txt
  fi
done

for a in *aligned.sam
  do
  if [ ! -f `basename ${a} .aligned.sam`.raw_variants.vcf.gz ]
    then echo `basename ${a} .aligned.sam`.raw_variants.vcf.gz >> missing_files.txt
  fi
done

for a in *aligned.sam
	do
	if [ ! -f `basename ${a} .aligned.sam`.recalibrated.boot.bam ]
	 then echo `basename ${a} .aligned.sam`.recalibrated.boot.bam >> missing_files.txt
	fi
done

for a in *aligned.sam
    do
    if [ ! -f `basename ${a} .aligned.sam`.raw_variants.boot.vcf ]
     then echo `basename ${a} .aligned.sam`.raw_variants.boot.vcf >> missing_files.txt
    fi
done

