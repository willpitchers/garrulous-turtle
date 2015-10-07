#! /bin/bash

touch missing_files.txt

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
