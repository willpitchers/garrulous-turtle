# 2015_genomic_data folder
---
*opened by WRP on Wednesday, 12 August 2015*

Wednesday, 12 August 2015

 - all contents of `20150707_DNASeq_PE`, `20150721_DASeq_PE` & `BioNano` directories uploaded from ftp://efishbackup.zoology.msu.edu
 - md5sum checked on all ...fasta.tar.gz files
 - git repo. initiated
 - added script for initial FastQC quality assessment: `step1_fastQC.qsub`
 - run with:
	- `sed s/WhichDir/20150707_DNASeq_PE/g step1_fastQC.qsub | qsub`
	- `sed s/WhichDir/20150721_DNASeq_PE/g step1_fastQC.qsub | qsub`

Thursday, 13 August 2015

 - wrote script to extract fastqc stats from reports: 512 reports is too many to check by eye so i intend to script the process


Friday, 14 August 2015

 - wrote an R script to make plots of per-bp quality score across all the read-files
 - run with:
	- `unzip 20150707_DNASeq_PE_qc.zip -d 20150707_DNASeq_PE_qc ; cd 20150707_DNASeq_PE_qc ; R --no-save < ../scripts/plot_bp_fastqc.R`
	- `unzip 20150721_DNASeq_PE_qc.zip -d 20150721_DNASeq_PE_qc ; cd 20150721_DNASeq_PE_qc ; R --no-save < ../scripts/plot_bp_fastqc.R`
 - wrote an R script to make GC-content plots
 - run with:
  - `unzip 20150707_DNASeq_PE_qc.zip -d 20150707_DNASeq_PE_qc ; cd 20150707_DNASeq_PE_qc ; R --no-save < ../scripts/plot_GC_content.R`
  - `unzip 20150721_DNASeq_PE_qc.zip -d 20150721_DNASeq_PE_qc ; cd 20150721_DNASeq_PE_qc ; R --no-save < ../scripts/plot_GC_content.R`


Monday, 17 August 2015

 - trying out 'SGA preprocess' from the 'SGA' toolkit (https://github.com/jts/sga/)
 - scripted a simple adapter trim with a 10bp headcrop after looking at the fastQC
 - NB: these Illumina files are coded in **phred33**
 - trimming scripts are crashing. Not clear why...


Tuesday, 18 August 2015

 - *Now* I see the problem -- sequences in the 20150707_DNASeq_PE now *fail* md5 check. Curious. I shall replace them with the originals from: `wget -r --user="pitchers" --password="" ftp://efishbackup.zoology.msu.edu:/Data/20150707_DNASeq_PE/*.gz`
 - I ran a test assembly of the 20150221 library with ABySS (using the default setting of k=64)

