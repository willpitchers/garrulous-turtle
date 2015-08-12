# 2015_genomic_data folder opened by WRP on Wednesday, 12 August 2015
===

Wednesday, 12 August 2015

 - all contents of `20150707_DNASeq_PE`, `20150721_DASeq_PE` & `BioNano` directories uploaded from ftp://efishbackup.zoology.msu.edu
 - md5sum checked on all ...fasta.tar.gz files
 - git repo. initiated
 - added script for initial FastQC quality assessment: `step1_fastQC.qsub`
 - run with:
	- sed s/WhichDir/20150707_DNASeq_PE/g step1_fastQC.qsub | qsub
	- sed s/WhichDir/20150721_DNASeq_PE/g step1_fastQC.qsub | qsub



