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
 - I ran a test assembly of the 20150221 library with ABySS (using the default setting of k=64). Scripted with `abyss_assembly.qsub`. 


Wednesday, 19 August 2015

 - I ran the trimming script on the newly-replaced, md5-passing read files
 - I then ran an ABySS assembly on the 20150707_DNASeq_PE reads as per above


Thursday, 20 August 2015

 - using QUAST to get some metrics on the ABySS assemblies:
	- `python /opt/software/QUAST/2.3--GCC-4.4.5/quast.py -o /mnt/research/efish/2015_genomic_data/quast_results_ABySS64/ --eukaryote --gene-finding 20150707_DNASeq_PE/P_kings-unitigs.fa 20150721_DNASeq_PE/P_kings-unitigs.fa`
 - 1st indications suggest that we have quite different results from the two Illumina libraries: 814 vs. 267 contigs, 119 vs. 45 contigs >1kbp, but N50 of 737 vs. 800. I am going to double-check the pipeline upstream to convince myself that this is not *my* error...
	- this difference in information content is reflected in the size of the files (bytes or lines) in the two libraries (added 2 R scripts to plot these): ...21 seems to much smaller than ...07
	- I shall bring this up in my meeting with JG tomorrow


Tuesday, 1st September 2015

  - FTP-ed more genomic data from the NAS drive to the HPC; `MSUEFISHLAB_DROPBOX_NEW/pkings_genome/*` files all pass md5 check.
  - 3 files from `/original_data/mormyriformes/paramormyrops/` FAILED md5 check. Problem files are:
	- `gerald_C5A4TACXX_3_CTTGTA.bam`
	- `gerald_C5A4TACXX_3_GCCAAT.bam`
	- `gerald_C5A4TACXX_3_GTGAAA.bam`
	- This is *not* a local transmission problem -- these files match those in the efish research space, and those on the NAS, but *do not* match the checksums in their associated `...bam.md5` files. JG is aware.

  - I also renamed the trimming and QC scripts to match the numbering of the pipeline planned in my meeting with JG on 21-aug-15.


Wednesday, 2 September 2015

  - I'm going to pool this years' Illumina data into one folder, since there seems to be no benefit in trimming them separately.
  - Moved my previous attemped assemblies into `Aug_Assemblies`.
  - `Quast` tool is probably not the most appropriate for this data (qv. discussion with Lex at NGS camp) -- ought to be using `CEGMA` & `BUSCO`.
  - generated a bowtie2 index from the *P. kings* genome v0.1 using `bowtie2-build -f supercontigs.fasta Pkings`
  - set up some test runs to estimate timing.


Thursday, 3 September 2015

  - dummy run results suggest that the assembly ought to require somewhere between 23 & 43 hrs... so a bit variable. Going to try 48hrs first.
  - updated Alignment script appropriately.




