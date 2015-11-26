# 2015_genomic_data folder
---
*opened by WRP on Wednesday, 12 August 2015*

Wednesday, 12 August 2015

 - all contents of `20150707_DNASeq_PE`, `20150721_DASeq_PE` & `BioNano` directories uploaded from the lab NAS drive.
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
  - `Quast` tool is probably not the most appropriate for this data (qv. discussion with Lex at NGS camp) -- ought to be using [`CEGMA`](http://korflab.ucdavis.edu/datasets/cegma/#SCT5) & [`BUSCO`](http://busco.ezlab.org/#software).
  - generated a bowtie2 index from the *P. kings* genome v0.1 using `bowtie2-build -f supercontigs.fasta Pkings`
  - set up some test runs to estimate timing.


Thursday, 3 September 2015

  - dummy run results suggest that the assembly ought to require somewhere between 23 & 43 hrs... so a bit variable. Going to try 48hrs first.
  - updated Alignment script appropriately.


Friday, 4 September 2015

  - Alignment still running (no crash yet)
  - I built a second bowtie index from the 'superscaffold' BioNano file. I intend to try alignments using this map as a reference also.


Monday, 14 September 2015

  - I did not update this file as I ought to have done due to illness. During the week of the 7th--11th I ran scripts to align all the Illumina data to the BioNano 'superscaffold' bowtie index, and to sort and index both the `sam` alignment files, and to convert to bam format.
  - These scripts don't appear to have performed as expected. After poking at the logfiles, it seems that I had accidentally loaded an older version of SAMtools, and the options have changed subtly... I see the problem. Scripts updated and resubmitted.


Tuesday, 15 September 2015

  - The indexing scripts have errored-out overnight: I clearly **didn't** diagnose the problem correctly.
  - I *think* that I have located the problem. In addition to changing the output flags in `samtools view` from SAMtools/0.* to SAMtools/1.2, the devs have also made the -T flag (prefix for temp file names) complusory in `samtools sort`. Scripts are re-submitted.


Thursday, 17 September 2015

  - Indexing the (huge) sam files is taking even longer than I projected -- the last scripts ran out of walltime. I have resubmitted with more time requested.


Friday, 18 September 2015

  - In light of my repeated underestimates for how long the SAMTools processes ought to take, I am going to split up the bam -> sam conversion, bam sorting and bam indexing steps into their own scripts so that I can request walltime for each process separately.
  - Sam -> Bam conversions have run.
  - Bam sorting is underway


Monday, 21 September 2015

  - Bam sorting completed. Bam indexing running.
  - Bam indexing ran surprisingly fast – completed in ~2hrs, vs. ~30hrs for sam->bam conversion and another ~30hrs for bam sorting.
  - Reorganized SAMTools scripts as above.


22-23 September 2015

  - read documentation for [GATK](https://www.broadinstitute.org/gatk/guide/) and worked through their suggested 'best practice workflow' with 2 of our datafiles
  - met with JG to discuss how GATK might help us reach our immediate goals... we decided to alter the pipeline to get our data into the most GATK-compatible formats
  - I made an attempt to get some data on how a GATK-generated `.vcf` file compares with one output by SAMTools `mpileup`, but VCFtools is fighting me on this...


24-25 September 2015

  - Thursday was mostly eaten by a course at Beacon/iCER and I spent Friday moving equipment from Ian's old lab...


Monday, 28 September 2015

  - planning out new pipeline to fit around the GATK. req.s...
  - started writing array scripts to handle over-large arrays and auto-resubmission.


Tuesday, 29 September 2015

  - testing smart array script for trimming steps
  - success with self-(re)-submitting script!


Thursday, 1st October 2015

  - Re-downloaded the `gerald_*.bam`s from the Genomics Centre to (hopefully) fix the files-with-non-matching-checksums problem.


Saturday, 3rd October 2015

  - Fixed up *another* bug in my alignment script -- now re-running.


Sunday, 4 October 2015

  - Found a *bizarre* bug whereeby the colour-coding in my shell settings was causing `bwa mem` to add colour-code escaped characters to the outputed `..sam` files (but only sometimes because it wants me to doubt my sanity). I have edited my script to pass the call to named variables in the `bwa mem` command through `sed` to trip out any color-codes that make it that far.
  - Now re-re-running.


Wednesday, 7 October 2015

  - FTP-ed the replacement `.bam` files from the replacement 'Original_data' set to the lab NAS drive for cold storage.
  - pipe has run unsupervised through trimming, alignment, indexing, sorting, deduping, indel realignment and base recalibration! ...but some files seem to have been dropped on the way.
  - I wrote `find_missing_files.sh` to list what's missing... 3 fastq files are weird in a way that breaks Picard. Investigating...


Thursday, 8 October 2015

  - I think that I have found the source of (some) of the problems -- I have truncated `.fastq` files in some cases, with the EOF appearing in the middle of a 4-line fastq entry.
  - I'm investigating whether I've introduced this error and if so when..


Monday, 12 October 2015

  -  Over the weekend I ran a couple of the suspect files through the pipeline again, step-by-step and running interactively to make sure that I could catch any errors. It seems that for the libraries that were 'leaking' the problem was occurring right at the start of the pipe; something caused the `gunzip` process to stop part-way through the extraction from the gzip archives.
  - In order to be sure that there aren't other failure points in the pipe I have thrown out all my output files thus far and started from scratch.


Tuesday, 13 October 2015

  - Finished the from-scratch re-run. All output files now appear to be present and correct!


Thursday, 15 October 2015

  - We have `.vcf` files!
  - Some of the larger (PE) files overran my 4hr estimate at the `HaplotypeCaller` stage, upped the resource request to 6hrs


Friday, 16 October 2015

  - merging the `.vcf` files using `vcf-merge` from the `VCFtools` toolkit isn't working as expected. Investigating...


Tuesday, 20 October 2015

  - I think that I've finally found the bug that was causing the merged vcf file to have 0 lines... testing fixed script.


Wednesday, 21 October 2015

  - I'm going to split the vcf merging into two steps for speed, since I can parallelize the first step.


Thursday, 22 October 2015

  - Irritatingly, the way I split the vcf-merging job caused the 'broken pipe' bug to return.
  - I have also realized that it would be smarter to use the 'remove-duplicates' flag for the merge *within* samples (i.e. across libraries & lanes). Then the second merge step *among* samples can be run as is since the duplicates between samples are informative... To this end I have re-written the vcf-merge scripts *again*.


Monday, 26 October 2015

  - I'm giving up on `vcf-merge` from the `VCFtools`... re-scripting using `CombineVariants` tool from GATK.


Wednesday, 28 October 2015

  - Added script to run Fst calculations for all possible pairwise population comparisons


Thursday, 29th October 2015

  - Adjusted the walltime and mem requirements of a few scripts based on how much time/memory they actually used -- keeping these demands realistically low will help me spend less time queueing.

Week of 2nd-6th Nov.

  + Minimal progress due to admin. issues not related to the project.
  + Preliminary analysis of Fst output -- `Fst_plots.Rmd` added to project folder


Week of 9th-13th Nov.

  - adjusted Fst script to enable window_size & window_step adjustment -- submit with e.g. `sed s/WSIZE/x/ calc_Fst_array.qsub | sed s/WSTEP/y/ | qsub` for size & step == x & y respectively (in bp).
    - now run at:
      + 100kb window -- 50kb step
      + 100kb window -- 100kb step
      + 500kb window -- 100kb step
      + 500kb window -- 50kb step
      + 200kb window -- 100kb step
      + 200kb window -- 50kb step
  - re-scripted the base_score_recalibration_array as base_score_recal-boot_array to take the final, merged .vcf file as an input for the 'known SNPs' field in `BaseRecalibrator`. This forms the loop that we'll need to bootstrap ourselves a high-confidence variant call.
    + looped once... waiting for scripts to run...
      - not all scripts have finished properly... 1st order of business on Monday will be to work out why.
  - added `AnalysisPipeline.md` file to document the pipeline in order of operation (rather than chronologically as here) for ease of reference.


Week of 16th-22nd November

  - restarted base_score_recal-boot_array as some jobs appear to have over-run their walltime
    + added walltime may bot have been sufficient for all the jobs in the array... investigating.
  - updated `Fst_plots.Rmd` to show effect of using different window/step sizes
  - writing scripts to implement GWAS with plink & SKAT...
  - the weather has turned cold -- this appears to have broken the HPC as is traditional


Week of 24th-28th November

  - NB: this week includes Thanksgiving & an HPC shutdown...
  - Downloaded the VecScreen database of bacterial sequences to check for potential contaminants in our genome...
    - built a blast DB for P_kings-supercontigs & blasted all ~5000 VecScreen sequences against it
    - 
