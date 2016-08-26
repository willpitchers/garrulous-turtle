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
  - Moved my previous attempted assemblies into `Aug_Assemblies`.
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
      + 100kbp window -- 50kbp step
      + 100kbp window -- 100kbp step
      + 500kbp window -- 100kbp step
      + 500kbp window -- 50kbp step
      + 200kbp window -- 100kbp step
      + 200kbp window -- 50kbp step
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
    - ~5500 results files; I used `grep -L  "No hits found" results* > vec_hits.txt` to list just those VecScreen sequences with hits in P_kings-supercontigs.
    - hits on Scaffolds 1910, 3315, 3294 & 24, all with tiny 'E' values.
      + Scaffolds 1910, 3315 & 3294 comprise only 4-5kbp each (and bear no SNPs), so I advocate simply discarding them
      + Scaffold 24 is *much* larger at 3.8Mbp, so it'd be better not to lose it (or the 36k SNPs it contains) completely... lets look deeper:
        + the matches on scaf1910 are both long (>1kbp) and near-perfect (99% & 100% with 0% gaps!)
          - ID= gnl|uv|J02482.1:1-5386-49 Coliphage phi-X174
        + the matches on scaf3315 & scaf3294 are shorter (both 145bp) and somewhat less convincing (80% with 6% gaps)
          - ID= gnl|uv|JN874480.1:6795-9434 Cloning vector pHUGE-Red
        + the match on scaf24 is the shortest (83bp) and similarly strong (85% with 2% gaps)
          - ID= gnl|uv|AY390769.1:1-1000  Synthetic construct BigDye Terminator Cycle Sequencing Standard sequence... so pretty likely contamination then.


Week of 30th November - 4th December

  - 'scriptified' the analysis (above) of VecScan results to ensure reproducibility
  - further work on Fst analysis


Week of 5th - 11th December

  - this week my laptop's HDD died: backups should prevent any loss of data, but it wasted a lot of time
  - I was able to do *some* work on the [Fst analysis](./Fst_plots.Rmd)


Week of 14th – 18th December

  - I have persuaded the base quality score recalibration scripts to run for the 1st iteration of the bootstrap. The sequence files needed dramatically different amounts of time/memory for the run, but that doesn't explain all of the crashes...
  - After discussion with JG, we're going to take a look at Fst at smaller window sizes. Now running at:
    + 5kbp window -- 5kbp step
    + 10kbp window -- 10kbp step
    + 10kbp window -- 5kbp step
  - some bug fixes for vcf-merge post bootstrap


21st December – 1st January

  - Xmas break week -- conveniently enough the iCER staff have partially shut down the HPC to fix some problem with scratch... inconveniently, that is where all our data lives.


Week of 4th – 8th January

  - Fst results at smaller window sizes look to be giving us a much more high-resolution picture of the genotypic differences – as expected...
  - updated `Fst_plots.Rmd` to include comparisons of newer vs. older results
  - after discussion with JG, filtered Fst results based on:
    - bins that have high Fst in *both* comparisons between fixed penetrating & fixed non-penetrating pop.s...
    - ...and have *low* Fst in the comparison between the phenotypically similar pops...
    - ...and have elevated Fst in the comparison between the *mixed* populations.
  - `Fst_plots.Rmd` updated to include tables of the top 10% and 5% of the filtered 'hits' from this schema


Week of 11th – 15th January

  - Looking into running an GWAS with the [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/index.shtml) toolkit.
  - data needs to be converted into the right formats – scripted array jobs for this using vcftools
  - problem: phenotype data needs to be added to `.ped` files *before* converting to `.bim`/`.fam`, otherwise PLINK throws an error message... wrote an AWK program to do this, based on pop-level phenotype summary in `pheno_pop.txt`


Week of 18th – 22nd January

  - ran an genome-wide association with Fisher's exact test using PLINK v.1.07 (NB: preliminary, since phenotypes are coded on a pop. level but the data is for *individuals*)
  - wrote `PLINK_assoc_results.Rmd` report


Week of 25th – 29th January

  - re-ran the genome-wide association with the bootstrapped `.vcf` file
  - wrote `PLINK_assoc_results_PostBS.Rmd` report


Week of 1st – 5th February

  - rewrote `PLINK_assoc_results_PostBS.Rmd` report
  - worked out how to get `BreakDancer` running...
  - `BreakDancer` run appears not to have finished cleanly... the problem seems to be in the `.bam` files. Investigating...


Week of 8th – 12th February

  - so BreakDancer seems to fail due to a missing end of file (EOF) marker on one of the bam files.
      - I'm backtracking along my pipeline to see when this got dropped...
      - I now suspect this may the result of PICARD & SAMTools disagreeing about what a `.bam` file *ought* to look like...
      - I'm running a test of BreakDancer on a `.bam` file that SAMTools hasn't touched.
  - BreakDancer repeatedly fails to finish the run, still blaming EOF marker problems...
  - while I'm trying to discover the root of the BAM EOF problem, I'm also setting up jobs for a 2nd GATK vcf-bootstrap


Week of 15 - 19th February

  - 2nd GATK vcf-bootstrap continues
  - recreating the `.bam` files from the start of the pipeline to try and fix the EOF problem...
    - no luck... time to come up with plan E.
  - split up the operations involved in the plink analysis to make the jobs shorter and easier to schedule
  - downloaded Gar genome from [ensembl](http://useast.ensembl.org/Lepisosteus_oculatus/Info/Index?redirect=no) for coding sequences to BLAST against our scaffolds... hopefully identification of syntenic regions will help us bridge some scaffold gaps


Week of 22nd – 26th February

  - JG provided spreadsheet with individual phenotypes: `Dropbox/Projects/Paramormyrops_Reseq/specimens_for_genome_reseq.xlsx` saved as `Dropbox/WILL/genome_assembly/garrulous-turtle/specimens_for_genome_reseq.txt`
  - discovered a suspected [typo](https://msuefishlab.slack.com/files/willpitchers/F0NLE0963/specimens_for_genome_reseq_wp.xlsx) in the phenotype file.
    - confirmed by [JG](https://msuefishlab.slack.com/archives/D0G6C7RJM/p1456247717000003)
    - simplest fix is to add a line to `plink_ind_pheno.qsub` to use the (incorrect) ID number to associate the correct phenotype with the mis-named sequence files. (`sed -i s/4923/3923/ individual.list`)
  - added `plink_ind_pheno.qsub` and `pheno_lookup.py` to scripts to read the individual-level phenotype data from `specimens_for_genome_reseq.txt` into the existing PLINK GWAS analysis.
    - set GWAS mk.2 running 23/02/2016
      - this run failed seemingly because of a *2nd* [typo](https://msuefishlab.slack.com/files/willpitchers/F0NS88DHU/specimens_for_genome_reseq_wp.xlsx)...
    - set GWAS mk.3 running 24/02/2016
      - *this* run failed with `FEXACT error 3` message... message appears unknown to the FAQ.
  - trying plan E re: BreakDancer `.bam` EOF problem...


Week of 29th Feb. – 4th March

  - BreakDancer still failing... got as far as sample BAVA_6623. I am going to back-track to the recalibration step to (hopefully) find where the file got truncated/lost its EOF.
    - there seem to be quite a few of the original pool of 768 bam files that left the recalibration step with a malformed EOF (but not severely enough to break the variant caller)... checking them all.
  - the HPC is replacing a rack of storage, so there's been much waiting/crashing.


Week of 7th – 11th March

  - Traveling...
  - HPC is *barely* functioning...


Week of 14th – 18th March

  - HPC is still barely functioning...


Week of 21st - 25th March

  - prepping for ICN2016
  - Traveling...


28th March – 7th April

  - ICN2016, efish satellite meeting, and trapped in Uruguay


Week of 11th – 15th April

  - Make subset `.vcf` for JG:
    - output list of 'hits' from `PLINK_assoc_results.Rmd` -> `top_hit_variants`
    - `vcftools --gzvcf all_variants_merged_27_10_2015.vcf.gz --out "tophits_" --positions ${DIR}top_hit_variants --recode`
  - rescripting the vcf-merging step to use the `CatVariants` tool in place of the `CombineVariants` tool (both from GATK) to try and get the vcf to list *individuals* rather than *libraries*.
  - having so much difficulty with the HPC that I'm going to try and set up an EC2 instance...
    - `wget https://www.dropbox.com/s/q5bwd3n54sed4cq/EC2_setup_script.sh` and run to set up machine
    - trying on an c3.4xlarge, Ubuntu 14
  - alternate solution: trying to `BLCR longjob` my stuff to help it slide into small queue gaps, then snapshot itself and hop back in the queue...


Weeks of 18th - 22nd & 25th - 29th April

  - focused on manuscript preparation


3rd - 11th May

  - Will in Seattle


Week of 16th - 20th May

  - Extracted genotypes from the 'hits' vcf file with `vcftools --vcf top100.recode.vcf --extract-FORMAT-info GT --out GT_top100`
    - `sed s/"\.\/\."/"0\/0"/g GT_top100.GT.FORMAT | sed s/"0\/0"/"R"/g | sed s/"0\/1"/"H"/g | sed s/"1\/1"/"A"/g > hits.geno`
  - trying to test PLINK sensitivity using simulated data:
    - build fake vcf rows -> append to subset of real vcf file
    - use vcftools to convert altered vcf file -> plink files
    - use awk to edit phenotypes in doctored plink files –> run plink
  - We think that a proportion of the apparently missing-data codes in the multi-individual vcf files may in fact represent homozygous loci for the reference allele... [seqanswers](http://seqanswers.com/forums/archive/index.php/t-28325.html) suggests that this is a know 'feature' of GATK.
    - a few hours of poking through manuals/googling suggests that the best way to fix this is to detour around the problem by moving the coalescence-of-individuals point further up the pipeline.
    - 1st attempt: coalesce individuals as recalibrated BAM files immediately prior to vcf discovery run..
      - NB: using `samtools merge` sequentially to merge e.g. 1.bam + 2.bam + 3.bam gives identical (i.e. no `diff`) results to merging 1.bam + 2.bam, followed by merging 1+2.bam + 3.bam (tested on Shockly with SAMtools/0.1.19)


Week of 23rd – 27th May

  - with reference to [AnalysisPipeline.md](./AnalysisPipeline.md):
    - the first build of the pipeline was pleasantly parallel for steps 1–8 (mostly 768-job arrays) before merging at the vcf stage...
    - updated version is going to merge bam files immediately after recalibration (step 5), merging in stages; libraries within individuals, then individuals within populations, then all populations
  - re-testing PLINK sensitivity – last week the results I was getting made no sense...
    - rebuilt fake vcf rows -> re-append to subset vcf -> edit vcf to replace incorrect version code (HPCC Ticket #20345 I think) with `sed -i s/"fileformat=VCFv4.2"/"fileformat=VCFv4.0"/ pedrows.vcf` -> covert vcf to ped using `vcftools --vcf pedrows.vcf --out pedrows --plink` -> insert phenotype data as per [plink_pheno](scripts/plink_pheno.qsub) -> run `plink --file pedrows --allow-no-sex --assoc fisher --out pedtest --allow-extra-chr`
    - this procedure now seems to work... I may need to simulate some more fake data at closer allele frequency intervals...
  - long drawn-out wrestling with the HPC trying to get the "longjob" checkpoint/restart tool to work...


Week of 30th May - 3rd June

  - talked with JG, decided that the appropriate 'level' at which to run the HaplotypeCaller is the individual, i.e. `HaplotypeCaller -I APA_6675_all_libraries.bam.g.vcf`, then run `GenotypeGVCFs` on all 63 `...g.vcfs`
  - progress made towards useable `longjob` script for checkpoint-&-restarting
    - successfully checkpointed/restarted a pointless command
    - successfully checkpointed at GATK job, but the restart does not work consistently...
    - wrote to iCER to ask for help, also contacted Dirk Colbry (who originally wrote the `longjob` powertool)
  - Due to continued HPCC problems, I've copied the individual-level BAM files to Shockly and started running GATK HaplotypeCaller (with the GVCF option) there as a back-up...


Week of 6–10th June

  - The restart problem is temp-file related. Dirk and I have exchanged a number of emails about this...
  - tried editing script to search for temp-files on restart...
  - tried editing script to mkdir its own folder to hold temp-files...
  - we've reached the point where I, Chun-Min *&* Dirk are all stumped as to why we cannot get a clean restart :-(
  - in order to get the GVCFs made in good time I'm playing with time and memory parameters and GATK's built-in multi-threading...


Week of 13-17th June

  - most of the week in crunch-mode on the efish special issue MS...
  - HaplotypeCaller still running on Shockly. 4 of 63 GVCFs finished so far.
  - As installed on the HPC, GATK can run multi-threaded at the *core* but not the *node* level, i.e. I can multi-thread across a single multi-core node.
  - testing suggests that the jobs GVCF jobs should run for ~15hrs with 20GB memory on a quad-core node
    – array submitted to the queue for 20hr jobs – all surviving nodes have >=4 cores and >=20GB so this ought to minimise queue-time
    – all jobs still queueing after 2 days. `showq` and `showstart` appear to have been broken by the HPC rebuild


Week of 20-24th June

  - most of the week in crunch-mode on the efish special issue MS.
  - 50 of 63 20hr jobs complete, re-submitted the remaining 13 for 30hrs with 4cores and 20GB
  - 30hr jobs queueing for days... new (intel16 'Laconia') nodes are all big-memory: now that they're online does the trade-off of requesting more cores/memory vs. time change?
    - testing suggests that with double memory these jobs should only take ~6hrs, so shouldn't queue long - testing 60GB on 4 cores for 8hrs
	- more testing suggests that with double memory *and* threads these jobs should only take 3-4hrs - testing 60G on 8 cores for 4hrs
		- both these tests run out of time for the 'awkward 13' individuals, but the time remaining on the GATK clock when the jobs get killed is >1hr in most cases...
		- paging through output files suggests that rate of progress plateaus part-way through the job – this limits the value of my testing
	- set 8hr, 8core, 100GB jobs running on the awkward 13
  - HaplotypeCaller still running on Shockly. 9 of 63 GVCFs finished so far.


Week of 27th June – 1st July

  - 8hr, 8core, 100GB jobs completed properly for 9 of 13 awkward individuals – still had 4 individuals run out of time, but with *seconds* left on the clock in 3 cases, the 4th had ~20mins
  - re-submitted 10hr, 8core, 100GB jobs for the final 4.
  - meanwhile, testing shows that `GenotypeVCFs` crashes out in ~10mins with currently-available GVCFs... investigating.
    - 1 individual that I thought had run correctly did not – possibly due to kernel panic and reboot of the machine that houses our research space? ...re-running.
  - testing suggests that `GenotypeGVCFs` may take longer than a week! ...I'm working on getting that number down.
	- preparation of alignment rate stats:
      - `for i in *bam ; do echo ${i} >> alignment_rates.txt  ; samtools flagstat ${i} >> alignment_rates.txt  ; done`
      - `grep bam alignment_rates.txt > ARnames.txt`
      - `grep --color='never' "mapped (" alignment_rates.txt > ARnums.txt`
      - run through `P_kingsleyae_Alignment_Rates.Rmd` to produce summary.


Week of 4th-8th July

  - `GenotypeVCFs` tests still erroring out...
    - tested GVCF file individually – `APA_6680..`, `APA_6683..` & `APA_6684..` seem to be the sources of the problem
    - `APA_6680..` -> `ERROR MESSAGE: Line 98481725`. what happens if I excise the problem line?
      - running the first 9481720 rows only seems to work, but there is some weird non-human-readable stuff around the problem line... I can only assume that this is scrambled output thanks to one of the scheduler hiccups
      - are there other potential problems? I have a copy of `APA_6680..` that was calculated on Shockly – I am running a `diff` between the 2 now...
  - NB: `sed -i '0,/VCFv4.2/s//VCFv4.1/' ..bam.g.vcf` required before `igvtools index ..bam.g.vcf`


Week of 11th-15th

  - A **plethora** of HPCC problems...
    - `GenotypeVCFs` job crashing out with "malformed GVCF file" error message
      - validating the files with `GenomeAnalysisTK.jar -T ValidateVariants -R supercontigs.fasta --validationTypeToExclude ALL -V ..g.vcf`
    - I/O problem with "ERROR MESSAGE: Timeout of 30000 milliseconds was reached while trying to acquire a lock on file" while running as a scheduled job, but I cannot recreate this interactively.
      - running interactively on a dev-node suggests that adding the `--disable_auto_index_creation_and_locking_when_reading_rods` option suggests this might fix the issue...
    - `GenotypeVCFs` job ran in ~6hrs, but dropped 3 `..g.vcf` files... not sure why. However...
    - ...re-running with all 63 files timed out after 7hrs with ~18days left on the GATK clock(!)
    - re-re-running `GenotypeVCFs` job with longer walltime, but it crashes out with `file 'BAVA_6627_all_libraries.bam.g.vcf' does not exist`. This file **definitely** exists. :-(
    - waiting on my turn at running this job interactively (in debug mode) so that I can maybe diagnose the error...
  - copying all the `..g.vcf` files over to Shockly so I can run there...
    - in addition to md5sum-checking the transfer, I'm double-checking files first with `java -Xmx30g -jar /home/GATK/GenomeAnalysisTK.jar -T ValidateVariants -R supercontigs.fasta --validationTypeToExclude ALL -V APA_6675_all_libraries.bam.g.vcf`
  - on Tuesday 12th I ran out of patience and went to iCER to ask for help. Chun-Min Chang sat with me and we re-ran my tests together. He said that the most recent failure of my job was "probably due to a problem we had with scratch last night"
  - BREAKTHROUGH: the job appears to have run correctly in the small hours of the morning Wednesday 13th!
    - running PLINK pipeline: "File contains 31896675 entries and 63 individuals"
      - re-running the analysis
      - `vcftools --vcf all_individuals_12_07_16.out.vcf --out "tophits_all_individuals_12_07_16.out" --positions ${DIR}top_hit_variants --recode`
      - `vcftools --vcf tophits_all_individuals_12_07_16.out.recode.vcf --out "tophits_all_individuals_12_07_16" --freq`
      - output the genotype-by-individual (GT) format using `vcftools --vcf tophits_all_individuals_12_07_16.out.recode.vcf --extract-FORMAT-info GT --out all_individuals_12_07_16`
      - reformat this output using `cat all_individuals_12_07_16.GT.FORMAT | sed s/1\\/1/2/g | sed s/0\\/0/0/g | sed s/0\\/1/1/g | sed s/\\.\\/\\./NA/g > all_individuals_12_07_16.GT.FORMAT.recode`

`grep "Scaffold81" all_individuals_12_07_16.20x4.out.vcf > scaf81.list`

`for i in {0..4667} ; do grep "Scaffold${i}[[:space:]]" all_individuals_12_07_16.out.vcf | tail -1 > scaf.test ; done`


  - 15 July vcf/PLINK testing: association analysis seems to have some weirdness – SNPs called that are at coords > length of scaffold...
    - `head -9000 all_individuals_12_07_16.out.vcf > plinky/test.vcf && tail -1000 all_individuals_12_07_16.out.vcf >> plinky/test.vcf` gives me a 10000 row vcf file with 5230 sites and all 63 individuals
    - `java -Xmx20g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T ValidateVariants -R ../../P_kings_genome/supercontigs.fasta --validationTypeToExclude ALL -V test.vcf` confirms that the vcf is valid from GATK's POV
    - `vcftools --vcf test.vcf --out test --plink` makes .ped and .map files. These seem to be correctly formed.
    - `test.map` has four cols as it ought to: Chromosome, 'rs#', genetic distance, BP coord... I'm going to try recoding the 'rs#' col to include the scaffold ID too: `awk -F "[\t ]" '{ $2=$1 "_" $2; print $1 "\t" $2 "\t" $3 "\t" $4 }' test.map > taco.map && mv taco.map test.map`
    - `pheno_lookup.py` seems to be working as intended, made an edit to ensure that output is tab-delim.
    - after making binary files: `plink --file test --out test --make-bed`, tried plinking: `plink --bfile test --allow-no-sex --maf 0.1 --geno 0.5 --fisher --out test`...
    - the `test.assoc.fisher` output seems to be the right format – 'CHR' has been coded as '0' throughout, but the SNP IDs read e.g. 'Scaffold4662_25' – now testing with the *real* all_individuals_12_07_16.out.vcf
    - both the Shockly & HPC versions of the `..assoc.fisher` file *do not* contain either of the errors that JG identified (!?!)   <scratches head>
    - the output files (after plotting and labelling etc. on my laptop) *also* do not contain the errors – evidence points to the problem happening on the way into JG's copy of excel maybe?
    - `for i in {0..4667} ; do grep "Scaffold${i}[[:space:]]" all_individuals_12_07_16.out.vcf | tail -1 | cut -f1,2 >> scaf.test ; echo Scaffold${i} ; done` followed by s



Week of 18th-22nd July

  - without COB: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ../P_kings_genome/supercontigs.fasta -V all_individuals_12_07_16.20x4.out.vcf -o without_COB_12_07_16.20x4.out.vcf -sn 6675 -sn 6676 -sn 6677 -sn 6678 -sn 6679 -sn 6680 -sn 6681 -sn 6682 -sn 6683 -sn 6684 -sn 6685 -sn 6737 -sn 6494 -sn 6496 -sn 6497 -sn 6498 -sn 6499 -sn 6500 -sn 6501 -sn 6502 -sn 6597 -sn 6598 -sn 6599 -sn 6602 -sn 6603 -sn 6604 -sn 6605 -sn 6619 -sn 6620 -sn 6621 -sn 6622 -sn 6623 -sn 6624 -sn 6625 -sn 6626 -sn 6627 -sn 3923 -sn 4816 -sn 4832 -sn 4834 -sn 4893 -sn 4894 -sn 4895 -sn 4896 -sn 4897 -sn 4921 -sn 4925 -sn 6716 -sn 6717 -sn 6718 -sn 6719 -sn 6720 -sn 6721 -sn 6722 -sn 6723 -sn 6724 -sn 6725`
  - just APA/BAM: `java -Xmx20g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ../P_kings_genome/supercontigs.fasta -V all_individuals_12_07_16.20x4.out.vcf -o APA_BAM_only_12_07_16.20x4.out.vcf -sn 6675 -sn 6676 -sn 6677 -sn 6678 -sn 6679 -sn 6680 -sn 6681 -sn 6682 -sn 6683 -sn 6684 -sn 6685 -sn 6737 -sn 6494 -sn 6496 -sn 6497 -sn 6498 -sn 6499 -sn 6500 -sn 6501 -sn 6502 -sn 6597 -sn 6598 -sn 6599 -sn 6602 -sn 6603 -sn 6604 -sn 6605`


Week of 25–29th July

  -


Week of 1st-5th August

  - The goal of this week is to make sense of the lack of overlap among the lists of 'top hits' from the four different PLINK association runs...
  - bootstrap things...
  -

Week of 8th August

  - is the problem with the sam files?
    - `java -jar $PICARD/ValidateSamFile.jar I=APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.aligned.sam MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.aligned.sam MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6722_GAATTCGT-GGCTCTG_L003_R1_pe.aligned.sam MODE=SUMMARY` -> "no errors found"
  - ...the `..dedup.bam`?
    - `java -jar $PICARD/ValidateSamFile.jar I=APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.dedup.bam MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.dedup.bam MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6722_GAATTCGT-GGCTCTG_L003_R1_pe.dedup.bam MODE=SUMMARY` -> "no errors found"
  - ...the `..realigned.bam`? AHA! problem found.
    - `java -jar $PICARD/ValidateSamFile.jar I=APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.realigned.bam MODE=SUMMARY` -> "ERROR:INVALID_VERSION_NUMBER    1"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.realigned.bam MODE=SUMMARY` -> "ERROR:INVALID_VERSION_NUMBER    1"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6722_GAATTCGT-GGCTCTG_L003_R1_pe.realigned.bam MODE=SUMMARY` -> "ERROR:INVALID_INDEX_FILE_POINTER        1"
  - did the HPC have another I/O glitch maybe? let's try re-running these by hand? -- no dice: same problem.
  - GATK help forums suggest that this can happen as a result of older versions of PICARD vs. newer versions of GATK...
    - this doesn't seem to be the problem however :-(
    - I have now gone back to the start and run these 3 libraries from the raw reads and I'm *still* getting the same errors(!)...
    - `APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.realigned.bam` & `MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.realigned.bam` *do* have different version numbers than the other libraries in their `..bam` headers: "@HD     VN:1.5" rather than "@HD     VN:1.4" (using `samtools view -H`).
    - `samtools view -H APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.realigned.bam | sed 's/VN:1.5/VN:1.4/' | samtools reheader - APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.realigned.bam > APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.realigned.RH.bam` & `samtools view -H MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.realigned.bam | sed 's/VN:1.5/VN:1.4/' | samtools reheader - MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.realigned.bam > MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.realigned.RH.bam` *do* fix the version number in the headers and they *do* then return 'no errors' from `ValidateSamFile.jar`.
    - bizzarely, `SortSam`ing and `BuildBamIndex`ing fixes the problem with `MOV_6722_GAATTCGT-GGCTCTG_L003_R1_pe.realigned.bam`
  - I'm going to manually run the fixed versions through the rest of the pipeline so that I can run PLINK today, then circle back around to fixing the reproducibility problem...
    - while that's running, I'm going to clear out my 'Analysis' directory – I'll keep the raw reads, `..vcf`s and individual-level `..bam`s, but ditch everything else.
  - Now that I've found workarounds, re-running `bam_merge_samples_array.qsub` and `vcf_discovery_array.qsub` and `genotype_gvcf.qsub`...
  - Wednesday, 10 August 2016 – another I/O crashed the `vcf_discovery_array` jobs. The error logs claim that the reference file is absent, but list the full path to the reference, which bash can see... re-running.
  - ...aaaaand another crash. puzzled...
  - OK. finally got some GVCF jobs to run... but the post-bootstrap BAM files seem to have a more serious case of the non-linear-progression problem... testing indicates that they ought to run in ~2hrs, but *none* of them finish in 4hrs!
    - tried with more threads - up to 14 cores (largest available on a single node [NB: the `HaplotypeCaller` tool can't handle multi-threading across nodes] ) and still not fast enough...
    - tried with more memory - up to 500GB (largest available with 14 cores) and most still don't finish in 4hrs
    - https://msuefishlab.slack.com/files/willpitchers/F20RN8CTV/infinity.tiff <sigh>
    - due to the way the PBS is set up, I can demand the max available memory & #cores for 4 hrs and the jobs will start immediately, but sadly I'm going to have to run for longer and that means queueing for *days*.


Week of August 15th

  - meeting with JG about comparison between V0.2 & V0.1 pipelines
    - we have identified a potential problem with the specification of readgroups from BWA being passed along to GATK - different fields are required and the meaning of the fields is *very* specific... working to update this now.
    - [multi-sampled & multiplexed read tags in GATK](https://software.broadinstitute.org/gatk/guide/article?id=6472) are best extracted from the read *names*. Our read names look like this: `@HWI-D00731:65:C6UU2ANXX:1:1101:2046:1852 1:N:0:GAGATTCCTATAGCC`, wherein:
      - `HWI-D00731` is the instrument name
      - `65` is the run ID
      - `C6UU2ANXX` is the flowcell ID
      - `1` is the flowcell lane
      - `1101` is the tile number within the lane
      - `2046` x-coord within the tile
      - `1852` y-coord within the tile (followed by a space)
      - `1` is the member of a pair (if paired end)
      - `N` is filtered Y/N?
      - `0` "0 when none of the control bits are on, otherwise it is an even number"??
      - `GAGATTCCTATAGCC` is the "index sequence"
    - to fix this it they [recommend](http://gatkforums.broadinstitute.org/gatk/discussion/2909) using picardTools' `AddOrReplaceReadGroups`... added `fix_readgroup_array.qsub` script to run between the deduplication and the indel realignment steps.
      - but *of course* the scratch space on the HPCC is doing the I/O error thing where it can't glob files correctly... delays.
      - I sent a list of the node IDs for the failed jobs to the iCER team (using `grep Owner Dedup_* | cut -d '@' -f 2 | uniq > failed.job.nodes`)
    -

Week of August 22nd

  - things!
  -
  - [ABySS](http://computing.bio.cam.ac.uk/local/doc/abyss.html#scaffolding) is now able to do super-scaffolding with long reads. I'm going to try it with our reference genome and the Bionano map:
    - `abyss-pe np=8 k=64 name=new_P_kings lib='pe1' pe1='../../P_kings_genome/supercontigs.fasta' long='../super_scaffold/Para_king_2015_013_20_40_15_90_3_superscaffold.fasta_contig.fasta'`
    - then to compare the old vs. new assemblies using BUSCO...
    - I installed BUSCO using `brew` and downloaded the [vertebrate database](http://busco.ezlab.org/files/vertebrata_buscos.tar.gz)
    - to run BUSCO, need to load `BLAST+`, `augustus`, `HMMER` modules
    - `python ~/BUSCO_v1.22/BUSCO_v1.1b.py -o P_k_genome_test -in /mnt/scratch/pitchers/eFISH/P_kings_genome/supercontigs.fasta -l /mnt/home/pitchers/vertebrata/ –m genome`
    - `python ~/BUSCO_v1.22/BUSCO_v1.1b.py -o P_k_genome_test -in /mnt/scratch/pitchers/eFISH/bionano/Abyss_scaffolded_assembly/new_P_kings-1.fa -l /mnt/home/pitchers/vertebrata/ –m genome`
