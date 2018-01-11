# 2015_genomic_data folder
---
*opened by WRP on Wednesday, 12 August 2015*

Wednesday, 12 August 2015

 - all contents of `20150707_DNASeq_PE`, `20150721_DASeq_PE` & `BioNano` directories uploaded from the lab NAS drive.
 - md5sum checked on all `...fasta.tar.gz` files
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

  - Found a *bizarre* bug whereby the colour-coding in my shell settings was causing `bwa mem` to add colour-code escaped characters to the outputed `..sam` files (but only sometimes because it wants me to doubt my sanity). I have edited my script to pass the call to named variables in the `bwa mem` command through `sed` to trip out any color-codes that make it that far.
  - Now re-re-running.


Wednesday, 7 October 2015

  - FTP-ed the replacement `.bam` files from the replacement 'Original_data' set to the lab NAS drive for cold storage.
  - pipe has run unsupervised through trimming, alignment, indexing, sorting, deduping, indel realignment and base recalibration! ...but some files seem to have been dropped on the way.
  - I wrote `find_missing_files.sh` to list what's missing... 3 fastq files are weird in a way that s Picard. Investigating...


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

  - Xmas  week -- conveniently enough the iCER staff have partially shut down the HPC to fix some problem with scratch... inconveniently, that is where all our data lives.

---

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
  - worked out how to get `Dancer` running...
  - `Dancer` run appears not to have finished cleanly... the problem seems to be in the `.bam` files. Investigating...


Week of 8th – 12th February

  - so Dancer seems to fail due to a missing end of file (EOF) marker on one of the bam files.
      - I'm backtracking along my pipeline to see when this got dropped...
      - I now suspect this may the result of PICARD & SAMTools disagreeing about what a `.bam` file *ought* to look like...
      - I'm running a test of Dancer on a `.bam` file that SAMTools hasn't touched.
  - Dancer repeatedly fails to finish the run, still blaming EOF marker problems...
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
  - trying plan E re: Dancer `.bam` EOF problem...


Week of 29th Feb. – 4th March

  - Dancer still failing... got as far as sample BAVA_6623. I am going to back-track to the recalibration step to (hopefully) find where the file got truncated/lost its EOF.
    - there seem to be quite a few of the original pool of 768 bam files that left the recalibration step with a malformed EOF (but not severely enough to  the variant caller)... checking them all.
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
  - We think that a proportion of the apparently missing-data codes in the multi-individual vcf files may in fact represent homozygous loci for the reference allele... [seqanswers](http://seqanswers.com/forums/archive/index.php/t-28325.html) suggests that this is a known 'feature' of GATK.
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
  - THROUGH: the job appears to have run correctly in the small hours of the morning Wednesday 13th!
    - running PLINK pipeline: "File contains 31896675 entries and 63 individuals"
      - re-running the analysis
      - `vcftools --vcf all_individuals_12_07_16.out.vcf --out "tophits_all_individuals_12_07_16.out" --positions ${DIR}top_hit_variants --recode`
      - `vcftools --vcf tophits_all_individuals_12_07_16.out.recode.vcf --out "tophits_all_individuals_12_07_16" --freq`
      - output the genotype-by-individual (GT) format using `vcftools --vcf tophits_all_individuals_12_07_16.out.recode.vcf --extract-FORMAT-info GT --out all_individuals_12_07_16`
      - reformat this output using `cat all_individuals_12_07_16.GT.FORMAT | sed s/1\\/1/2/g | sed s/0\\/0/0/g | sed s/0\\/1/1/g | sed s/\\.\\/\\./NA/g > all_individuals_12_07_16.GT.FORMAT.recode`
      - `grep "Scaffold81" all_individuals_12_07_16.20x4.out.vcf > scaf81.list`
      - `for i in {0..4667} ; do grep "Scaffold${i}[[:space:]]" all_individuals_12_07_16.out.vcf | tail -1 > scaf.test ; done`


  - 15 July

    - vcf/PLINK testing: association analysis seems to have some weirdness – SNPs called that are at coords > length of scaffold...
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


Week of 29th July - 5th August

  - The goal of this week is to make sense of the lack of overlap among the lists of 'top hits' from the four different PLINK association runs...
  - bootstrap things...


Week of 8th August

  - is the problem with the sam files?
    - `java -jar $PICARD/ValidateSamFile.jar I=APA_6675_GAGATTCC-TATAGCC_L001_R1_pe.fq MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6724_GAATTCGT-TAATCTT_L008_R1_pe.fq MODE=SUMMARY` -> "no errors found"
    - `java -jar $PICARD/ValidateSamFile.jar I=MOV_6722_GAATTCGT-GGCTCTG_L003_R1_pe.fq MODE=SUMMARY` -> "no errors found"
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


Week of August 22nd

  - [ABySS](http://computing.bio.cam.ac.uk/local/doc/abyss.html#scaffolding) is now able to do super-scaffolding with long reads. I'm going to try it with our reference genome and the Bionano map:
    - `abyss-pe np=8 k=64 name=new_P_kings lib='pe1' pe1='../../P_kings_genome/supercontigs.fasta' long='../super_scaffold/Para_king_2015_013_20_40_15_90_3_superscaffold.fasta_contig.fasta'`
    - then to compare the old vs. new assemblies using BUSCO...
    - I installed BUSCO using `brew` and downloaded the [vertebrate database](http://busco.ezlab.org/files/vertebrata_buscos.tar.gz)
    - to run BUSCO, need to load `BLAST+`, `augustus`, `HMMER` modules
    - `python ~/BUSCO_v1.22/BUSCO_v1.1b.py -o P_k_genome_test -in /mnt/scratch/pitchers/eFISH/P_kings_genome/supercontigs.fasta -l /mnt/home/pitchers/vertebrata/ –m genome`
    - `python ~/BUSCO_v1.22/BUSCO_v1.1b.py -o P_k_genome_test -in /mnt/scratch/pitchers/eFISH/bionano/Abyss_scaffolded_assembly/new_P_kings-1.fa -l /mnt/home/pitchers/vertebrata/ –m genome`


Week of August 29th

  - check completeness of ..g.vcf files...
    - at the shell `for i in *all_libraries.bam.g.vcf ; do echo ${i} >> GVCF_counts_1Sept.txt ; tail -1 ${i} | cut -f 1 >> GVCF_counts_1Sept.txt ; done`
    - ...then in vim `:%s/vcf\n/vcf\t/g`
    - ...then in **R** `tbl_df(read.table( "GVCF_counts_1Sept.txt", header=FALSE )) %>% transmute( ID=V1, Scaf=as.integer( sub( "Scaffold", "", V2 )) ) %>% arrange( Scaf ) %>% write.table( "truncatedGvcfs.txt", quote=FALSE, row.names=FALSE )`
  - Make V.3 top-SNPs subset `.vcf` for JG:
    - output list of 'hits' from `Old_vs_New.Rmd` -> `top_hit_variants3`
    - `vcftools --vcf all_individuals_29_08_16.noloc.out.vcf --out "tophits_all_individuals_29_08_16.out" --positions top_hit_variants3 --recode`
  -
    - make vcf subset *without* APA & BAM `bcftools view -Ov -s ^APA_6675,APA_6676,APA_6677,APA_6678,APA_6679,APA_6680,APA_6681,APA_6682,APA_6683,APA_6684,APA_6685,APA_6737,BAM_6494,BAM_6496,BAM_6497,BAM_6498,BAM_6499,BAM_6500,BAM_6501,BAM_6502,BAM_6597,BAM_6598,BAM_6599,BAM_6602,BAM_6603,BAM_6604,BAM_6605 all_individuals_29_08_16.noloc.out.vcf > without_APA_and_BAM_29_08_16.noloc.out.vcf`
    - make vcf subset with *only* APA & BAM `bcftools view -Ov -s APA_6675,APA_6676,APA_6677,APA_6678,APA_6679,APA_6680,APA_6681,APA_6682,APA_6683,APA_6684,APA_6685,APA_6737,BAM_6494,BAM_6496,BAM_6497,BAM_6498,BAM_6499,BAM_6500,BAM_6501,BAM_6502,BAM_6597,BAM_6598,BAM_6599,BAM_6602,BAM_6603,BAM_6604,BAM_6605 all_individuals_29_08_16.noloc.out.vcf > only_APA_and_BAM_29_08_16.noloc.out.vcf`
    - make vcf subset *without* COB `bcftools view -Ov -s ^COB_4004,COB_4006,COB_4018,COB_4019,COB_4027,COB_4029 all_individuals_29_08_16.noloc.out.vcf > without_COB_29_08_16.noloc.out.vcf`
  - run `GATK/GenotypeGVCFs` *in parallel* on all 63 fish for the sake of comparison


Week of 6th September

  - we now have 2 GVCF-derived all-fish.vcf files – the goal is to see how they differ:
    - `wc`: 36359475  2617544446 23191271366 `all_individuals_post_GVCF_merged_04_09_2016.vcf`
      - 36354704 non-header rows
    - `wc`: 31671518  2280011585 56405095271 `all_individuals_29_08_16.noloc.out.vcf`
      - 31666749 non-header rows
    - running `all_individuals_post_GVCF_merged_04_09_2016.vcf` through the PLINK pipeline... `all_individuals_post_GVCF_merged_04_09_2016_maf10_geno50.assoc.fisher` generated.
    - continuing comparisons in [`Old_vs_New.Rmd`](./Old_vs_New.Rmd)...
    - running a vcf-diff with `vcftools --vcf all_individuals_29_08_16.noloc.out.vcf --diff all_individuals_post_GVCF_merged_04_09_2016.vcf`
      - `out.diff.indv_in_files` list all individuals as included in both files.
      - `out.diff.sites_in_files` is long and complex...
        - how many variants differ? `wc -l out.diff.sites_in_files` – 36075225 - NOPE! all sites are listed!
        - `cat out.diff.sites_in_files | cut -f 3 | grep -c B` > 31318405 sites present in both files
        - `cat out.diff.sites_in_files | cut -f 3 | grep -c 1` > 348343 sites present only in v.3
        - `cat out.diff.sites_in_files | cut -f 3 | grep -c 2` > 4408476 sites present only in v.4
        - make a file of *just* unshared sites `grep --color='never' -v B out.diff.sites_in_files > out.diff.sites_not_in_both_files`
        - `wc -l out.diff.sites_not_in_both_files` – 4756820 (as it should be)
        - latest scaf included seems to be Scaffold999... because string-ordering! (4534 unique scaf no.s present)
        - `out.diff.sites_not_in_both_files` is small enough to be easier to deal with in **R** -> [`VCF_version_3vs4.Rmd`](./VCF_version_3vs4.Rmd) for details...
          - need a table of scaffold lengths: these can be extracted from `supercontigs.fasta.fai`
          - I'm becoming suspicious after seeing some plots of the distribution of unshared SNPs by scaffold... which scaffolds appear in version3? -> `grep -v ^#  all_individuals_29_08_16.noloc.out.vcf | cut -f 1 | uniq` gives me a list.
        - wondering which sites appear in both files, but with different alternate allele... `awk '{ if ( $3 == "B" && $5 != $6 )  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' out.diff.sites_in_files >> out.diff.sites_varible_between_files`
  - also added `pop_level_vcf_discovery.qsub` script to call genotypes on a per-population basis...


Week of 12-16th September

  - ran a few tests of `pop_level_vcf_discovery.qsub` over the weekend. These are going to be long-running jobs. I have them queued on the HPC to run for 24hrs – currently due to start on the 18/19th. 1st job is also running on Shockly as a full-scale test. (these jobs started early as it happened: 15 September 2016)
    - I'm going to run the (incomplete) `pop_XXX_all_fish_10_09_16.vcf` files from these tests through `GATK/CombineVariants` -- I can probably learn something by making comparisons with the other vcf versions over just the first few scaffolds.
      - running the vcf-diff command on the 2 new combinations: `vcftools --vcf all_individuals_29_08_16.noloc.out.vcf --diff all_fish_genocalled_as_pops_12_09_16.vcf` & `vcftools --vcf all_individuals_post_GVCF_merged_04_09_2016.vcf --diff all_fish_genocalled_as_pops_12_09_16.vcf`
  - posted on GATK forums about our pipeline indecision – got an answer really quickly:

  > "There are a few technical obstacles that arise if you don't have a prior collection of known variants, as is often the case for non-model organisms. Ie the recalibration tools require known variants (though this can be bootstrapped for BQSR). For variant filtering, you'll need to use hard filters as the machine learning tools won't work on your data.

  > Besides that there's not really anything that is assumed that would be human specific. The GVCF workflow does assume that variants are more likely to be real if they're observed in more than one sample -- but given your experimental design I don't think that would be a problem. You're looking for variants that are shared within sub populations, right?

  > The low coverage is probably the biggest weakness; you may want want to do some tests comparing sensitivity between traditional multi sample calling (giving all samples simultaneously to the HaplotypeCaller) versus applying the two-step GVCF workflow."

  - Following the above advice entails a slight adjustment from pipeline version 1.; firstly to ensure that the readgroups are available for BSQR and then to funnel all the fish into one HaplotypeCaller run... making adjustments.
    - re-running the pipeline from the top to move the `fix_readgroup_array` step up.
      - an odd assortment of ~10 `fix_readgroup_array` jobs failed, followed by 26 `deduplication_array` jobs...?
      - I'm running a few of these interactively to see if I can ID the problem(s)... as far as I can tell the initial failures appear to be I/O errors on the older nodes; I cannot replicate them. This has become so tedious that I'm going to specify that all my jobs should *only* run on the 'intel16' nodes. <grumblegrumble>
      - the `fix_readgroup_array` jobs that crash obviously propagate their error down to `deduplication_array`, but it seems that the 'fixed' readgroups may be *causing* an error in some cases too... `MOV_6725_GAATTCGT-CAGGACG_L007_R1_pe.fq` is my test case:
          - `..fq` passes `$PICARD/ValidateSamFile` test...
          - `..aligned.rg.sam` fails with 'MISMATCH_READ_LENGTH_AND_QUALS_LENGTH'...
          - but re-running the `$PICARD/AddOrReplaceReadGroups` command manually yields a `..aligned.rg.sam` that passes (!?) This is odd. Testing again with `BAM_6498_ATTACTCG-GGCTCTG_L002_R1_pe.fq`...
          - these files produce `..aligned.rg.sam` outputs that pass, but contain only a header!
          - the header s (at least in this test case) because the `$PICARD/AddOrReplaceReadGroups` tool meets an empty/missing read with an incomplete header... the plot thickens..


Week of 19-23rd September

  - first priority – readgroup problem...
    - test 1: can I have `bwa-mem` fix the readgroups?
      - wrote `Alt-Align.qsub` script: output `BAM_6498_ATTACTCG-GGCTCTG_L002_R1_pe.fq.sorted.dedup.taco`
        - this seems to work nicely, so going with this for the time being...
    - test 2: can I use `samtools` to fix the readgroups instead of `PICARD`?
      - if we upgrade to `SAMTools/1.3.1` we can use `addreplacerg`: output `BAM_6498_ATTACTCG-GGCTCTG_L002_R1_pe.aligned.dedup.bam`...?
    - `samtools view -b MOV_6725_all_libraries.bam Scaffold0 > MOV_6725_all_libraries_scaf0.bam`
      - need to do `samtools index` before split, then `$PICARD/BuildBamIndex` after
        - Scaffold0, 1 fish; `HaplotypeCaller` takes only 9mins
      - testing with `all_bam_merged_recal_24_05_2016.bam` (all 63 fishes, 592GB)
        - samtools indexing takes 2hrs
        - splitting takes only ~2mins
        - `HaplotypeCaller` takes ~7mins for 1kb, ~10mins for 10kb, ~21mins for 100kb, 1Mb for ~2.5hrs, ~4hrs for 1.5Mb
          - memory-wise: 1Mb allocated 300GB – actually used ~=50GB
        - PICARD indexing takes ~7mins
    - OK, so we want to aim for bam subsets of ~1.25Mb, full genome is 799,426,128bp, divided by 1.5Mb is 639.54, so let's go with as nice round number of 650 chunks...
      - no. lines in full bam: 1807387006 / 650 ~= 2780596...
      - the plan should be something like: `samtools view all_fish.bam | split -d -l 2780596 - BamChunks/all_fish_`
  - OK. change of plans... using `split` is much faster than using `samtools view` for the files I was testing on, but to chunk a merged `..bam` file covering the full genome and all 63 fish requires it to run for ~12 hrs. I made 3 attempts but in every case the job encountered an I/O error long before nearing completion... `samtools view` it is!
    - the way I've set this up is to build a script - `write_scaf_indices.sh` - that generates a list - `indices.list` - of scaffolds and coordinates that can be read into `vcf_disco_chunk_array.qsub`.
    - `write_scaf_indices.sh` has a setting for the max. length of chunk, and writes out scaf-subsets where scaf-length > chunk-length, and doubles-up or octuples-up scafs where scaf-length << chunk-length...

====

Week of 26-30th September

  - HPC nodes are going up & down like yo-yos today... I am firefighting numerous small problems.
'malformed'
  - missing sam files: `APA_6682_ATTCAGAA-TATAGCC_L003_R1_pe.trimmed.fq`, `BAM_6597_TCCGGAGA-TATAGCC_L001_R1_pe.trimmed.fq` & `MOV_6724_GAATTCGT-TAATCTT_L007_R1_pe.trimmed.fq`
    - rerunning the appropriate instance of the alignment array manually, then testing with `$PICARD/ValidateSamFile`. No errors found.
  - then re-running the deduplication... huh. somehow, just having PICARD `MarkDuplicates` seems to cause these files to fail `ValidateSamFile`... this is *weird*...
    - these files seem to be being changed when they didn't ought to be... can I make them read-only?
      - I can make them read-only. This does not solve the problem...
    - I'm going to try the preceding operation – `bwa-mem` – with the `-T 40` flag to prevent it outputting lower-quality alignments to see if this allows PICARD to work with them...
    - this seems to be working... I have the `APA_6682..` file as far as the realignment step without error...
    - `BAM_6597..` is OK so far...
    - `MOV_6724..` still fails test after deduplication (!)
    - ...
  - OK. Change of direction...


Week of 3-7th October

  - splitting the pipeline to output version 5.1 and version 5.2 VCF files
    - with reference to http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode we have decided that we need to know what difference it makes to the output whether the GATK's haplotype calling and genotyping steps are run separately or together...
    - both versions 5 will run as arrays so that the GATK walkers can deal with small subsets of sequence at a time
    - testing suggests that GATK can get through ~1.5MBp of sequence takes in 4hrs, so I intend to have the array work on 1MBp slices to keep the jobs short and easy on the scheduler
    - wrote `write_scaf_indices.sh` script to calculate the 'intervals' needed to pass to the walkers
  - **Excellent News** – we apparently have a research group scratch space!
    - I have moved the end of the pipeline to dump its effluent into the new space...

Week of 10-14th October

  - 124820 V5-2 vcf files
  - find duplicate files with ``ls all_fish_*slice_*.vcf > foo && for i in `cat foo` ; do slice=`echo ${i} | cut -d'_' -f 7`;  if [ `ls *slice_${slice} | wc -l` -gt 1 ]; then echo ${i}; fi; done``
  - `java -Xmx10g -cp /opt/software/GATK/3.5.0/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${ref} ${samples} -out all_fish_version_5-2_HPC.vcf`
  - **Stopping points**
    - V5-2_chunks/IVI_4897/IVI_4897_10_2016_slice_598.g.vcf x
    - V5-2_chunks/MOV_6725/MOV_6725_10_2016_slice_1211.g.vcf x
    - V5-2_chunks/MOV_6724/MOV_6724_10_2016_slice_1397.g.vcf x
    - V5-2_chunks/BAVA_6621/BAVA_6621_10_2016_slice_1423.g.vcf x
    - V5-2_chunks/MOV_6724/MOV_6724_10_2016_slice_1382.g.vcf x
    - V5-2_chunks/BAVA_6627/BAVA_6627_10_2016_slice_1367.g.vcf x
    - V5-2_chunks/COB_4029/COB_4029_10_2016_slice_1375.g.vcf x
    - V5-2_chunks/MOV_6716/MOV_6716_10_2016_slice_1404.g.vcf x
    - V5-2_chunks/MOV_6718/MOV_6718_10_2016_slice_1365.g.vcf x
    - V5-2_chunks/MOV_6717/MOV_6717_10_2016_slice_1374.g.vcf x
    - V5-2_chunks/MOV_6722/MOV_6722_10_2016_slice_1432.g.vcf x
    - V5-2_chunks/MOV_6720/MOV_6720_10_2016_slice_1202.g.vcf x


Week of 17th-21st October

setdiff( expected, taco$V1 )
698  798  898  998 1098 1198 1298 1398 1498 |1598
1302 1402  1502 1602 1702 |1802
1311 1411 1511 1611 1711 |1811
1465 1475 1565 1665 |1765
1467 1567 1667 |1767
1474 1574 |1674
1482 1582 1682 |1782
1497 1597 1697 |1798
1504 1604 1704 1804
1523 1623 1723 1823
1532 1632 1732 1832
1575 1675 1775
1797
1698
1774

dir=/mnt/ls15/scratch/groups/efish/WILL/
ref=/mnt/scratch/pitchers/eFISH/P_kings_genome/supercontigs.fasta

java -Xmx60g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T VariantsToBinaryPed \
             -R ${ref}  -V ${input_data} \
             -m ${meta} -mgq 0 \
             -bed ${input_data}.bed \
             -bim ${input_data}.bim \
             -fam ${input_data}.fam

GATK docs lie about this function!

Wednesday, 19 October 2016 Scaffold133:1400019 malformed line crashed job. remaking vcf slice 635
Thursday, 20 October 2016 Scaffold154:1400004 malformed line crashed job. remaking vcf slice 698
  JG points out that this behaviour is very I/O-error-reminiscent...

Thursday, 20 October 2016 Scaffold167:15 malformed line crashing job. remaking vcf slice 735


vcftools `vcf-sort -c infile.vcf > outfile.vcf`

comparing `all_fish_version_5-1_shockly.vcf` & `all_fish_version_5-1_HPC.vcf`

https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php

> Output
> Either a VCF or gVCF file with raw, unfiltered SNP and indel calls. Regular VCFs must be filtered either by variant recalibration (best) or hard-filtering before use in downstream analyses. If using the reference-confidence model workflow for cohort analysis, the output is a GVCF file that must first be run through GenotypeGVCFs and then filtering before further analysis.


Week of 24th-28th October

  -



Week of 31st Oct. - 4th Nov.

  - Implemented extra error-checking `if` statements in the VCF_5.2 as suggested by Andrew Keen...
  - latest 5.2 discovery array seems to be running well, but after **3** complete runs, there are still **293** `..g.vcf` files missing.
    - I'm re-running a couple of these interactively to see if I can find out what's going on...
      - 1 thing to note: the `readarray` command seems to not be stripping whitespace. Added the `-t` flag to force this behaviour...
      - 2ndarily, one of the java error-checks suggested by AK seems to be ing my error-checking `if` statements... I'm commenting it out for the nonce
      - so there also appear to be a pair of version-specific GATK issues to work around:
        - the "v3.6-0-g89b7209" version cannot handle >20 variants at a locus. This appears to be the cause of some of the failed `HaplotypeCaller` jobs, specifically e.g. slice 1807 which failed for all fish. Presumably these regions are just particularly variant-rich
        - the "vnightly-2016-11-01-gaca5d7b" nightly build handles many multiple alleles handily, but does not include (strangely) the ability to use the '-gvcf' flag for `ValidateVariants`.
  - after kludging through the last few `..g.vcf` files, qsub-ed the vcf genotyping array...
    - new problems!
      - 1534 of 2765 jobs have `job.o..` reports
      - some reports claim that GATK failed to detect the `#CHROM` header line in the `..g.vcf` files – however `grep` has no problem finding these lines (in `..g.vcf`s that have already *passed* `ValidateVariants` look you!)
      - 12 of the output `..vcf` file are empty.
    - interactively running one of the failed jobs with the stable-release (3.6.0) version of GATK completes without problems. I have adjusted the vcf genotyping script and re-submitted.
  - 3.6.0 version re-run of vcf genotyping leaves 100 files missing...
    - errors all seem to be traceable to an *empty* (i.e. 0 byte) `..g.vcf` file in each case
    - re-running the relevant vcf_disco_chunk_5-2.. job seems to fix matters – after which I'm re-submitting the vcf_genotype_chunk_5-2.. job
    - ...but 8/100 of the broken jobs must be missing >1 `..g.vcf` file, because I need to re-re-run 8 jobs
    - after which I ran `merge_all_vcfs`, which was stopped by an error in `all_fish_slice_1699.vcf`.
      - I fixed this by deleting the 63 `..slice_1699.g.vcf` files, re building them and re-running `vcf_genotype_chunk_5-2` on that slice.
    - after which I re-ran `merge_all_vcfs`, which was stopped by an error in `all_fish_slice_2749.vcf`..
      - again I deleted the offending `..vcf`, rebuilt the constituent `..g.vcf`s, and re-ran the genotyping step, before re-submitting the merge script
    - `merge_all_vcfs`, was then stopped by an error in `all_fish_slice_2756.vcf`... <howls of impotent rage>
      - *weirdness*; `all_fish_slice_2756.vcf` *passes* validation, even without the `validationTypeToExclude` flags. GATK only seems to be able to find a problem with it when running `CatVariants`(!)
      - as a backup plan; `scp`-ing all the `..vcf`s to Shockly again

Week of 5th – 11th Nov.

  - returning to `CatVariants` problem...
    - 2756 is the point of failure again... same error message:
        >  Line 4955: there aren't enough columns for line  (we expected 9 tokens, and saw 1 ), for input source: /mnt/ls15/scratch/groups/efish/WILL/V5-2_chunks/all_fish_slice_2756.vcf

    - this is the last line in the file and is possibly truncated...?
    - same point of failure occurs with `CatVariants` on Shockly.
    - make new idx and diff the 2 – check age of indices?
  - tempdir issue

  - 5.1 version: 9032114 rows – 5.2 version: 7902104 ... 1130010 fewer variants in `g.vcf` version
    - 7360384 variants are present in both versions: 81% of v5.1 & 93% of v5.2
    - R^2 between v5.1 & v5.2 p-values = 0.96, OR = 0.93
  - rerunning `calc_Fst_array` on both the new-shiny `..vcf`s...


Week of 14 – 18th Nov.

  - PLINK association analyses
    - as previously, with a filter at MAF of 10% and a genotyping rate of 50%, with both v5.1 and v5.2 versions
    - without an MAF filter, and with a genotyping rate filter of 50%, with both v5.1 and v5.2 versions
  - comparisons made between v5.1 & v5.2 in .


Week of 21 – 25th Nov.

  - HPC is giving me error messages about having no free space (won't even let me delete old files), so I'm scp-ing to shockly today (Monday)... [*Update*](http://icer.msu.edu/service-status)
    - subsetting out the APA & BAM individuals from the v5.1 vcf file with: `java -Xmx30g -jar /home/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R ../individual_bams/supercontigs.fasta -V all_fish_version_5-1_HPC.vcf -o APA_BAM_only_Nov_HPC.vcf -sn APA_6675 -sn APA_6676 -sn APA_6677 -sn APA_6678 -sn APA_6679 -sn APA_6680 -sn APA_6681 -sn APA_6682 -sn APA_6683 -sn APA_6684 -sn APA_6685 -sn APA_6737 -sn BAM_6494 -sn BAM_6496 -sn BAM_6497 -sn BAM_6498 -sn BAM_6499 -sn BAM_6500 -sn BAM_6501 -sn BAM_6502 -sn BAM_6597 -sn BAM_6598 -sn BAM_6599 -sn BAM_6602 -sn BAM_6603 -sn BAM_6604 -sn BAM_6605`
    - subsetting out the IVI & MOV individuals from the v5.1 vcf file with: `java -Xmx30g -jar /home/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R ../individual_bams/supercontigs.fasta -V all_fish_version_5-1_HPC.vcf -o IVI_MOV_only_Nov_HPC.vcf -sn IVI_3923 -sn IVI_4816 -sn IVI_4832 -sn IVI_4834 -sn IVI_4893 -sn IVI_4894 -sn IVI_4895 -sn IVI_4896 -sn IVI_4897 -sn IVI_4921 -sn IVI_4925 -sn MOV_6716 -sn MOV_6717 -sn MOV_6718 -sn MOV_6719 -sn MOV_6720 -sn MOV_6721 -sn MOV_6722 -sn MOV_6723 -sn MOV_6724 -sn MOV_6725`
    - subsetting out the IVI & APA individuals from the v5.1 vcf file with: `java -Xmx30g -jar /home/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R ../individual_bams/supercontigs.fasta -V all_fish_version_5-1_HPC.vcf -o APA_IVI_only_Nov_HPC.vcf -sn APA_6675 -sn APA_6676 -sn APA_6677 -sn APA_6678 -sn APA_6679 -sn APA_6680 -sn APA_6681 -sn APA_6682 -sn APA_6683 -sn APA_6684 -sn APA_6685 -sn APA_6737 -sn IVI_3923 -sn IVI_4816 -sn IVI_4832 -sn IVI_4834 -sn IVI_4893 -sn IVI_4894 -sn IVI_4895 -sn IVI_4896 -sn IVI_4897 -sn IVI_4921 -sn IVI_4925`
    - subsetting out the BAM & BAVA individuals from the v5.1 vcf file with: `java -Xmx30g -jar /home/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T SelectVariants -R ../individual_bams/supercontigs.fasta -V all_fish_version_5-1_HPC.vcf -o BAVA_BAM_only_Nov_HPC.vcf -sn BAM_6494 -sn BAM_6496 -sn BAM_6497 -sn BAM_6498 -sn BAM_6499 -sn BAM_6500 -sn BAM_6501 -sn BAM_6502 -sn BAM_6597 -sn BAM_6598 -sn BAM_6599 -sn BAM_6602 -sn BAM_6603 -sn BAM_6604 -sn BAM_6605 -sn BAVA_6619 -sn BAVA_6620 -sn BAVA_6621 -sn BAVA_6622 -sn BAVA_6623 -sn BAVA_6624 -sn BAVA_6625 -sn BAVA_6626 -sn BAVA_6627`
  - re-ran the PLINK association on the subset vcfs. IVI vs MOV is a comparison that we can't make with this data as all IVI & MOV fish are triphasic...
  - Updating the installation or R on Shockly... some weirdness occurring...


Week 28th Nov. – 2nd Dec.

  - Discussion with JG about best way to proceed vis-a-vis understanding our effect sizes... we
  - new PLINK command:
    - `plink --bfile ${input_data} --allow-no-sex --allow-extra-chr --logistic --ci 0.95 --pfilter 1 --out ${input_data} --covar ${pheno} --covar-number 3`
      - `--bfile`, ` --allow-no-sex`, `--allow-extra-chr`, & `--out` are as before...
      - `--logistic` fits a logistic regression model instead of the Fisher's exact test
      - `--ci 0.95` outputs 95% confidence intervals on the estimates of OR
      - `--pfilter 1` is the command to filter by p-value, with the option '1' this just removes NAs
      - `--covar ...` specifies a path to a file with covariate data and `--covar-number` gives the column number of the covariate to include... the covariate-containing file needs to look like this:
      >FID	IID	POP
      >APA	APA_6675	APA
      >...

    - `plink --bfile ${input_data} --allow-no-sex --allow-extra-chr --model --fisher --ci 0.95 --pfilter 1 --out ${input_data} --test-missing`

      - choice of test must be appropriate for small our sample size
      - use --model --fisher and add permutation? ...and add --test-missing? & counts?
      then after – get permutation p-values
  - plan:– sort by exclusivity in `.model` output and then find/compare P-values with...
    - `hypoth=( GENO ALLELIC DOM REC TREND )`, then ``for i in ${hypoth[@]} ; do grep ${i} ${input_data} >> `basename ${input_data} .model`.${i}.model ; done``
    - checking the coding in genotype counts:
      - `..GENO..` file has homo/het/homo in the order ref.ref : ref.alt : alt.alt
      - `..ALLELIC..` file has n.ref / n.alt
      - `..DOM..` file has n.homo.ref+n.het / n.homo.alt
      - `..REC..` file has n.homo.ref / n.het+n.homo.alt
      - `..TREND..` file has n.ref / n.alt
      -
      - `awk -F "[\t ]+" ' NR>1 { split( $6, tri, "/") ; split( $7, bi, "/") ; print $0 "\t" tri[1] "\t" tri[2] "\t" tri[3] "\t" bi[1] "\t" bi[2] "\t" bi[3] "\t"   }' all_fish_version_5-1_HPC.GENO.model >> all_fish_version_5-1_HPC.GENO.model.reformatted` followed by `echo -e "CHR SNP A1 A2 TEST AFF UNAFF P p_Hr p_Het p_Ha np_Hr np_Het np_Ha" | cat - all_fish_version_5-1_HPC.GENO.model.reformatted | sponge all_fish_version_5-1_HPC.GENO.model.reformatted` to fix the header row and `awk '{$1=$1}1' OFS="," all_fish_version_5-1_HPC.GENO.model.reformatted > all_fish_version_5-1_HPC.GENO.model.reformatted.csv` to strip out the weird multi-space delimiters.
  - I'm going to start with the `..GENO..` file and sort to find max & min homozygotes...


Week of Dec. 5-10th

  - re-making presentation-quality Fst plot for JG's talk...

  - JG confirms mis-classification of 3 BAM fish in the phenotype 'source' file  `/mnt/research/efish/2015_genomic_data/specimens_for_genome_reseq.txt`... I manually edited this file to switch the phenotype coding from 'P0 present' to 'P0 absent' for BAM_6494, BAM_6497 & BAM_6500
  -

`assoc <- mutate( taco, CHR=factor( gsub( "Var-", "", gsub( "\\-\\d+", "", taco$SNP ))), SNP=as.integer(gsub( "Var\\-Scaffold\\d+\\-", "", taco$SNP )))`
  - for the `ALLELIC` file:


 –––  for the Holidays: 12th Dec. - 9th Jan. –––


Week of 16-20th Jan.

  - can we test biallelic-esque ref vs. not-ref as a way of shoehorning multi-allelic-ness
  - `bcftools norm all_fish_version_5-1_HPC.vcf.gz -o all_fish_version_5-1_HPC.NORMALIZED.vcf.gz` ?
  - `plink --file all_fish_version_5-1_HPC --freq --out all_fish_version_5-1_HPC.freq`

  - Analyses to compare after incorporating 'ultimate phenotypes':
    - all fish, all hypotheses, real ref alleles
    - APA & BAM only, all hypotheses, real ref alleles
    - not-COB, all hypotheses, real ref alleles

  - need to get the group scratch directory tidied up...
    - keeping: individual BAMs, output vcfs
    - tossing: slice-wise vcfs, old figures, check-if-done dummy files/folders

  - testing the effect of population as a covariate –
  - testing for agreement between PLINK and R:fdrtool...?
    - `--adjust` in PLINK is only available for trend-only models...


Week of 30th Jan. – 3rd Feb.

  - Meeting w/JG:
    - signal:noise ratio between candidates & non-candidates... p-val distribution
    - models we want (as per Asana), a direrctory containing:
      1. Fexact-test w/ genotypes spreadsheet top 250
      2. - 5. as above for other 4 hypothesis-wise tests
      6. comparison/ranking between p-vals from files 1-5
      7. analysis of p-val distribution
  - making a vcf subset using: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V all_fish_version_5-1_HPC.vcf -o without_COB_5-1_HPC.2017.vcf -sn APA_6675 -sn APA_6676 -sn APA_6677 -sn APA_6678 -sn APA_6679 -sn APA_6680 -sn APA_6681 -sn APA_6682 -sn APA_6683 -sn APA_6684 -sn APA_6685 -sn APA_6737 -sn BAM_6494 -sn BAM_6496 -sn BAM_6497 -sn BAM_6498 -sn BAM_6499 -sn BAM_6500 -sn BAM_6501 -sn BAM_6502 -sn BAM_6597 -sn BAM_6598 -sn BAM_6599 -sn BAM_6602 -sn BAM_6603 -sn BAM_6604 -sn BAM_6605 -sn BAVA_6619 -sn BAVA_6620 -sn BAVA_6621 -sn BAVA_6622 -sn BAVA_6623 -sn BAVA_6624 -sn BAVA_6625 -sn BAVA_6626 -sn BAVA_6627 -sn IVI_3923 -sn IVI_4816 -sn IVI_4832 -sn IVI_4834 -sn IVI_4893 -sn IVI_4894 -sn IVI_4895 -sn IVI_4896 -sn IVI_4897 -sn IVI_4921 -sn IVI_4925 -sn MOV_6716 -sn MOV_6717 -sn MOV_6718 -sn MOV_6719 -sn MOV_6720 -sn MOV_6721 -sn MOV_6722 -sn MOV_6723 -sn MOV_6724 -sn MOV_6725`


Week of 6-10th Feb.

  - quick-pic manhattan code: `assoc %>% filter( Scaf == "Scaffold858" ) %>% ggplot( aes( x=SNP, y=-log10( P )) ) + geom_point() + ylab( "-log(10) p value" ) + xlab( "BP coord" ) + ggtitle( "Scaffold 858" )`
  - maf settings and link to 'bands' in mahattan plots?
  - what happens if we bin *all* P0 as present?
  - what happens if we drop the mid-P0 individuals?
    - making a vcf subset using: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V all_fish_version_5-1_HPC.vcf -o no_small_P0_5-1_HPC_2017.vcf -sn APA_6675 -sn APA_6676 -sn APA_6677 -sn APA_6678 -sn APA_6679 -sn APA_6680 -sn APA_6682 -sn APA_6683 -sn APA_6684 -sn APA_6685 -sn APA_6737 -sn BAM_6496 -sn BAM_6497 -sn BAM_6498 -sn BAM_6500 -sn BAM_6501 -sn BAM_6502 -sn BAM_6597 -sn BAM_6598 -sn BAM_6599 -sn BAM_6602 -sn BAM_6603 -sn BAM_6605 -sn BAVA_6619 -sn BAVA_6620 -sn BAVA_6621 -sn BAVA_6622 -sn BAVA_6623 -sn BAVA_6626 -sn BAVA_6627 -sn IVI_3923 -sn IVI_4816 -sn IVI_4832 -sn IVI_4834 -sn IVI_4893 -sn IVI_4894 -sn IVI_4895 -sn IVI_4896 -sn IVI_4897 -sn IVI_4921 -sn IVI_4925 -sn MOV_6716 -sn MOV_6717 -sn MOV_6718 -sn MOV_6719 -sn MOV_6720 -sn MOV_6721 -sn MOV_6722 -sn MOV_6723 -sn MOV_6724 -sn MOV_6725`
  - can we filter out SNPs that are 'private' to populations?
  - NB: microsat scaffs - 4,8,22,60,77

  - vcftools to output *all* the stats:
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --out all_fish_version_5-1_HPC_out --freq`
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --out all_fish_version_5-1_HPC_out --site-mean-depth`
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --geno-depth`
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --out all_fish_version_5-1_HPC_out --site-quality`
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --out all_fish_version_5-1_HPC_out --SNPdensity 1000`
    - `vcftools --vcf all_fish_version_5-1_HPC.vcf --out all_fish_version_5-1_HPC_out --hardy`
      - outputs are `all_fish_version_5-1_HPC_out.frq`, `...ldpeth.mean`, `...ldepth`, `...gdepth`, `...lqual`, `...snpden` & `...hwe` respectively
    - `depth_...` script to read these outputs and generate some plots summaries...
  - 'QD' over 'QUAL'
  - follow hard-filter guidelines unless there's a good reason not to – then throw into PLINK again
    - `java -jar GenomeAnalysisTK.jar -R reference.fasta -T VariantsToTable -V file.vcf -F CHROM -F POS -F ID -F QUAL -F AC -o results.table`
    - `sed 's/0_//g' all_fish_version_5-1_HPC.assoc | sed s/":"/"_"/g > all_fish_version_5-1_HPC.assoc.in`
    - `awk '{ print $0 "\t" $1 "_" $2 }' all_fish_version_5-1_HPC.vcf.stats.table > all_fish_version_5-1_HPC.stats.in`
    - `sort -g -k 2 all_fish_version_5-1_HPC.assoc.in > all_fish_version_5-1_HPC.assoc.in.sorted`
    - `sort -n -k 15 all_fish_version_5-1_HPC.stats.in > all_fish_version_5-1_HPC.vcf.stats.in.sorted`
    - `sort -n -k 2 all_fish_version_5-1_HPC.assoc.in > all_fish_version_5-1_HPC.assoc.in.sorted`
    - `join -1 2 -2 15 all_fish_version_5-1_HPC.assoc.in.sorted all_fish_version_5-1_HPC.vcf.stats.in.sorted > tacopark`


Week of 13-17th

  - before filtering: 27871297 snps, & 5860821 indels --- after filtering: 26345714 snps, & 5572344 indels
  - filtered out 1525583 snps, & 288477 indels --- ~5.5% snps, & ~5% indels
  - loop over fish to produce quality density plots? ...wrote `loop_over_fishes.qsub` that calls `plot_depth_by_fish.R`
  - extract individual-fish quality scores from vcf with `vcftools --vcf all_fish_version_5-1_HPC.filtered.snps.vcf --extract-FORMAT-info GQ --out filtered_SNPs`
  - estimate of heterozygosity in *P. kingsleyae* 2.7e-2 – compare wild stickleback 1.43e−3, wild Medaka 1.5e-3...

    > Coefficients:                  Estimate Std. Error t value Pr(>|t|)
    > (Intercept)                   2.234e-01  2.125e-03 105.134  < 2e-16
    > as.numeric(assoc_means$Scaf) -5.889e-06  8.254e-07  -7.134 1.13e-12
    > Residual standard error: 0.07092 on 4456 degrees of freedom
    > Multiple R-squared:  0.01129,   Adjusted R-squared:  0.01107
    > F-statistic:  50.9 on 1 and 4456 DF,  p-value: 1.129e-12

  - after looking at GQ distributions for each fish, decided to run association analyses with cut-offs at ≤20 and ≤25
    - to that end, wrote `plink_prep_w_filter.qsub` script – using vcftools rather than GATK to do the filtering, because the GATK `VariantsToBinaryPed` tool skips from `...vcf` to `...bed` etc. and I need a `...ped` (i.e. flat text file) in order to be able to add in the phenotype values


Week of 20-25th

  - `vcftools --minGQ` does not appear to have filtered out any variants...
    - solution – `vcftools` can't filter on `INFO` fields that aren't present in *all* rows, and rows with non-variant genotypes don't get `GQ` scores... we fix this by implementing the `GATK` recommended hard-filter *before* the `GQ` filter
      - this ^ step takes ~2hrs, then renaming takes a few seconds, then `vcftools --plink` takes ~20mins, then awk to fix the `..map` rownames takes ~1min.... safe to leave it in 1 script w/4hr walltime.
    - genotypes do not seem to match phenoptyes in `filtered_output_GQ20.xlsx`...
      - the workflow goes: vcf -> filter -> plink files -> add pheno -> bed/bim -> plink assoc -> scrape ID's -> full join -> Excel
      - filtering step – individual rows that pass are unchanged – filtering removed 12175530 or 18562241 rows!?
        - indels/snps filtered separately
      - are genotypes being scraped from the right file? – makes no diff <– filtered or no. Good!
      - making plink files seems to work as it should
      - logically, 'F_A' and 'F_U' from the `PLINK assoc` output ought to total to the frequency of 'allele 1', which should be equal to one of the allele frequencies from the `vcftools ..frq` output, assuming that we have a handle on how these programs are behaving...
      - are the phenotypes being correctly assigned? - `ultimate_phenotype.txt` file is unchanged... - `cut -f 1,6 all_fish_version_5-1_HPC.filtered.snps_GQ20.ped > phenotest.txt` pulls out the phenotypes from the ped file (where plink actually sees them) – these phenotypes are correct
    - maybe I don't *need* to work out where the problem happens if I can route around it... wrote `vcf_to_bed.qsub` to use the `GATK VariantsToBinaryPed` tool instead of using vcftools->python->plink pipeline...
      - NB: the VariantsToBinaryPed tool [docs](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToBinaryPed.php) requires the `--minGenotypeQuality` flag, so there's no reason not to combine this with the filtering step
      - ...aaaaand the missingness-weirdness of vcftools was being compounded by my drawing the genotypes from the main vcf file rather than the filtered vcf file... AT LAST our genotypes match up to our statistics!!
  - moving forward:
    - dial down genotyping-rate filter on plink-assoc call -> try 100%, 90%, 80%
    - play with the GQ filter? -> reduce to phred 15 (95%)?


Week of 27th Feb – 3rd March

  - re-writing `vcf_to_bed.qsub` and `vcf_filter.qsub` to allow for multiple levels of filtering
  - new R notebook `MakeGenotypeTable.Rmd` to try and automate the process of scraping variant genotypes by P-value and making easy-to-read tables to help us make decisions about filtering params etc...


Week of 6-10th  March

  - rewrote `plink_fisher.qsub` to allow for multiple levels of filtering
  - bug fixes for  `vcf_to_bed.qsub` and `vcf_filter.qsub`
  - `MakeGenotypeTable.Rmd` replaced with `MakeTopSNPsWidget.qsub` & `MakeTopSNPsWidget.R` to use the power/speed of the HPC for building `DT::DataTable` widgets...
	- leading to lots of problems... knitr and pandoc are set up weirdly on the HPC and I'm going to change my mind about the extent to which extra speed is needed...
	- the job is now going to get done locally, by `MakeTopSNPsWidget.R` & `MakeTopIndelsWidget.R`


Week of 13-17th March

  - some issues are apparent with the output tables:
    - APA_6679 - NOT mixed, but marked as mixed. Fixed in `ultimate_phenotype.txt` and `ultimate_phenotype.xlsx`
    - p-vals fluctuate among hypothesis (!?)
      - missingness mis-classification?
      - ID of top candidates vary among gq15 & gq0 in a non-intuitive way
      - is there a plink bug with missingness?
  - genotypes visualised wrongly? COB_4004 & BAM_6496 – vcf files match md5sum from server to HPC: either
    - genotypes scraped wrong from vcf?
    - genotypes parsed wrong by R?
    - GT's being scraped from once-filtered VCF, not twice-filtered VCF!
  - 1 vcf -> snps & indels separated -> 2 filtering thresholds -> 3 genotyping thresholds -> 6 hypothesis tests = 72 output files in 1 'set'...
    - dropping 50% genotyped analyses == 48 files now a 'set'
  - demonstrate the difference between GQ15 & GQ20? ...venn sets? p-val overlap by threshold?
  - Make manhattan plots for scaf 23 – panel by dom. geno. rec. trend, & fisher tests at GQ20, missing 75...
  - re-run assoc with mixed ind.s flipped to cases at miss 75 or 100 and GQ20
    - `plink --bfile all_fish_version_5-1_HPC.filtered.snps.GQ20.mix --allow-no-sex --geno 0.00 --allow-extra-chr --assoc fisher --ci 0.95 --pfilter 1 --test-missing --out all_fish_version_5-1_HPC.filtered.snps.GQ20.mix_geno100 --a2-allele all_fish_version_5-1_HPC.filtered.snps.GQ20.vcf 4 3 '#' `


Week of 20-25th March

  - made new manhattan plots...
  - does GQ vary among individuals/pop.s?
  - check out:
      - SNP Scaffold23:1042779
      - Deletion Scaffold23:1042759
  -   


----

1. which test to use?
2. why not haz variant – coverage too low?
                      - some weird repetitive code?

run APA/BAM/BAVA only, with filtering & mix – Bed/Bim files are cooking

look up how to set up a power analysis for GWAS

GATK depth of coverage analysis?


n. cases – 35 (41 mix)
n. ctrls – 28 (22 mix)
prevalence – 0.6


can I collapse the SNP & indel at Scaf 23?

1. can I find BAM contamination elsewhere?
    - went to `/mnt/scratch/pitchers/eFISH/Analysis` and ran ``for i in `ls *dedup.bam` ; do echo $i >> bamgrep ; samtools view -h ${i} | grep '^@PG' >> bamgrep ; done``
    - went to `/home/willpitchers/individual_bams` on Shockly and ran ``for i in `ls *bam` ; do samtools view -h $i | grep '^@PG' >> ${i}_headrows ; done``
  - change of plan: both on HPC and Shockly, bams tested with ```for i in `ls *all_libraries.bam` ; do echo ${i} ; samtools view -H ${i} | grep --color="never" -Eo "SM:[0-9A-Z_]{3,9}" | uniq ; done``` – misplaced samples are easily spotted in this output...
  - contamination seems limited to:
    - 2 instances of APA_6676 found in `APA_6681_all_libraries.bam`
    - 2 instances of APA_6682 mislabelled as 6682 in `APA_6682_all_libraries.bam`

2. suspicious that readgroup assignment might be source of problem...
    - reminder: our read names look like [this](https://support.illumina.com/help/SequencingAnalysisWorkflow/Content/Vault/Informatics/Sequencing_Analysis/CASAVA/swSEQ_mCA_FASTQFiles.htm): `@HWI-D00731:65:C6UU2ANXX:1:1101:2046:1852 1:N:0:GAGATTCCTATAGCC`, which means;
    - @ instrument_name : runID : flowcellID : lane# : tile# : tile_x-coord : tile_y-coord <space> pair_member# : filtered? : control_bits? : index_sequence

  -


Week of 27th-31st March

  - rerunning pipeline to regenerate all the things...
  - taking this opportunity to update the pipeline to work with latest version of GATK (3.7.0):
    - allows for dropping the indel realignment step
    - speeds up filtering
    -
  -


Week 3rd-7th April

  - **data storage mishap!**
    - the current analysis directory (`/mnt/ls15/scratch/groups/efish/WILL/Pipeline6`) was accidentally `rm -rf`-ed incompletely... Since non-missing files *might* be truncated/decapitated (command was ctrl-c-ed), for safety I elected to nuke the folder, reimport the raw `..fastq.gz` files and restart the pipeline.
    - using this as an opportunity to streamline/check for software updates etc.
  -


Week of 10-15th April

  - some weirdness seems to have been occurring...
  - attach JG's permutation command to the end of pipeline


  - nperms?
  - n. PCs to incorporate?
  - repeat/understand EMP2 calculation
  - try logistic+ model flags?


plan for methods/results section?
  -

  - http://meeting.spsp.org/2016/sites/default/files/Lane,%20Hennes,%20West%20SPSP%20Power%20Workshop%202016.pdf

read about taqman assay?


Week of 17th-21st April

  -



Scaffold4236      Var-Scaffold4236-574        574    T        ADD       39        inf        11.58    4.974e-31
Scaffold4236      Var-Scaffold4236-574        574    T     DOMDEV       39          0       -11.58    4.982e-31
Scaffold4236      Var-Scaffold4236-574        574    T   GENO_2DF       39         NA        134.3    6.993e-30
Scaffold322    Var-Scaffold322-616692     616692    A        ADD       32        inf        5.777    7.593e-09
Scaffold322    Var-Scaffold322-616692     616692    A     DOMDEV       32          0       -5.777    7.608e-09
Scaffold79    Var-Scaffold79-1141895    1141895    A     DOMDEV       44      82.43        3.263     0.001102
Scaffold375    Var-Scaffold375-123643     123643    G        ADD       45      21.62        3.184     0.001452
Scaffold84      Var-Scaffold84-66482      66482    G     DOMDEV       38      141.5        3.174     0.001504
Scaffold2     Var-Scaffold2-5834893    5834893    G     DOMDEV       48  1.485e+04        3.169     0.001531





odds ratio for biphasic het/biphasic homo. alt.
odds ratio for biphasic homo. ref./biphasic homo. alt.

my plan is to list all <0.05 (permuted-p-value) markers, then bin them into these categories and simulate to encompass the range of OR values present in our output...
  - subset p<0.05 tests: `cat all_fish_version_6-1.filtered.snps.GQ20_pruned_GENO_assoc_mperm.assoc.logistic.mperm | awk '{if ($3 < 0.5) print $0}' > all_fish_version_6-1.filtered.snps.GQ20_pruned_GENO_assoc_mperm.assoc.logistic.sigP`
  - pull out SNP indices: `cat all_fish_version_6-1.filtered.snps.GQ20_pruned_GENO_assoc_mperm.assoc.logistic.sigP | awk '{print $1 ":" $2}' | sed s/Var-Scaffold[0-9]\\+-//g > sigP.intervals`
  - pull out genotypes: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -V all_fish_version_6-1.filtered.snps.GQ20.vcf -o sigP.vcf -R ${ref} -L sigP.intervals`
  - convert genotypes to table: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T VariantsToTable -V sigP.vcf -F CHROM -F POS -GF GT -R ${ref} -o sigP.genotypes.table` (then reformat with `cat sigP.genotypes.table | tr '/' '_' > temp && mv temp sigP.genotypes.table`)
  - now I need to pull out variants with the right combinations... the 1st biphasic fish in order is APA_6681, but since this is a mixed fish I'll also use the first non-mixed biphasic fish, BAM_6496... I can combine the SNP table with the assoc output to get the ID of the ref allele...
    - NB: I may have misunderstood this... I think that I need to calculate the odds ratio of the combinations analysis-wide... OK. I can read the `sigP.genotypes.table`

  - plink `--simulate` can take a missingness argument: `plink --bfile all_fish_version_6-1.filtered.snps.GQ20_pruned --missing --allow-extra-chr` to get missingness report for our data. Overall genotyping rate is 0.48022

Week of 1st-5th May

  - mean p-value `awk '{ total += $9 } END { print total/NR }' simulation.assoc.logistic`, and proportion of p-vals < 0.05 `awk '$9 < 0.05 { count ++ } END { print count/NR }' simulation.assoc.logistic`
  - some tests:
      - 1000 sims, mean p-value 0.444829, prop. <0.05 = 0.0767441
      - with none missingness, 1000 sims, mean-P =0.399, prop. <0.05 = 0.12538
      - with double no. of cases/controls; 1000 sims, mean p 0.319102, prop. <0.05 = 0.220276
  - so... what sets of simulations are needed?
      - vary sample size
          - straightforward – can just alter call at the command line
      - vary OR for class membership... calculate using only extreme pop.s vs. only mixed pop.s?
          - IVI, BAVA, MOV are the consistent populations... BAVA is the only 100% P0-absent pop.
          -


          threshold! where and why!?
          LD pattern
          graph power by n. , OR?, missingness?,


Week of 8-12th May

  -


Week of 15-19th May

  - Scaffold0:309 – tri C_T 20 & T_T 8 vs. bi C_T 25  

  - Scaffold 23 is weird...
    1. pull out all the reads mapped to Scaffold23 – ```for i in ${bamfiles[@]} ; do samtools view ${i} "Scaffold23" > `basename ${i} .bam`.Scaf23.sam ; done```
      - change of plan; sticking closer to the problem region... ```for i in ${bamfiles[@]} ; do samtools view ${i} "Scaffold23:1042000-1043500" > `basename ${i} .bam`.Scaf23.sam ; done```
    2.

  - COB-free analysis: start with:
    `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants \
        -V all_fish_version_6-1.filtered.all_variants.GQ20.vcf \
        -o all_fish_version_6-1.filtered.all_variants_noCOB.GQ20.vcf \
        -R /mnt/ls15/scratch/groups/efish/P_kings_genome/supercontigs.fasta -xl_se 'COB_[0-9]{4}' `
  - things


Week of 22nd-26th May

  - setting up to run 100,000 permutation on Shockly
	- older version of PLINK is installed - updating...
    - brew recipe for plink/1.9 not working out of the box...
    - troubleshooting not successful after ~1hr => just `wget` precompiled version from website
		- this seems to run just fine.
  -


Week of 29th May – 2nd June

  - at this point, we have tried so many different permutations (ha-ha!) of the association analysis that I'm getting confused... I'm going to need to spend some time tidying up
    - tidying/rationalization is timely – after talking with JG we agree that the best candidate list that we can obtain at this point is going to come from organized comparison of candidate intersects between analysis...
  -

Week of 3th – 9th June

  - Added sections on comparison between analyses to `GWAS_methods_summary.Rmd`:
    - intersection plots to show how few SNPs are 'candidates' from multiple analyses
    - EMP2 distributions to show the power drop of subset analyses
    - we are *really* at full stretch to support any candidates when the results are this unstable...
  - Meeting with JG -- plan is to switch focus to pattern of pop. gen. stats for the nonce, plan for further sequencing
    - to that end, updated my `calc_Fst_array` script
    - ...and the `calc_TajD` script
    - new notebook -- `Fst_Analysis_P6.Rmd` -- added to keep pop gen type analyses together

    - need to calculate pairwise LD too
    - can I set up LD & Fst calculations by p0p vs. p0a?
    - visualize the set-based logic


Week of 12th - 16th June

  - reanalysed Fst with fish grouped by *phenotype* and ignore pop.
  - searching for correlation between Fst and variant p-values
    - correlation is there, but not super-strong (~0.28)
  - re-examining the distributions of missingness, depth and quality in our data...
  - this leads us back to alignment rates. grabbing a summary of how these look with: `bamtools stats -insert -in THISFISH_all_libraries.bam`
  - *much* missingness concentrated in IVI pop... going to re-run assoc analysis without IVI fish
    - to that end, new script `vcf_to_bed_noIVI.qsub`
  - JG found realtionship between coverage & missingness...
    - are we filtering too agressively?
    - the GQ-genotype filter seems to cut deeper in IVI than other pops...
    -
  - ...
  - So, with a new canonical VCF **all** the downstream analyses will need to be re-run... also, indels need to be re-incorporated...
      - old VCFs zipped & cold-stored in `research/2015_genomic_data/P_kings_VCFs/`
      -


Week 19th - 23rd June

  - tidying up old files in dropbox and HPC research space...
    - make sure all current scripts are git-logged
    - zip up and stash obsolete scripts/outputs...
  - trying again to get `Dancer` to detect [structural variants...](http://gmt.genome.wustl.edu/packages/dancer/documentation.html)
    - looks like loooong run times for the 2nd step
    - attempted to install on Shockly... no success as yet (sadface)
      - dep.s on samtools and 'boost'... both of these installed without problem
      - seems to be some issue with the compiler. Tried 2 alternate compilers without success...
    - `dancer_config.qsub` and `dancer_max.qsub` scripts added
      -
      - running with incrementally increasing runtimes to estimate hrs needed


Week 26 - 30th June

  - dancer workflow still has problems;
    - I'm taking a VCF-chunk-type approach, slicing the single-fish BAM files into smaller bamfiles, then combining these to get ~1200 all-fish slices
    - this process takes a good while, as for each slice I need to loop over 63 fish, then loop again to index these 1-fish-bam-slices, then they can be merged into an all-fish-bam-slice, then *that* needs to be indexed, then that can be fed into `dancer-config`
-  Following JG's suggestion of scraping the literature for fish gwas to compare to


Week 3rd - 7th July

  - dancer finishes run after ~187 hrs (!) – thankful for the use of Shockly!
    - suggests that the genome contains: 2021 inversions, 20879 deletions, 160964 CTX?, 917415 ITX? and 0 insertions(?)
        - CTX and ITX are intrachromosomal and interchromosomal translocations respectively (?)
        - output format detailed [here](/Users/willpitchers/Dropbox\ \(MSU\ Efish\ Lab\/WILL/genome_assembly/garrulous-turtle/dancer_output_format.txt)
    -

    empirically what OR do we care about - more depth vs. more fish?
    make lit. results table for JG
    LAB GURU!! talk to Colin about auto-slurping JSON format notes? (check how much data LG can handle in a lump)
    JG targeted 8x coverage but we maybe need more like 12-20x ?
    genome quality score by individual/position
    what can we say about LD? decay functions for LD differ among pop.s?


Week 10–14th July

  - ?


Week 17th–21st July

  - ?

`vcftools --vcf all_fish_version_6-1.stringent.filtered.missingrm.snps.recode.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 10000 --keep apa_bam.txt --remove-filtered-all --max-missing 0 --out apa_bam_` takes ~1min

`...--ld-window-bp 5000...` takes ~35secs
`...--ld-window-bp 100000...` takes ~5.5mins
`...--ld-window-bp 50000...` takes ~3mins


Week 24–28th July

  - QQ-plots!! XX
  - model varying sample size and case:ctrl ratio XX
  - GIF lambda thing...?


Week 31st July - 4th Aug.

  - writing up analyses on labguru



Week 7–11th Aug.

  - PHENOTYPES!


Week 14-18th Aug.

  - mostly just recovering from surgery...


Week 21-25th Aug.

  - Selection of re-re-sequencing candidates
  - more labguru-ing
  -


vcftools --vcf all_fish_version_6-1.filtered.all_variants.GQ20.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 5000 --keep apa.txt --remove-filtered-all --max-missing 0 --out apa

After filtering, kept 86,716 out of a possible 9,335,131 Sites

vcftools --vcf all_fish_version_6-1.filtered.all_variants.GQ20.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 5000 --keep apa.txt --remove-filtered-all --max-missing 0 --out apa --remove-indels

After filtering, kept 73,205 out of a possible 9335131 Sites
Run Time = 214.00 seconds

vcftools --vcf all_fish_version_6-1.filtered.all_variants.GQ20.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 5000 --keep apa.txt --remove-filtered-all --max-missing 0 --out apa --remove-indels --maf 0.05

After filtering, kept 36,514 out of a possible 9,335,131 Sites

Run Time = 111.00 seconds

vcftools --vcf all_fish_version_6-1.filtered.all_variants.GQ20.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 5000 --keep apa.txt --remove-filtered-all --max-missing 0 --out apa --remove-indels --maf 0.05 --minDP 3

After filtering, kept 36,352 out of a possible 9,335,131 Sites

Run Time = 117.00 seconds

vcftools --vcf all_fish_version_6-1.filtered.all_variants.GQ20.vcf --chr "Scaffold0" --geno-r2 --ld-window-bp 5000 --keep apa.txt --maf 0.05 --remove-indels --minDP 3 --max-missing 0 --out apa

After filtering, kept 37,597 out of a possible 9,335,131 Sites

---

http://doc.goldenhelix.com/SVS/latest/svsmanual/genotype_association_tests.html



look to shared work dir to try to categorise EODs – definite absent/present... what to do about (possibly two categories of) intermediates
starts with 6 == 2009 (newest freshest)
3 == older


Week of 28th Aug - 1st Sept.

  - JG & I put our heads together and decided on our strategy for re-re-sequencing, we want to:
      1. break up the pop/EOD confounding as far as possible while
      2. maintain approximate balance in phenotype for max power...
  - problems! JG found that not all of the samples in the freezer have aged sufficiently well to be reused...
    - we'll have to make a slight adjustment to the picking-list
  - full-blown pipeline hands-off run-through a qualified success...
    - problem only appear at the `HaplotypeCaller` stage, but this matched up with some weird HPCC behaviour so locating bug proving tricksy...


Week of 4-8th Sept.

  - Labor Day.
  - some of the slices where `HaplotypeCaller` jobs failed seem to have overun, which is a little odd, but as we are about to push a lot more data through the pipeline I opted to decrease the size of the 'slices' that get fed in...
  - Started variant-calling anew... 2 runs of the `07_vcf_disco_chunk_6-1_array.qsub` script seems to have produced a complete set of output slices this time...
  - merging the slices took just under an hour – this should mean that it has time to finish ~2x the data in <4hrs.
  - a fresh-start run of `07_vcf_disco_chunk_6-1_array` resulted in only 5 'problem slices'... *however* there was an unannounced node reboot during the run...
  - slice 1295 is weird!


Week of 11th-15th Sept.

  - `HaplotypeCaller` encountered a weird error message (new one to me!!) : `/opt/software/Java/jre1.8.0_31/bin/java: Permission denied` ...900 slices failed this way!
    - running a few of these jobs interactively encounters no problems ***sigh***
    - re-re-running... seems to be working... nope. Still 36 jobs dropped!
      - *all* failed jerbs have same error – all were running on `intel16` – submitted a ticket to iCER
        - iCER recommend using a different java version... testing now...
        - the permissions problem still occurs in ~20% of jobs – waiting on iCER response.
      - meanwhile, working on the assumption that iCER may not come through for me, working on a revised array script that will check for broken jobs and re-submit them selectively... just in case I have to resort to brute-force approach...
      - Progress! X, C & P at iCER seem to have fixed the permissions issue! ...sadly not 100% error-free though:
        - 135 slices fail with 'insufficient memory' - 1 failed with 'malformed bam' - 1 failed with a 'GATK runtime error'
        - running GATK interactively seems to avoid these errors
      - brute-force might be the way to go anyway...?



https://github.com/JustinChu/JupiterPlot

https://www.biomedcentral.com/collections/evolutionarygenomics




Week of 18th-22nd Sept.

  - Mostly working on how-to GWAS MS, but...
  - emails back & forth with Chun-Min, and Patrick, *and* Xiaoge
    - tried `Java/jdk1.8.0`
    - tried `Java/1.7.0_51`
    - tried `Java/1.8.0_31_trail`


Week of 25-29th Sept.

  - early indications are the variant calling pipeline ran *without issue* <engage happy mode>
  - I SPOKE TOO SOON <sadface>
    - 100's of jobs failed to start because java failed to load **using the version they said was fixed**
      - OFFS `Java/1.8.0_31_trail` is now no longer present in `module spider`!?!
    - switching from `Java/1.8.0_31_trail` to `Java/1.8.0_31` and re-re-re-running <sigh>
  - vcf file lengths are inconsistent:
      - 9339834	/mnt/research/efish/2015_genomic_data/P_kings_VCFs/all_fish_version_6-1.filtered.all_variants.vcf
      - 10137913	JulyAssoc/all_fish_version_6-1.vcf
      - 10225321	Assoc/all_fish_version_6-1.vcf
      - 10229992	SeptAssoc/all_fish_version_6-1.vcf
      - 10238363	AugAssoc/all_fish_version_6-1.vcf
      - 10239061	/mnt/research/efish/2015_genomic_data/P_kings_VCFs/all_fish_version_6-1.vcf
      - 10239191	LateSeptAssoc/all_fish_version_6-1.vcf
      - 10239271	LaterSeptAssoc/all_fish_version_6-1.vcf
  - looking into the lengths of vcf *slices*...
    - 70, 484, 1716 are only the header, 2280 overran walltime
    - 2201, 2238, 2283, 1003, 2116, 1633, 2355  think they completed successfully...?
    - these numbers aren't making sense to me...
      - `bash crapsearch`
      - `cat crap | cut -d '-' -f 3 > crapnums`
      - `uniq crapnums > unicrap`
      - `wc -l *crap*`


Week of 2nd-6th Oct.

    - I have tried shortening the length of the slices to avoid the (rare) over-runs, and adding an `if` statement to check that the output slice continues past the header... deletes output if not.
    - This ^ wasn't sufficient – adding another check to the `slice_check.qsub` script to search the stdout file for the string "Done." which (with the capitalisation and the fullstop) is indicative of GATK having finished the job.
      - 10137913 JulyAssoc/all_fish_version_6-1.vcf
      - 10225321 OctAssoc/all_fish_version_6-1.vcf
      - 10229992 SeptAssoc/all_fish_version_6-1.vcf
      - 10239191 LateSeptAssoc/all_fish_version_6-1.vcf
      - 10239271 LaterSeptAssoc/all_fish_version_6-1.vcf
      - 10238048 Oct2Assoc/all_fish_version_6-1.vcf
      - 10238172 Assoc/all_fish_version_6-1.vcf
      - 10238363 AugAssoc/all_fish_version_6-1.vcf
    - there is no pattern that I can detect in which nodes are running when jobs mystery-fail...
      - 10238239 04/10/2017 version !?!?!
      - How many rows *aren't* header?
        - 04/10/2017 10233544, July  10133218, August 10233668, late Sept. 10234496, early Sept.  10234576, October #1 10220626, October #2 10233353, October #3  10233477
        - All headers are 4695 lines long, so we're definitely losing lines from the end...
  - **NEW DATA** has arrived (1st lump of 136 `..fastq.gz` – 2nd run data expected on 10/10/2017)
    - `scp`-ed to the research-space-scratch – all `md5sum` checks passed
      - quality checks running... added a little Rscript to use the `fastqcr::qc_aggregate` function to create a summary output
        -
      - trimmomatic array running...
        -
      - alignment array running...


Week of 9th-13th Oct.

  - **NEW DATA**
    - deduplication array running...
      -
    -
  - ** *NEW* NEW DATA** has dropped - 72 new `..fastq.gz` archives
    - `sftp`-ing from `titan`... all `md5sum` checks passed!
    - trimming – run without error msg
    - aligning... a few over-runs
    BAM_6558_Extract_S17_L007_R2_se.aligned.dedup.bai missing pe files
    - trying just to let pipeline run on new data in hopes of getting a semi-cromulent new VCF for JG so that hew can add to his talk for 14/10/17
      - output ("big") VCF ends up at 506M in size... <sigh> Obviously more than just the `BAM_6558_Extract_S17...` data didn't make it through... So much for fast and dirty.


Week of 16th-22nd Oct.

  - OK, follow the workflow through:
    - `..fastq.gz` - 432 archives
    - `..fastq` - 432 files
    - `..trimmed.fq` - 864 files (432 SE & 432 PE)
    - `..aligned.sam` - 648 files (432 SE & 216 PE; `..R1_pe..` & `..R2_pe..` are co-aligned)
    - `..dedup.bam` - 646 files =( `BAM_6559_S60_L007_R1_pe.aligned.sam` & `BAM_6501_ATTACTCG-CAGGACGT_L001_R1_se.aligned.sam` failed. Looks like an error in the alignment step caused `Picard` to fail out...
      - seems to be a node glitching, as both run fine interactively...?
      - re-running the entire array seem to work fine (?)
    - `..dedup.bam` - 648 files now
    - `..recalibrated.bam` – 636 files now. Should be 648... problem seems to be resource over-run.
      - `..recal_data.table` - 637 then `..post_recal_data.table` - 636 , so the problem happens *early*
        - problems at 355, 379, 404, 544, 571, 574, 583, 598, 604, 607, 610, 616, 622, 647... going for a full-array re-run
    - `..recalibrated.bam` - 647 files now... #headdesk
    - GODSDAMMIT THE HPC WILL NOT STOP DROPPING NODES!!?!
  - OK. <deep breath> Back to the start because I'm confused...
    - `01_trimmomatic...` 432 archives in (216 pairs)
      - 864 `..trimmed.fq` out - (432 SE & 432 PE)
    - `02_alignment...` 864 fastq's in
      - 648 `..aligned.sam` out - files (432 SE & 216 PE)
      - checking formatting with `ValidateSamFile`...
      - I have also added this check to the alignment script so that errors will be recorded in the logfiles
    - `03_deduplication...` 648 sam files in
      - running full array...
      - 603 of 648 jobs report completing with no errors.
      - `BAM_6568_S22_L007_R1_se.aligned.sam` didn't write a logfile, but ran without error interactively (!)
      - 44 failed jobs resubmitted – failed again...
        - different error messages this time – suggests that there's error in the sam format
          - trying to re-run the alignment interactively... no errors before they get killed.
          - `BAM_6604_TCCGGAGA-TAATCTT_L006_R2_se..` & `BAM_6871_S27_L007_R1_pe..` as test cases:
            - sam files both pass interactive `ValidateSamFile` - "No errors found"

    - The problems seem to start with bwa-mem, but don't occur when I run interactively on intel16... so:
      - restricting the array to intel16 nodes
      - these all have 128GB available, so may as well use all of it (testing suggests that this is 10X overkill for most cases)
      - keeping to a 4hr time-slot to reduce queue-time, but multi-threading at 4X for speed
      - `02_alignment...` - 648 files *with no errors* at last!
      - `03_deduplication...` <sigh>
        - only 643 logfiles!
          - 639 report "No errors" - but 642 dedup-ed bam files exist (!?)
          - 2 report "killed" - 348,408
          - 2 with PICARD error messages - 384,607
          - 9 jobs disappeared without trace - missing array values are 219,419,584,608,619,648
          - bam files missing for: `APA_6676_GAGATTCC-ATAGAGGC_L002_R1_se..`, `APA_6737_ATTCAGAA-AGGCGAAG_L002_R1_pe..`, `BAM_6561_Extract_S65_L008_R2_se..`, `BAM_6840_S11_L007_R2_se..`, `BAM_6865_S25_L007_R2_se..`, `BAM_6867_S54_L006_R1_se..`
        - killed job re-submitted with +1hr walltime...


Week of 23rd-27th Oct.

  - continuing to smack face against wall:
    - `APA_6737_ATTCAGAA-AGGCGAAG_L002_R1_pe..`, `BAM_6865_S25_L007_R2_se..`, `BAM_6840_S11_L007_R2_se..`, `BAM_6867_S54_L006_R1_se..`, and `BAM_6561_Extract_S65_L008_R2_se..` all completed successfully on re-submission
      - these jobs failed on `ifi` (intel14) nodes, which is suggestive, although not all jobs on `ifi` failed – precautionarily limiting future jobs to `lac` (intel16) nodes
    - `APA_6676_GAGATTCC-ATAGAGGC_L002_R1_se..` job **still** didn't make an output file on re-re-run... but running interactively (on intel16) gives "No errors found"...
      - I'm not really sure where to go with this without logfiles to look at because I can't reproduce the error!?
    - I'm submitting `04_base_score_recal..` now that I have the full set of checked-out bam files
      - the reference folder has been emptied – I'm presuming by a timed scratch-cleaning script – so repopulating from the research space... also need to regenerate `..fai` & `..dict` index files...
    - `04_base_score_recal..` has written 636/648 logfiles, but all contain "no errors" messages
      - missing are 227,394,427,458,545,572,575,584,599,608,617,620,623,627,648 (why too many!? re-subbed just in case)
      - only 630 output bams exist (!?)
      - missing bams: `APA_6675_GAGATTCC-TATAGCC_L001_R1_pe..`,`APA_6675_GAGATTCC-TATAGCC_L002_R1_pe.. `, `APA_6675_GAGATTCC-TATAGCCT_L001_R1_pe..`, `BAM_6865_S25_L007_R1_se..`, `BAM_6865_S25_L007_R2_se..`, `BAM_6865_S62_L008_R2_se..`, `BAM_6867_S35_L007_R2_se..`, `BAM_6867_S54_L006_R2_se..`, `BAM_6868_S36_L007_R2_se..`, `BAM_6869_S68_L008_R1_pe..`
      -
    - testing write-outs to both temp scratch and the ffs17 filesystem...
      - both of these options work well when run interactively – the temp-scratch would be the better option if all else is equal as there are *tight* limits on space on the ffs17 filesystem...
      - subbed full-arrays... queue seems pretty full so might have a delayed start...
        - pleasingly, these arrays produce bam files with no `diff`
        - some over-running jobs in the temp-scratch array...
        - first run of the ff17 array hit problems with over-filling the quota -- going to try and fix this by explicitly `rm`-ing files right after they're copied back to the working directory
      - ffs17 *still* over-filling the quota, even limited to 50 jobs at a time..
      - TMPDIR seems to disappear between writing out the `..bam` files and `mv`ing then back to scratch for some jobs(!?)


Week of 30th Oct.-3rd Nov.

  - Moving forward with the BSQR problems...
    - JG suggested that permissions mismatch between scratch and TMPDIR might be to blame...
      - https://wiki.hpcc.msu.edu/display/hpccdocs/Permissions+on+HPCC+File+Systems :
        - Scratch "Directories are created as world-readable by default. If you do not wish this to be the case, you will need to alter the permissions of any top-level directories you create in scratch."
        - TMPDIR space "Directories are creates as world-readable by default, but the scheduler deletes the contents of $TMPDIR after a job exits."
      - running interactively I see the same permissions `-rw-rw-r--` for both scratch & TMPDIR...
    - **Breakthrough** TMPDIR has some *weird* behaviour with copying – `cp TMPDIR/* ./` doesn't work, but `cd TMPDIR && cp * there` does work.
      - This makes none sense to me, but this behaviour is repeatable interactively, and when qsub-ed, so it hopefully fixes the array...
      - as these jobs start to complete, I've tested 3 (haphazardly chosen) against interactively built versions – no `diff` (woo!)
      - 31 jobs have over-run (on first submission)...
    - another good suggestion from JG: now that I'm restricting array jobs to reun on the Laconia (intel16) nodes, it makes sense to optimise the use of memory and multi-threading such that each job uses a whole node: i.e. `nodes=1:ppn=28,walltime=04:00:00,mem=128gb`... making this the new standard should give me a small speed boost
      - OK. Somehow there are *still* some files missing...
        - e.g. `BAM_6867_S54_L006_R2_se.aligned.dedup.bam` didn't get a `..recalibrated.bam`... but all 4 commands run interactively with "no warn messages" and validates with "No errors found" !
  - Met with Pat Bills from iCER – we have decided on a 2-pronged approach:
    - I am going to make a clean run of the array with some added reporting statements added in order that we can build a more complete picture of any cascading errors...
    - meanwhile Pat Bills is going to attempt to reproduce our errors/outputs independently – I've made copies of the input files on scratch for him to use...


Week of 6-10th November

  - starting by building a job-tracking spreadsheet as discussed with JG & PB
    - took this info. to iCER office hours; met with Michelle Szidik, who was able to **actual** diagnose the problem!
    - the scheduler (PBS) was running the array from a version of the script that contained git merge annotations that the shell was trying and failing to interpret as commands...
    - base recal. is done after a rerun and selected re-re-run for walltime over-runs
      - to avoid this in future I've split the 4 GATK commands in the `04_base_score_recalibration_array` script over 2 scripts, effectively doubling the max. runtime. This needs to be tested, but I shall get back to that once I have a VCF from the new fish and can put tests in the queue without slowing down anything mission-critial
      - the remaining mystery error – seems to be a down-stream consequence of partial re-run. See: https://gatkforums.broadinstitute.org/gatk/discussion/7057/baserecalibrator-stops-before-finishing-and-does-not-throw-an-error
  - running vcf discovery...
    - some errors popping up...
      - problem #1 - `07_vcf_disco_chunk_7_array.qsub` resubmits each version of itself if it sees "warn" messages from GATK `ValidateVariants`. This is a problem where there are many (>15) potential alleles at a locus GATK warns us, and the script gets trapped in an infinite looped
      - problem #2 - `BAM_6559_all_libraries.bai` is reported as truncated. My suspicion is that multiple infinite-looping I/O has broken this file. I shall re-index manually with `samtools index -b BAM_6559_all_libraries.bam`.
      - problem #3 – a few jobs complained that the file they were trying to write could not be written because the output directory wasn't writeable... this stems from a typo in the `if` statement that re-submits failing jobs – I fix.
  - additionally, new sequencing data archives copied to the backup-backup NAS drive because I'm paranoid


Week 13-17th November

  -
  - maybe pre-make the bam slices?
    - 2min, 1.5min, 1.5min, 1.2min, 1min per slice for the first 5 slices tested...
    - wrote `make_bam_slices.qsub` to spawn an array, each job to write out 100 slices – testing n=1 as an interactive job...
      - 100 single-slice BAMs built in 168mins
    - array completed, but the `BAM_6559_all_libraries.bam` file is again a problem...
      - reindexed and resubmitted – problem repeats
      - testing the individual-library BAMs for malformation etc....
      - `BAM_6599_all_libraries` *still* the problem... regenerating it from the single-library bams...
    - re-running `make_bam_slices`
      - New problem fish: `BAM_6562` is now causing the same errors
      - this may indicate just how much simultaneous I/O scratch can take:
        - `make_bam_chunks.o49416415-51` failed 62/100 times
        - `make_bam_chunks.o49416415-42` failed 70/100 times
        - `make_bam_chunks.o49416415-43` failed 100/100 times, as did `make_bam_chunks.o49416415-44`, and 45-50
      - for now I'm going to remake `BAM_6562` and try again, but I'll need to come up with a way to avoid this in future
        - aaaaand `BAM_6562` is now the problem again... re-re-making `BAM_6562_all_libraries`
        - 42, 43, 44, 45, 46, 47, 48, 49, 50, 51
          - testing a version of `make_bam_slices` that writes to the ${TMPDIR} to see if this streamlines the `..bai` corruption issue...
      - problem solved! (hopefully) by invoking the `--disable_bam_indexing` option in GATK
        - all 5162 bam slices are written!
        - edited `08_vcf_disco_chunk_7_array.qsub` to point each array instance at a single BAM slice rather than the list of `FISH_XXX_all_libraries.bam` files
        - new VCFs generated
          - slices 1500 & 2500 failed to index, and 3799,3989,4102 & 4152 couldn't *open* the index
          - explicitly indexing bam slices right after creation from now on!
          - 4102 timed out, but that explanation doesn't fit for 3989, 3799

          - 4 jobs errored out – they claim that the files are incomplete - re-made and re-indexed BAM slices


Weeks 20-24th November & 27th Nov. - 1st December

  - merging new WG VCF file...
    - GATK reports slices 3859 & 4097 are empty...?
      - re-running the calling fixes those slices
    - GATK now complaining about 2 variants that have gotten themselves in the wrong order on scaffold 17...
      - sorting the VCF fixes this problem

  moving forward; a wishlist of items...
  - re-ran vcf_filter with GQ filter turned up to 30 – remade bed/bim/fam – reran plink_fisher

  <--- THANKSGIVING --->

  - candidate list == loci where association is perfect given number of genotypes we have at that locus...
    - biallelic script from JG -> this should *fixed* differences between phenotypes
      - `biallelic.pl` needs a group assignment file – building one with:
        - `awk -F "[\t ]+" '{ P=$6 - 1 ; if ( P=="0" )  N="absent" ; else  N="present" ; print $1, " ", P, " ", N }' all_fish_version_7.fam > fish_pheno_for_biallele.txt`
      - run with `perl /mnt/research/efish/2015_genomic_data/biallelic.pl all_fish_version_7.sorted.filtered.snps.GQ30.filt_pass.vcf.gz fish_pheno_for_biallele.txt`... problems.
        - perl script is looking for a command called `vcf-subset`... is this from VCFtools?
          - it *does* seem to want VCFtools... but is running for a good while. Might need to qsub it?
    - knock this nonsense on the head <- using PLINK instead
      - PLINK with: `plink --bfile ../all_fish_version_7.sorted.filtered.snps.GQ20 --allow-no-sex --geno 0.25 --allow-extra-chr --freq case-control --out scaffold137_geno75_jrg_test --chr Scaffold137 --a2-allele ../all_fish_version_7.sorted.filtered.snps.GQ20.vcf 4 3 '#'  --keep-fam ../wave2.fam`
      - then munge with `awk '{print $2, $5 * $7,(1-$5)*$7, $6 * $8,(1-$6)*$8}' scaffold137_geno75_jrg_test.frq.cc > counts_and_freqs.txt` and `paste scaffold137_geno75_jrg_test.frq.cc counts_and_freqs.txt >  analyze_this.txt`

  - look into rerunning pipeline with 2015 fish? (use v.6.1 VCF)
    - combined the vcfs with `java -Xmx120g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R ${ref} --variant all_fish_version_6-1.vcf --variant all_fish_version_7.sorted.filtered.snps.GQ20.vcf  -o all_fish_versions_6-1_plus_7-filtered.vcf --genotypemergeoption UNIQUIFY`
    - submitted `vcf_to_bed` X problems!! old vcf has individuals labeled differently than new one (APA_6675.variant vs. APA_6675)
      - using sed to fix this... then retrying `vcf_to_bed`... nope.
    - instead trying `java -Xmx120g -Djava.io.tmpdir=${TMPDIR} -cp /mnt/home/pitchers/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ${ref} --variant all_fish_version_7.sorted.filtered.snps.GQ20.filter_pass.vcf --variant all_fish_version_6-1.vcf -out all_fish_versions_6-1_plus_7-filtered.vcf`
      - initial failure because GATK won't combine filtered with unfiltered VCFs...
      - used `all_fish_version_7.sorted.vcf`
      - vcf_filter doesn't work on all_fish_versions_6-1_plus_7.vcf because it's not sorted...  sorting with `java -jar ~/picard.jar SortVcf I=all_fish_versions_6-1_plus_7.vcf O=all_fish_versions_6-1_plus_7.sorted.vcf`
      - trying `vcf_filter` and `vcf_to_bed` with `all_fish_versions_6-1_plus_7.sorted.vcf`...
      - also problems. GATK finds the VCF malformed because the header specifies 76 genotypes but some positions contain only 63 genotypes... do I need to merge at the BAM file stage to get this to work...?
    - OK. attempt the third: `java -Xmx120g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T CombineVariants -R ${ref} --variant all_fish_version_6-1.vcf --variant  all_fish_version_7.sorted.vcf -o all_fish_versions_6-1_plus_7.vcf --genotypemergeoption UNIQUIFY`
      - no warn messages. trying to make a set of plink files... success! Making second set from filtered output

  - It would be helpful to me to clarify the group memberships of suspect-EOD fish...

  - set threshold WRT highest *possible* pvals given nfish
    - looking into adjusting threshold based on missingness – ran `plink --missing --bfile all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass --allow-extra-chr --out all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.bed --keep-fam wave2.fam`
      - output (`all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.bed.lmiss`) can be joined to the assoc output...
      - re-used chi-square simulation from `PowerNotebook.Rmd` to build a tbale of min. pvals by missingness (I made the assumption that missingness does not alter balance, which is probably not true at every locus, but will make comparison against empirical pvals inherently conservative, so I'm OK with it)

  - run a simple APA vs. BAM GWAS
    - all this needs is an alternate `..fam` file to reclassify POP as PHENO...
    - Ugh. I have  malformed VCF again... remaking with `java -Xmx120g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V all_fish_version_7.sorted.filtered.snps.GQ20.vcf --setFilteredGtToNocall -o all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.vcf`
    -


    -

    - plot alternate hypotheses!!
    - fix black-grey deal
  - LD and Fst calcs need to be rerun


  - where are the pvals better for EOD than for pop?
    - pass pop into plink as covariate
    -


Week of 4-8th December

  - Is missingness still problematic?
    - missingness report with `plink --bfile ${input_data} --allow-no-sex --geno 0.25 --allow-extra-chr --keep-fam wave2.fam --out ${input_data} --missing`
      - freq. missingness per fish *not* predictable by phenotype (p-val=0.50)
      - freq. missingness per fish *not* predictable by population (p-val=0.84)
      - freq. missingness per scaffold gets larger with smaller scaffolds (p-val<0.001), but only slightly (coef=2.9e-5)
      - var. missingness per scaffold also gets larger with smaller scaffolds (p-val<0.001), but only slightly (coef=2.5e-6)
      - these results hold (qualitatively) true for input_data that's been filtered at GQ30 as well as GQ20.
      - missingness per fish is >0.05 & <0.3 at GQ20 but jumps to >0.1 & <0.4 at GQ30
    - is JG seeing very different results with `vcftools`-generated missing reports...?
      - he used; `vcftools --vcf ${input_data}.vcf --max-missing 0.5 --minQ 20 --recode --recode-INFO-all --remove-filtered-all --keep ../wave2.fam --out ${input_data}.wave2_GQ20_hard_filt`, and then; `vcftools --vcf ${input_data}.wave2_GQ20_hard_filt.recode.vcf --missing-site`, and then; `vcftools --vcf ${input_data}.wave2_GQ20_hard_filt.recode.vcf --missing-indv`

      - `vcftools --vcf ${input_data}.vcf --remove-filtered-geno-all --remove-filtered-all --keep wave2.fam --max-missing 0.5 --minQ 20 --recode --recode-INFO-all`

      - plink vs. VCFtools
        - first attempt finds that they *do not* agree...
        - details in `Missingness2017.Rmd`
        -


      estimate pop missingness by alleles
      where is the sweet spot for filtering?
      big pops – plenty of scope for private alleles to exist *per fish*

      explicit list of filters on way in to PLINK:

    - `GQx` -- GATK option in `vcf_to_bed`
      -
    - filterExpression -- GATK options in `vcf_filter`
      - QD < 2.0  --  'QualByDepth' normalizes call confidence by depth of sample reads supporting a variant
      - FS > 60.0  --  'FisherStrand' is the strand bias estimated using Fisher's Exact Test; the output is a Phred-scaled p-value.
      - MQ < 40.0  --  'MappingQuality' is a threshold for the mapping qualities of the reads supporting the variant
      - MQRankSum < -12.5  --  'MappingQualityRankSumTest' is a test for for mapping qualities of REF versus ALT reads. The ideal result is a value close to zero, which indicates there is little to no difference. A negative value indicates that the reads supporting the alternate allele have lower mapping quality scores than those supporting the reference allele, and vice versa.
      - ReadPosRankSum < -8.0  --  'ReadPosRankSumTest' is a test for relative positioning of REF versus ALT alleles within reads, since sequencers make the most errors at the ends of reads.  A negative value indicates that the alternate allele is found at the ends of reads more often than the reference allele and vice versa.
    - geno 0.05  --


Week of 11-15th December

  - rebuilding colour-coded genotype tables of 'top snps'
    - weirdness about allelic, dom, trend, geno topsnps lists... lists seem incomplete 5000ish when coded for 10000
    - best cases Fisher:
      - 121	Scaffold704	64680
      - 205	Scaffold31	2433337
      - 216	Scaffold229	387059
      - 772	Scaffold40	1654521
      - 917	Scaffold162	1029026
  - running associations with no missingness requirement:
    - `plink --bfile ${input_data} --allow-no-sex --geno 0 --allow-extra-chr --assoc fisher --pfilter 1 --test-missing --out ${input_data}_geno00 --a2-allele ${input_data}.vcf 4 3 '#'  --keep-fam wave2.fam`
  - running associations in logistic mode with pop. as covariate:
    - `plink --bfile ${input_data} --allow-extra-chr --geno 0.10 --allow-no-sex --logistic mperm=10000 --mperm-save --out ${input_data}_geno90_mperm --pfilter 1 --a2-allele ${input_data}.vcf 4 3 '#'  --keep-fam wave2.fam`
    - `plink --bfile ${input_data} --allow-extra-chr --geno 0.10 --allow-no-sex --logistic mperm=10000 --mperm-save --out ${input_data}_geno90_mperm --pfilter 1 --a2-allele ${input_data}.vcf 4 3 '#'  --keep-fam wave2.fam --covar ../population_metadata.fam --covar-number 4`
  - what do the associations look like *within* Bambomo?
    - make a Bambomo-only plink-file set: `java -Xmx30g -cp $GATK -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R ${ref} -V ${input_data}.vcf --sample_expressions BAM -o ${input_data}_BAMonly.vcf`, then passing the vcf into `vcf_to_bed.qsub`

  - remake both waves plink files...
  - outputs of:
    - assoc 90% genotyping filtered GQ20 wave2 logistic
      - `all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_geno90_logi.assoc.logistic`
    - assoc 90% genotyping filtered GQ20 wave2 logistic with pop covar
      - `all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_geno90_popCV.assoc.logistic`
    - assoc 90% genotyping filtered GQ20 BAM logistic
      - `all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_BAMonly_geno90.assoc.logistic`


18th-31st December

  # The File Lose-ening of 2017

  ## Non-issues

    - 13 files gone from the `FOR_PAT` directory - I'm just gonna nuke the folder
    - 697 VCF/sam 'slices' from the `Pipeline6` directory... these are likely defunkt anyway
    - 60 files from `Pipeline7/vcf_chunks_attempt1` - a testing folder that I probably ought to clear out anyway
    - 7 FastQC files that are trivially easy to regenerate
    - 46 process files spun off by GATK (`..table`, `..metrics.txt`, `..stats` & `..recal_plots.pdf`)
    - 5 `..log` plink output files

  ## Minor issues

    - 211 version 7 bam 'slices' these will need to be regenerated if I want to re-do variant calling...
    - 256 vcf chunks...
    - 36 `..aligned.sam` files from current pipeline
    - 17 `..trimmed.fq` files - trivially east to regenerate
    - 630 various index files

  ## Issues (sadface)

    - files that need to be re-made:
      - all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_BAMonly_geno90.missing
      - all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_geno90_PopCV.assoc.fisher
      - all_fish_version_7.filtered.snps.GQ20.vcf
      - all_fish_versions_6-1_plus_7_geno00.REC.model.reformatted.csv
      - all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass.filtered.snps.GQ20_BAMonly_geno00.assoc.logistic.topsnps.csv
    - files that *may* need to be re-made:
      - pop_gwas/all_fish_version_7.sorted.filtered.snps.GQ30.filt_pass_pop_geno100.TREND.model.reformatted.csv
      - pop_gwas/all_fish_version_7.sorted.filtered.snps.GQ20.filt_pass_pop_geno100.REC.model
    - samples.list2 – I shall have to check why I made this variant on the `samples` file

  ## Raw sequence files lost from scratch

    - BAM_6603_TCCGGAGA-AGGCGAAG_L001_R2_001.fastq.gz
    - APA_6684_ATTCAGAA-CCTATCC_L004_R1_001.fastq.gz
    - BAM_6598_TCCGGAGA-ATAGAGGC_L002_R1_001.fastq.gz
    - BAM_6604_TCCGGAGA-TAATCTT_L005_R2_001.fastq.gz
    - BAM_6499_ATTACTCG-AGGCGAA_L005_R2_001.fastq.gz
    - BAM_6603_TCCGGAGA-AGGCGAA_L005_R2_001.fastq.gz
    - 20170929_DNASeq_PE/BAM_6564_S17_L002_R2_001.fastq.gz
      - this is the most critical; I've sftp-ed it from the backup drive

  - Phew. After that...
    - ` samtools view -b input.bam "Scaffold0:1-7257134" > output.bam ` to pull scaf0 reads out of bam files

---

## 2018

1st-12th January

- GATK forums on parallelism:
  - found [this article](https://gatkforums.broadinstitute.org/gatk/discussion/4026/parallel-running-in-gatk) which ends "Hah, this thread was old! We've learned a lot since then..." and links to:
  - [this article](https://gatkforums.broadinstitute.org/gatk/discussion/9659/what-are-the-smallest-units-i-can-break-whole-human-genomes-into-for-scatter-gather#latest) which says "If you can wait ... to use GATK4 it will be a massively better and less painful solution [than scatter-gather]"...
    - scatter-gather with smaller-than-chromosome "isn't wrong in principle", but Geraldine warns that "It's not so much a question of size of region as it is a question of what you risk interrupting -- namely, any multi-nucleotide variant the spans a boundary between segments."
    - this problem can be avoided by using intervals only "...bounded by stretches of Ns in the genome as well as by regions where we know we can never make confident calls." but of course they have this kind of high-resolution knowledge for just the human genome ATM.
  - I posted [here](https://gatkforums.broadinstitute.org/gatk/discussion/9659/what-are-the-smallest-units-i-can-break-whole-human-genomes-into-for-scatter-gather) and will update these notes when one of the Broad team responds.

  - replaced missing `..fastq.gz` files from `efishbackup`
  - Re-running md5sum checks in case of corruption in the files that I *didn't* replace...
    - `APA_6684_ATTCAGAA-CCTATCC_L004_R2_001.fastq.gz` & `APA_6684_ATTCAGAA-CCTATCC_L004_R2_001.fastq.gz` failed and were replaced
      - ...then trimmed, aligned, recalibrated
      - alignment issues seem to be happening with some of the single-end read files... The ancestor gzips are still good; the gunzipping seems to have glitched (shrug) -> re-extraction fixes errors
    -

  - aiming for a number of comparable analyses, all restricted to scaf0 for speed:
    - Pipeline8 standard
    - Pipeline8 gVCF
    - Pipeline8 standard, aligned to BioNano-super-scaffolded ref.
    - Pipeline8 gVCF, aligned to BioNano-super-scaffolded ref.


  - aligning with super-scaffolded reference... using `/mnt/research/efish/lab_data/assemblies/genomes/Para_king_2015_013/super_scaffold`
    - copied all files to `/mnt/ls15/scratch/groups/efish/P_kings_SuperScafGenome`, and ran `bwa index` on `Para_king_2015_013_20_40_15_90_3_superscaffold.fasta`
    - added `alt_ref_alignment_array.qsub` script – modified only to point at alternate ref. and send output to a new folder


samples=`for i in ${GVCFs[@]} ; do echo --variant $i ; done | xargs`

`java -Xmx60g -cp $GATK -jar /opt/software/GATK/3.7.0/GenomeAnalysisTK.jar -T CombineGVCFs -R ${ref}  ${samples}  -L Scaffold0:1-7257134 > wave2_all_fish_scaf0.g.vcf`


`time java -Xmx60g -cp $GATK -jar /opt/software/GATK/3.7.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R ${ref}  ${samples}  -L Scaffold0:1-7257134 > taco.vcf`
