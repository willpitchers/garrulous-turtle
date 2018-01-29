# Analysis Pipeline for Variant Discovery and GWAS of *Paramormyrops kingsleyae* EOD Structure

This file is intended as a map of the pipeline; it should grow with the pipeline as I build it to act as aide-memoire for me and JG.

## Overview

The rough outline goes like this:

**QC** --> **Alignment** --> **Sorting/Indexing** --> **Base Score Recalibration** --> **Sample Merging** --> **Variant Calling** --> **Statistical Analysis**

One added level of complexity is to be found in the **Variant Calling** stage given the need to bootstrap ourselves a SNP database in order for GATK to work optimally...

Once data reaches the **Statistical Analysis** stage it may become necessary for the pipeline to branch as it seems likely that I'll need to use more than one program to generate the statistics that we're aiming at...

***UPDATE: it is not yet completely clear to us which of the variations on our pipeline is most appropriate...***

## Version 1 (completely parallel)

   |  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
1. | **QC** | `trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and trims
2. | **Alignment** | `alignment_array.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Reference.fa` | `Sample_library_Lane_RX_Ye.aligned.sam` | pulls details from filenames to make an `@RG` tag and then runs alignment
3. | **Sorting/Indexing** | `deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.aligned.sam` | `Sample_library_Lane_RX_Ye.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
4. | | `indel_realign_array.qsub` | GATK/3.4.46 | | `Sample_library_Lane_RX_Ye.dedup.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.realignment_targets.list` & `Sample_library_Lane_RX_Ye.realigned.bam` | RealignerTargetCreator ID's targets, IndelRealigner realigns
5. | | `base_score_recalibration_array.qsub` | GATK/3.4.46 | '--run_without_dbsnp_potentially_ruining_quality' | `Sample_library_Lane_RX_Ye.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` & `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.recalibrated.bam` | 2 passes with BaseRecalibrator (with the infamous no-dbsnp flag), then AnalyzeCovariates prints the plots/stats, then PrintReads writes the output bam
6. | **Variant Calling** | `vcf_discovery_array.qsub` | GATK/3.4.46 | `--genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30` | `Sample_library_Lane_RX_Ye.recalibrated.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.raw_variants.vcf` | HaplotypeCaller does what it says on the tin
7. | | `compress_vcfs.qsub` | tabix/0.2.6 & vcftools/4.2 | | `Sample_library_Lane_RX_Ye.raw_variants.vcf` | `Sample_library_Lane_RX_Ye.raw_variants.vcf.gz` | each vcf file is compressed (for vcftools compatibility) and then tabix-indexed
8. | | `vcf_merge_samples_array.qsub` | tabix/0.2.6 & vcftools/4.2 | `--remove-duplicates` | `Sample_library_Lane_RX_Ye.raw_variants.vcf.gz` | `Sample_all_libraries.vcf.gz` | vcf files merged at the level of population
9. | | `merge_all_vcfs.qsub` | tabix/0.2.6 & vcftools/4.2 | | `Sample_all_libraries.vcf.gz` | `all_variants_merged_$date$.vcf.gz` | this is a long step; merging all the sample-level vcfs into one file for passing to SKAT
10. | **Statistical Analysis** | `calc_Fst_array` | vcftools/4.2 | `--fst-window-size` & `--fst-window-step` | `all_variants_merged_${date}.vcf` | `Fst_POP1_vs_POP2.windowedXkb.stepYkb.weir.fst` | Fst stat. calculated for all pairwise between-pop comparisons. Options passed into output filenames.
11. | | `plink_vcf_convert.qsub`, `plink_prep.qsub`, & `plink_fisher.qsub` | vcftools/4.2 & plink/1.07 | `--maf` & `--geno` (missingness) | `all_variants_merged_${date}.vcf` | `all_variants_merged_${date}.assoc.fisher`, `*.ped`, `*.bed`, `*.bim`, `*.fam`, `*.log`, `*.map` & `*.nosex` | Association with Fisher's exact test and simplified presence/absence phenotype data.
12. | **Find Structural Variants** | `bam_merge_samples_array.qsub ` | picardTools/1.89 | | `*.bam` | `${individual}_all_libraries.bam` | ...
13. | | `breakdancer.qsub` | BreakDancer/1.1.2 & SAMTools/1.2 | `-t -q 10 -d` | `${individual}_all_libraries.bam` | `breakdancer_${date}_analysis.cfg` & `...ctx` | Detects structural variants...

---

## Version 2 (defunct)

---

## Version 3 (GVCF -> joint genotyping with GenotypeGVCFs)

   |  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
1. | **QC** | `trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and trims
2. | **Alignment** | `alignment_array.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Reference.fa` | `Sample_library_Lane_RX_Ye.aligned.sam` | pulls details from filenames to make an `@RG` tag and then runs alignment
3. | | `fix_readgroup_array.qsub` | picardTools/1.89 | `${PICARD}/AddOrReplaceReadGroups.jar` | `Sample_library_Lane_RX_Ye.aligned.sam` | `Sample_library_Lane_RX_Ye.aligned.rg.sam` | overwrites the `@RG` tag to make it *GATK-compatible*
4. | **Sorting/Indexing** | `deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.aligned.sam` | `Sample_library_Lane_RX_Ye.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
5. | | `indel_realign_array.qsub` | GATK/3.5.0 | | `Sample_library_Lane_RX_Ye.dedup.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.realignment_targets.list` & `Sample_library_Lane_RX_Ye.realigned.bam` | RealignerTargetCreator ID's targets, IndelRealigner realigns
6. | **Base Score Recalibration** | `base_score_recalibration_array.qsub` | GATK/3.5.0 | '--run_without_dbsnp_potentially_ruining_quality' | `Sample_library_Lane_RX_Ye.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` & `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.recalibrated.bam` | 2 passes with `BaseRecalibrator` (with the infamous no-dbsnp flag), then `AnalyzeCovariates` prints the plots/stats, then `PrintReads` writes the output bam
7. | **Sample Merging** | `bam_merge_samples_array.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_library_Lane_RX_Ye.recalibrated.bam` | `XXX_NNNN_all_libraries.bam` | I added this step so that we could use the GATK 'Joint Variant Calling' workflow; parallelizing over individuals. Merging files within individual, then re-index the resulting `..bam`
8. | **Variant Calling** | `vcf_discovery_array.qsub` | GATK/3.5.0 | `--genotyping_mode DISCOVERY --emitRefConfidence GVCF --output_mode EMIT_ALL_SITES` | `Sample_library_Lane_RX_Ye.recalibrated.bam` & `Reference.fa` | `XXX_NNNN_all_libraries.bam.g.vcf` | HaplotypeCaller...
9. | | `genotype_gvcf.qsub` | GATK/3.5.0 | `-stand_call_conf 30` & `-stand_emit_conf 30` | `*_all_libraries.bam` | `all_individuals_${date}.out.vcf` | This script uses `GenotypeGVCFs`for joint variant calling on the stack of all 63 `..bam.g.vcf` files
10. | | `vcf_merge_samples_array.qsub` | tabix/0.2.6 & vcftools/0.1.9 | `--remove-duplicates` | `Sample_library_Lane_RX_Ye.raw_variants.vcf.gz` | `Sample_all_libraries.vcf.gz` | vcf files merged at the level of population

11. | | `merge_all_vcfs.qsub` | tabix/0.2.6 & vcftools/0.1.9 | | `Sample_all_libraries.vcf.gz` | `all_variants_merged_$date$.vcf.gz` | this is a long step; merging all the sample-level vcfs into one file for passing to SKAT
12. | **Statistical Analysis** | `calc_Fst_array` | vcftools/0.1.9 | `--fst-window-size` & `--fst-window-step` | `all_variants_merged_${date}.vcf` | `Fst_POP1_vs_POP2.windowedXkb.stepYkb.weir.fst` | Fst stat. calculated for all pairwise between-pop comparisons. Options passed into output filenames.
13. | | `plink_vcf_convert.qsub`, `plink_prep.qsub`, & `plink_fisher.qsub` | vcftools/0.1.9 & plink/1.07 | `--maf` & `--geno` (missingness) | `all_variants_merged_${date}.vcf` | `all_variants_merged_${date}.assoc.fisher`, `*.ped`, `*.bed`, `*.bim`, `*.fam`, `*.log`, `*.map` & `*.nosex` | Association with Fisher's exact test and simplified presence/absence phenotype data.
14. | | script | R/3.2.0 & SKAT v1.1.2 | complex – see Rscript | `all_variants_merged_${date}.ped`, `*.bed`, `*.bim`, `*.fam`, & `*.map` | Sequence Kernal Association Test
15. | **Find Structural Variants** | `bam_merge_samples_array.qsub ` | picardTools/1.89 | | `*.bam` | `${individual}_all_libraries.bam` | ...
16. | | `breakdancer.qsub` | BreakDancer/1.1.2 & SAMTools/1.2 | `-t -q 10 -d` | `${individual}_all_libraries.bam` | `breakdancer_${date}_analysis.cfg` & `...ctx` | Detects structural variants...

---

## Version 4 (defunct)

---

## Version 5 (joint genotyping with HaplotypeCaller)

|  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
1. | **QC** | `trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and runs trimmomatic
2. | **Alignment** | `alignment_array.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Reference.fa` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | pulls details from fastq files to build a GATK-compatible `@RG` tag, and then runs alignment
3. | **Sorting/Indexing** | `deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
4. | | `indel_realign_array.qsub` | GATK/3.5.0 | | `Sample_library_Lane_RX_Ye.dedup.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.realignment_targets.list` & `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` | RealignerTargetCreator ID's targets, IndelRealigner realigns
5. | **Base Score Recalibration** | `base_score_recalibration_array.qsub` | GATK/3.5.0 | `--run_without_dbsnp_potentially_ruining_quality` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` & `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.recalibrated.bam` | 2 passes with `BaseRecalibrator` (with the infamous no-dbsnp flag), then `AnalyzeCovariates` prints the plots/stats, then `PrintReads` writes the output bam
6. | **Sample Merging** | `bam_merge_samples_array.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.recalibrated.bam` | `XXX_NNNN_all_libraries.bam` | I added this step so that we could use the GATK 'Joint Variant Calling' workflow; parallelizing over individuals. Merging files within individual, then re-index the resulting `..bam`
7. | | `merge_all_bams.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_all_libraries.bam` | `all_bam_merged_'+%d_%m_%Y'.bam` |
7. | | `merge_all_bams.qsub` | SAMTools/1.3.1 | `merge` | `XXX_NNNN_all_libraries.bam` | `all_bam_merged_'+%d_%m_%Y'.bam` |
8. | **Variant Calling** | `vcf_disco_chunk_array.qsub` | SAMTools/1.3.1, picardTools/1.89 & GATK/3.5.0 | `view -b`, `BuildBamIndex` & `HaplotypeCaller --genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES` | `all_bam_merged_'+%d_%m_%Y'.bam`, `indices.list` & `Reference.fa` | `all_bam_merged_'+%d_%m_%Y'.bam_chunkNNN.bam` | this script uses samtools to slice the all-fishes BAM file into chunks of ~1.5Mbp in length, passes that to picard in order to index it, then passes that to GATK-HaplotypeCaller
9. | | `reunite_vcf_chunks.qsub` | tabix/0.2.6 & vcftools/0.1.9 | `--remove-duplicates` |

---

## Version 6 (joint genotyping with HaplotypeCaller, indel realignment now deprecated, now with more paranoia)

|  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
0. | **QC** | `00_fastQC_array.qsub`  | FastQC/0.11.2 | | `Sample_library_Lane_RX_Ye.fastq.gz` | ??? | read quality assessment...
1. | **Trimming** | `01_trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and runs trimmomatic
2. | **Alignment** | `02_alignment_array.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Reference.fa` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | pulls details from fastq files to build a GATK-compatible `@RG` tag, and then runs alignment
3. | **Sorting/Indexing** | `03_deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
4. | **Base Score Recalibration** | `04_base_score_recalibration_array.qsub` | GATK/3.5.0 | `--run_without_dbsnp_potentially_ruining_quality` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` & `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | 2 passes with `BaseRecalibrator` (with the infamous no-dbsnp flag), then `AnalyzeCovariates` prints the plots/stats, then `PrintReads` writes the output bam
5. | **Sample Merging** | `05_bam_merge_samples_array.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | `XXX_NNNN_all_libraries.bam` | I added this step so that we could use the GATK 'Joint Variant Calling' workflow; parallelizing over individuals. Merging files within individual, then re-index the resulting `..bam`
6. | | `06_bam_file_check.qsub` | picardTools/1.89 | `ValidateSamFile` | `XXX_NNNN_all_libraries.bam` | `XXX_NNNN_samples`
7. | **Variant Calling** | `07_vcf_disco_chunk_6-1_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode EMIT_ALL_CONFIDENT_SITES` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | ```all_fish_`date '+%d_%m_%Y'`_slice_${n}.vcf``` | builds *n* 'slices' of `..vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`)
<!-- 8. | | `07_vcf_disco_chunk_6-2_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence GVCF` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | ```all_fish_`date '+%d_%m_%Y'`_slice_${n}.g.vcf``` | builds *n* 'slices' of `..g.vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`) -->

---


### File Flow
  1. 512 `..fastq.gz` files
  2. 1024 `..trimmed.fq` files
  3. 768 `..aligned.sam` files
  4. 768 `..aligned.sorted.sam` files
  5. 768 `..aligned.dedup.bam` files
  6. 768 `..dedup.realigned.recalibrated.bam` files
  7. 63 `..all_libraries.bam` files


### Tools & Version
  - [Trimmomatic/0.33]( http://www.usadellab.org/cms/?page=trimmomatic )
  - Java/1.8.0_31
  - [bwa/0.7.12.r1044]( http://bio-bwa.sourceforge.net/bwa.shtml )
  - picardTools/1.89
  - SAMTools/1.3.1
  - [GATK/3.7.0]( https://software.broadinstitute.org/gatk/ )
  - R/3.2.0
  - tabix/0.2.6
  - [vcftools/0.1.14]( https://vcftools.github.io/man_latest.html )
  - [plink/1.9]( https://www.cog-genomics.org/plink/1.9/ )


---


## Version 7 –– with **NEW**, deeper data, and including only Apassa & Bambomo fish

|  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
0. | **QC** | `00_fastQC_array.qsub`  | FastQC/0.11.2 | | `Sample_library_Lane_RX_Ye.fastq.gz` | ??? | read quality assessment...
1. | **Trimming** | `01_trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and runs trimmomatic
2. | **Alignment** | `02_alignment_array.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Reference.fa` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | pulls details from fastq files to build a GATK-compatible `@RG` tag, and then runs alignment
3. | **Sorting/Indexing** | `03_deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
4. | **Base Score Recalibration** | `04_base_score_recalibration_array.qsub` | GATK/3.5.0 | `--run_without_dbsnp_potentially_ruining_quality` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` & `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | 2 passes with `BaseRecalibrator` (with the infamous no-dbsnp flag), then `AnalyzeCovariates` prints the plots/stats, then `PrintReads` writes the output bam
5. | **Sample Merging** | `05_bam_merge_samples_array.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | `XXX_NNNN_all_libraries.bam` | I added this step so that we could use the GATK 'Joint Variant Calling' workflow; parallelizing over individuals. Merging files within individual, then re-index the resulting `..bam`
6. | | `06_bam_file_check.qsub` | picardTools/1.89 | `ValidateSamFile` | `XXX_NNNN_all_libraries.bam` | `XXX_NNNN_samples`
7. | **Variant Calling** | `07_vcf_disco_chunk_6-1_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode EMIT_ALL_CONFIDENT_SITES` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | ```all_fish_`date '+%d_%m_%Y'`_slice_${n}.vcf``` | builds *n* 'slices' of `..vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`)
<!-- 8. | | `07_vcf_disco_chunk_6-2_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode EMIT_ALL_CONFIDENT_SITES --emitRefConfidence GVCF` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | ```all_fish_`date '+%d_%m_%Y'`_slice_${n}.g.vcf``` | builds *n* 'slices' of `..g.vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`) -->


### File Flow
  1. 432 `..fastq.gz` files
    -
  2. 864 `..trimmed.fq` files
  3. 648 `..aligned.sam` files
  4. 648 `..aligned.sorted.sam` files
  5. 648 `..aligned.dedup.bam` files
  6. 648 `..dedup.realigned.recalibrated.bam` files
  7. 76 `..all_libraries.bam` files for 76 fish


### Tools & Version
  - [Trimmomatic/0.33]( http://www.usadellab.org/cms/?page=trimmomatic )
  - Java/1.8.0_31
  - [bwa/0.7.12.r1044]( http://bio-bwa.sourceforge.net/bwa.shtml )
  - picardTools/1.89
  - SAMTools/1.3.1
  - [GATK/3.7.0]( https://software.broadinstitute.org/gatk/ )
  - R/3.2.0
  - tabix/0.2.6
  - [vcftools/0.1.14]( https://vcftools.github.io/man_latest.html )
  - [plink/1.9]( https://www.cog-genomics.org/plink/1.9/ )



## Output Formats

  - **VCFtools ...weir.fst**
    1. CHROM      - Scaffold in our case
    2. BIN_START  - bp coordinate at start of window
    3. BIN_END    - bp coordinate at end of window
    4. N_VARIANTS - no. variants within the window
    5. WEIGHTED_FST
    6. MEAN_FST   - mean Fst across the window
  - **PLINK ..assoc.fisher**
    1. CHR - Scaffold in our case
    2. SNP - SNP ID
    3. BP - Physical position (base-pair)
    4. A1 - Minor allele name (based on whole sample)
    5. F_A - Frequency of this allele in cases
    6. F_U - Frequency of this allele in controls
    7. A2 - Major allele name
    8. P - Exact p-value for this test
    9. OR - Estimated odds ratio (for A1)
  - **BreakDancer**
    1. Chromosome               - Scaffold number in our case
    2. Position 1
    3. Orientation 1
    4. Chromosome 2
    5. Position 2
    6. Orientation 2
    7. Type of a SV             - DEL, INS, INV, ITX (intra-chromo. trans.), CTX (inter-chromo. trans.) or Unknown
    8. Size of a SV             - length of SV in bp (meaningless for CTX)
    9. Confidence Score
    10. No. supporting read pairs
    11. No. supporting read pairs from each map file
    12. Estimated allele freq.  - this estimate is not to be trusted at the current version
    13. Software version        -
    14. The run parameters


## Important Files to be preserved

  - Original read files - 512 x `${individual}_[GATC]_[GATC]_L00[1-8]_R[12]_00[12].fastq.gz`
  - First-pass variant calls for all individuals – `all_variants_merged_21_01_2016.vcf`
    - PLINK association results from the above `all_variants_merged_21_01_2016.assoc.fisher`
  - Second-pass variant calls for all individuals – `all_variants_merged_27_10_2015.vcf`
    - PLINK association results from the above `all_variants_merged_27_10_2015.assoc.fisher`
  -



  ===

  ## Version 8 –– with **NEW** data, including only Apassa & Bambomo fish, with 2 forks, and with extra paranoid self-documentation

|  Stage   |   Script    |   Tool(s)   |   Options   |   Input   |   Output    |   Description
---|----------|-------------|-------------|-------------|-----------|-------------|---------------
0. | **QC** | `00_fastQC_array.qsub`  | FastQC/0.11.2 | | `Sample_library_Lane_RX_Ye.fastq.gz` | ??? | read quality assessment...
1. | **Trimming** | `01_trimmomatic_array.qsub ` | Trimmomatic/0.32 | `ILLUMINACLIP:Illumina_adapters.fa:2:30:10 HEADCROP:10 MAXINFO:50:0.5` | `Sample_library_Lane_RX_Ye.fastq.gz` | `Sample_library_Lane_RX_Ye.trimmed.fq` | gunzips ..fastq.gz files and runs trimmomatic
2. | **Alignment** | `02_1_align_array_IL.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `supercontigs.fasta` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | pulls details from fastq files to build a GATK-compatible `@RG` tag, and then runs alignment against the P.kings genome
  | | `02_2_align_array_BN.qsub` | bwa/0.7.12.r1044 | `bwa mem -M -R ..` | `Sample_library_Lane_RX_Ye.trimmed.fq` & `Para_king_2015_013_20_40_15_90_3_superscaffold.fasta` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | pulls details from fastq files to build a GATK-compatible `@RG` tag, and then runs alignment against the *bionano-scaffolded* version of the P.kings genome
3. | **Sorting/Indexing** | `03_deduplication_array.qsub` | picardTools/1.89 | `SORT_ORDER=coordinate` | `Sample_library_Lane_RX_Ye.trimmed.aligned.sam` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.bam` | SortSam.jar sorts, MarkDuplicates.jar marks duplicates & BuildBamIndex.jar indexes
  |

4. | **Base Score Recalibration** | `04_base_score_recal_array1.qsub` | GATK/3.5.0 | `--run_without_dbsnp_potentially_ruining_quality` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_data.table` | 2 passes with `BaseRecalibrator` (with the infamous no-dbsnp flag), this writes recal. tables
5. | **Base Score Recalibration** | `05_base_score_recal_array2.qsub` | GATK/3.5.0 | `--run_without_dbsnp_potentially_ruining_quality` | `Sample_library_Lane_RX_Ye.trimmed.aligned.dedup.realigned.bam` & `Reference.fa` | `Sample_library_Lane_RX_Ye.recal_plots.pdf` & `Sample_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | `AnalyzeCovariates` prints the plots/stats, then `PrintReads` writes the output bam based on the recal. table
6. | **Sample Merging** | `05_bam_merge_samples_array.qsub` | picardTools/1.89 | `MergeSamFiles.jar` & `BuildBamIndex` | `XXX_NNNN_library_Lane_RX_Ye.aligned.dedup.realigned.recalibrated.bam` | `XXX_NNNN_all_libraries.bam` | I added this step so that we could use the GATK 'Joint Variant Calling' workflow; parallelizing over individuals. Merging files within individual, then re-index the resulting `..bam`
7. | | `06_bam_file_check.qsub` | picardTools/1.89 | `ValidateSamFile` | `XXX_NNNN_all_libraries.bam` | `XXX_NNNN_samples`

8. | **Variant Calling** | `07_vcf_disco_chunk_6-1_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode EMIT_ALL_CONFIDENT_SITES` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | ```all_fish_`date '+%d_%m_%Y'`_slice_${n}.vcf``` | builds *n* 'slices' of `..vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`)

8. | **Variant Calling** | `08_vcf_disco_chunk_8_array.qsub` | GATK/3.7.0 `-T HaplotypeCaller` | `--genotyping_mode DISCOVERY -stand_call_conf 30  -mbq 20 --output_mode GVCF` | 63 x`POP_ID##_all_libraries.bam` & `indices.list` | 63 x`POP_ID##_all_libraries.g.vcf` |
??? builds *n* 'slices' of `..vcf` file, where slice size is set by the file `indices.list` (itself build by `write_scaf_indices.sh`)


  ### File Flow
    0. 344  `..fastq.gz` files
    1. 688  `..trimmed.fq` files
    2. 516  `..aligned.sam` files
    3. 516  `..aligned.sorted.sam` & `..aligned.dedup.sam` files
    4. xxx  `..dedup.realigned.recal_data.table` & `..dedup.realigned.post_recal_data.table` files
    5. xxx  `..dedup.realigned.recalibrated.bam` & `..dedup.realigned.recal_plots.pdf` files
    6. 62   `..all_libraries.bam` files for 62 fish
    7. xxx  
    8. xxx  

===

  ----

  explicit list of filters on way in to PLINK:

    - GQx -
    - filterExpression
      - QD < 2.0
      - FS > 60.0
      - MQ < 40.0
      - MQRankSum < -12.5
      - ReadPosRankSum < -8.0
    -
