#!/bin/bash
# All possible pairwise comparisons between populations
#
# APA_vs_BAM  APA_vs_BAVA APA_vs_COB  APA_vs_IVI APA_vs_MOV 
# BAM_vs_BAVA BAM_vs_COB  BAM_vs_IVI  BAM_vs_MOV 
# BAVA_vs_COB BAVA_vs_IVI BAVA_vs_MOV 
# COB_vs_IVI  COB_vs_MOV 
# IVI_vs_MOV

# setting 2 variables thus:
export leftpop=( APA BAM BAVA COB IVI APA BAM BAVA COB APA BAM BAVA APA BAM APA )
export rightpop=( BAM BAVA COB IVI MOV BAVA COB IVI MOV COB IVI MOV IVI MOV MOV )

# ...will allow for comparisons to be called with:
# ..population_${leftpop[${n}]}.txt ..population_${rightpop[${n}]}.txt

