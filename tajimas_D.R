
require( ggplot2 )
require( dplyr )

setwd( "/mnt/lustre_scratch_2012/pitchers/eFISH/Illumina_2015/GATK-compatible" )

taj <- read.table( "all_pops.Tajima.D", header=TRUE )

tajD <- filter( taj, TajimaD!="NaN" ) %>% 
  group_by( CHROM ) %>% 
	summarise( mean_nSNPS=mean( N_SNPS ), mean_TajD=mean( TajimaD ), scaf_len_kb=(length( CHROM )*100) )

pdf( "TajD.pdf" )
	mutate( tajD, scaffold_length=scaf_len_kb ) %>%
	ggplot( aes( mean_nSNPS, mean_TajD )) +
		geom_point( aes(color=scaffold_length ), cex=3, pch=16 ) + 
		xlab( "mean no. SNPs per bin") + 
		ylab( "mean Tajima's D per bin" )
dev.off( )

write.csv( tajD, "tajimasD_by_scaffold.csv", quote=FALSE, row.names=FALSE )

