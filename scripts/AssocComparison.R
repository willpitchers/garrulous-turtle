rm( list=ls())

require( MASS )
require( fdrtool )
require( ggplot2 )
require( dplyr )
require( VennDiagram )
require( cowplot )
require( scales )
require( data.table )

R.Version()

filename <- c( "all_fish_version_5-1_HPC_geno50.assoc" )

assoc <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc, Scaf=factor( gsub( "0_(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "0_Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
                            A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), P=P, OR=OR )

assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "fisher_", "all_fish_version_5-1_HPC_geno50.assoc", ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )

# assoc$QVAL <- fdrtool( assoc$P, statistic="pvalue", plot=FALSE )$qval
# assoc$LFDR <- fdrtool( assoc$P, statistic="pvalue", plot=FALSE )$lfdr


# other modelled effects

filenames <- c( "all_fish_version_5-1_HPC_g50_CI.GENO.model.reformatted.csv", "all_fish_version_5-1_HPC_g50_CI.ALLELIC.model.reformatted.csv", "all_fish_version_5-1_HPC_g50_CI.DOM.model.reformatted.csv", "all_fish_version_5-1_HPC_g50_CI.REC.model.reformatted.csv", "all_fish_version_5-1_HPC_g50_CI.TREND.model.reformatted.csv" )

for( i in 1:length( filenames )) {
		filename <- filenames[ i ]
		assoc <- tbl_df( fread( filename,  header=TRUE ))
		assoc <- mutate( assoc, CHR=factor( gsub( "0_(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )), 
				                    SNP=as.integer( gsub( "0_Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )))
		assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "model_", filename, coll='', sep='' ), quote=FALSE, row.names=FALSE )
	}

###

filename <- c( "without_COB_5-1_HPC.2017_geno50.assoc" )

assoc <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc, Scaf=factor( gsub( "0_(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "0_Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
                            A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), P=P, OR=OR )

assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "fisher_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )

# assoc$QVAL <- fdrtool( assoc$P, statistic="pvalue", plot=FALSE )$qval
# assoc$LFDR <- fdrtool( assoc$P, statistic="pvalue", plot=FALSE )$lfdr


# other modelled effects

filenames <- c( "without_COB_5-1_HPC.2017_g50_CI.GENO.model.reformatted.csv", "without_COB_5-1_HPC.2017_g50_CI.ALLELIC.model.reformatted.csv",  "without_COB_5-1_HPC.2017_g50_CI.DOM.model.reformatted.csv", "without_COB_5-1_HPC.2017_g50_CI.REC.model.reformatted.csv", "without_COB_5-1_HPC.2017_g50_CI.TREND.model.reformatted.csv" )

for( i in 1:length( filenames )) {
        filename <- filenames[ i ]
        assoc <- tbl_df( fread( filename,  header=TRUE ))
        assoc <- mutate( assoc, CHR=factor( gsub( "0_(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                                    SNP=as.integer( gsub( "0_Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )))
        assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "model_", filename, coll='', sep='' ), quote=FALSE, row.names=FALSE )
    }


##### 3_level phenotype

filename <- c( "all_fish_version_5-1_HPC_3level.assoc" )

assoc3 <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc3, Scaf=factor( gsub( "(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
                            A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), P=P, OR=OR )


assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "3level_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )


##### 3-level phenotype with pop as covariate

filename <- c( "all_fish_version_5-1_HPC_3level_PopCov.assoc.logistic" )

assoc4 <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc4, Scaf=factor( gsub( "(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
                            A1=factor( A1 ), Test=factor( TEST ), N=NMISS, P=P, OR=OR, Tstat=STAT )


assoc %>% arrange( P ) %>% filter( Test=="ADD" ) %>% head( 250 ) %>% write.csv( paste( "PopCov_P0_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )

assoc %>% arrange( P ) %>% filter( Test=="POP" ) %>% head( 250 ) %>% write.csv( paste( "PopCov_pop_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )


##### 3-level phenotype with pop as covariate

filename <- c( "all_fish_version_5-1_HPC_3levelas2.assoc" )

assoc6 <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc6, Scaf=factor( gsub( "(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
							A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), chisq=CHISQ, P=P, OR=OR )


assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "3levelsas2_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )


##### 2-level phenotype stratified by pop.

filename <- c( "all_fish_version_5-1_HPC_strata.assoc" )

assoc5 <- tbl_df( fread( filename,  header=TRUE ))

assoc <- transmute( assoc5, Scaf=factor( gsub( "0_(Scaffold[0-9]{1,4}):[0-9]+", "\\1", SNP )),
                            SNP=as.integer( gsub( "0_Scaffold[0-9]{1,4}:([0-9]+)", "\\1", SNP )),
                            A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), chisq=CHISQ, P=P, OR=OR )


assoc %>% arrange( P ) %>% head( 250 ) %>% write.csv( paste( "Strata_", filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )


q( "no" )
