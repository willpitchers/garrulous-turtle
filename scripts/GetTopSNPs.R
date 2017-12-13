#!/usr/bin/env Rscript

args <- commandArgs( trailingOnly=TRUE )
require( dplyr )
require( data.table )

R.Version()
####

filename <- args[ 1 ]

assoc <- tbl_df( fread( filename,  header=TRUE ))

if ( grepl( "model", filename ) == TRUE ) {
		assoc %>% arrange( P ) %>% mutate( BP=as.numeric( gsub( "Var-Scaffold[0-9]+-", "", SNP ))) %>% head( 10000 ) %>%
				write.csv( paste( filename, ".topsnps.csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )
	}


#if ( grepl( "assoc.logistic", filename ) == TRUE ) {
#		perms <- tbl_df( fread( paste( filename, ".mperm", sep='', coll='' )))
#		full_join( assoc, perms, by=c( "SNP", "SNP" )) %>% 
#		arrange( EMP2 ) %>%
#		arrange( P ) %>%
#				head( 5000 ) %>% write.csv( paste( filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )
#	}


if ( grepl( "fisher", filename ) == TRUE ) {
		assoc <- transmute( assoc, Scaf=factor( CHR ), SNP=BP, A1=factor( A1 ), F_A=F_A, F_U=F_U, A2=factor( A2 ), P=P, OR=OR, SE=SE )
		assoc %>% arrange( P ) %>% head( 10000 ) %>% write.csv( paste( filename, ".topsnps.csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )
	}

q( 'no' )
