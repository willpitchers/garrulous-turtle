#!/usr/bin/env Rscript

args <- commandArgs( trailingOnly=TRUE )

require( MASS )
require( dplyr )
require( data.table )

# R.Version()

####

filename <- args[ 1 ]

assoc <- tbl_df( fread( filename,  header=TRUE ))


if ( grepl( "fisher", filename ) == TRUE ) {
		assoc <- transmute( assoc, Scaf=factor( CHR ), SNP=BP, A1=factor( A1 ), F_A=F_A, F_U=F_U, 
									A2=factor( A2 ), P=P, OR=OR, SE=SE )
		assoc %>% arrange( P ) %>% 
		head( 250 ) %>% 
		write.csv( paste( filename, ".csv", coll='', sep='' ), quote=FALSE, row.names=FALSE )
	} 


if ( grepl( "model", filename ) == TRUE ) {
		assoc %>% arrange( P ) %>% mutate( SNP=gsub( "Var-Scaffold[0-9]+-", "", SNP )) %>% 
		head( 250 ) %>% 
		write.csv( paste( "top_250_", filename, coll='', sep='' ), quote=FALSE, row.names=FALSE )
	}


if ( grepl( "mperm", filename ) == TRUE ) {
		perms <- tbl_df( fread( paste( filename, ".mperm", sep='', coll='' )))
		assoc <- full_join( assoc, perms, by="SNP" )
		assoc %>% arrange( EMP2 ) %>% head( 250 ) %>% 
				 write.csv( paste( "top_250_", filename, coll='', sep='' ), quote=FALSE, row.names=FALSE )
	}


q( 'no' )
