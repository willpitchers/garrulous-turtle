#!/bin/env python

### This script is written to be called by `plink_ind_pheno`.
### It takes 2 input files ( `individual.list` & `all_variants_merged_${date}.ped` ),
### and writes out a new `.ped` file with the phenotype codes specified in the
### `individual.list` file written in to the 6th field.

## scrape the arguments passed to the script
import sys
mapfile = sys.argv[1]
pedfile = sys.argv[2]
outfile = sys.argv[3]

## this function turns the mapfile into a dictionary
def write_map( mapfile ):
    pheno_dict = {}
    for line in open( mapfile, 'r' ):
        elements = line.strip().split()
        ID, pheno = elements[ 0 ], elements[ 1 ]
        pheno_dict[ ID ] = pheno
    return pheno_dict

## this function queries the dictionary
def return_pheno( ID_no, pheno_dict ):
    pheno = pheno_dict[ ID_no ]
    return pheno

## this function reads the .ped file, queries the dictionary,
## and writes the phenotype code to the output file
def read_ped( pedfile, outfile ):
    outped = open( outfile, 'w' )
    for line in open( pedfile, 'r' ):
        Line = line.strip().split()
        ID = line.strip().split()[ 0 ].split( '.' )[ 0 ]
        pheno = return_pheno( ID, pheno_dict )
        Line[ 6 ] = pheno
        Lineout = ' '.join( Line )
        outped.write( Lineout + '\n' )

## these lines run the input through the functions
pheno_dict = write_map( mapfile )
read_ped( pedfile, outfile )
