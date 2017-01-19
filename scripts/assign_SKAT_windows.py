#!/bin/env python

### This script is intended to be called with a stripped-down VCF file produced
### e.g. by: `grep '^Scaffold' myfile.vcf | cut -d' ' -f 1-2`

import sys, subprocess

scaf_file = sys.argv[1]
infile = sys.argv[2]
outfile = sys.argv[3]

mycommand = ''.join(( "wc -l ", infile ))
process = subprocess.Popen( mycommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

# pass these in as variables where possible
numSNPs=int( output.split()[0] )
winLen=10000
windows=range( 0, numSNPs, winLen )


def write_SetID_file( infile, outfile ):
	inFile = open( infile, 'r' )
	outFile = open( outfile, 'w' )
	for line in inFile:
		Line = line.strip().split()
		winID = GetWindowID( Line[0], Line[1] )
		LineW = [ winID ] + Line
		Lineout = '\t'.join( LineW )
		outFile.write( Lineout + '\n' )
	inFile.close()
	outFile.close()


def GetWindowID( scaf, snp ):
	scaflength = ScafDict[ scaf ]
	# windows = range( 0, scaflength, winLen )
	# numWin = len( windows )
	# IDs = {}
	i, j = 0, 0
	while i < int( snp ):
		i = i + winLen
		j += 1
	return ''.join(( "win_", str( j ) ))


def BuildScafDic( scaf_file ):
	D = {}
	with open( scaf_file ) as F:
		for line in F:
			(key, val) = line.strip().split()
			D[ key ] = val
	return D

########

# run the shit here

ScafDict = BuildScafDic( scaf_file )

write_SetID_file( infile, outfile )
