#!/bin/env python

### This script is intended to be called with a stripped-down VCF file produced
### e.g. by: `grep '^Scaffold' myfile.vcf | cut -d' ' -f 1-2`

import sys, subprocess

mycommand = "wc -l SNPlist"
process = subprocess.Popen( mycommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

# pass these in as variables where possible
numSNPs=int( output.split()[0] )
winLen=1000
windows=range( 0, numSNPs, winLen )

def write_SetID_file( infile, outfile ):
	inFile = open( infile, 'r' )
	outFile = open( outfile, 'w' )
	for line in inFile:
		Line = line.strip().split()
		winID = GetWindowID( Line[0], Line[1] )
		Line = Line.append( winID )
		Lineout = '\t'.join( Line )
		outFile.write( Lineout + '\n' )
	inFile.close()
	outFile.close()


def GetWindowID( scaf, snp ):
	# scaflength = D[ scaf ]
	# windows = range( 0, scaflength, winLen )
	# numWin = len( windows )
	# IDs = {}
	i, j = 0, 0
	while i < snp:
		i = i + winLen
		j += 1
	return ''.join(( "win_", j ))


def BuildScafDic( scaf_file ):
	D = {}
	with open( scaf_file ) as F:
		for line in F:
			(key, val) = line.strip.split()
			D[ key ] = val
	return D

########

# run the shit here
