
rm( list=ls())
require( dplyr )

geno <- read.table( "all_individuals_12_07_16.GT.FORMAT.recode", header=TRUE )

keepCols <- function( data, cols ){
	return( data[ , names(data) %in% cols] )
	}

pheno1 <- c( "X4004", "X4005", "X4006", "X4027", "X4029", "X4030", "X6496", "X6498", "X6499", "X6501", "X6502", "X6598", "X6599", "X6602", "X6603", "X6604", "X6605", "X6619", "X6620", "X6621", "X6622", "X6623", "X6624", "X6625", "X6626", "X6627", "X6681" )

pheno2 <- c( "X4018", "X4019", "X4816", "X4832", "X4834", "X4893", "X4894", "X4895", "X4896", "X4897", "X4921", "X3923", "X4925", "X6494", "X6497", "X6500", "X6597", "X6489", "X6675", "X6676", "X6677", "X6678", "X6679", "X6680", "X6682", "X6683", "X6684", "X6685", "X6716", "X6717", "X6718", "X6719", "X6720", "X6721", "X6722", "X6723", "X6724", "X6725", "X6737" )

geno <- mutate( geno,  ct_NA_biphase= rowSums( is.na( keepCols( geno, pheno1 ))), 
                  ct_NA_triphase= rowSums( is.na( keepCols( geno, pheno2 ))) )


# the contingency table:
#		control		case
#		biphasic	triphasic
# ----------------------------
# ref		A			C
# A1
# alt		B			D
#A2

numSNPs <- nrow( geno )
taco <- rep(NA, numSNPs)

for( i in 1:numSNPs ) {
	maxBi <- 50 - 2*geno$ct_NA_biphase[i]
	maxTri <- 76 - 2*geno$ct_NA_triphase[i]
    # biphasic, reference. Max possible 50 alleles
	A <- maxBi - sum( keepCols( geno[ i, ], pheno1 ), na.rm=TRUE )
	# biphasic, alternate. Max possible 50 alleles
	B <- sum( keepCols( geno[ i, ], pheno1 ), na.rm=TRUE )
	# triphasic, reference. Max possible 76 alleles
	C <- maxTri - sum( keepCols( geno[ i, ], pheno2 ), na.rm=TRUE )
	# triphasic, alternate.
	D <- sum( keepCols( geno[ i, ], pheno2 ), na.rm=TRUE )

	ConTab <- matrix( c(A,B,C,D), ncol=2 )
	FishTest <- fisher.test( ConTab )
	taco[i] <- FishTest$p.value
}


