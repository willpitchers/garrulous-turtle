# script to read in a folder of ..fastqc_data.txt files and plot their quality score by base pair

# these lines grab a list of the fastqc_data files in the directory and the name of the directory
files.list <- dir( pattern=".+fastqc\\_data\\.txt" )
nfiles <- length( files.list )
dirname <- unlist( strsplit( getwd(), "/" ))[ length( unlist( strsplit( getwd(), "/" ))) ]

# open the pdf device to plot
pdf( paste( dirname, "GC_content.pdf", sep='', coll='' ), width=7, height=7 )
par( las=2 )

# the first loop draws the plot of *median* quality score by base pair
for ( i in 1:nfiles ){
  lines <- readLines( files.list[i] )
  headrows <- grep( ">>Per base sequence content\t[warnpsfil]{4}", lines )+1
  tailrows <- grep( ">>END_MODULE", lines)[grep( ">>END_MODULE", lines) > headrows][1]-1

  basegc <- gsub( "[- #]", "_", lines[ (headrows+1):tailrows ] )
  basegcnames <- gsub( "[- #]", "_", lines[ headrows ] )
  base.gc <- data.frame( matrix( unlist( strsplit( basegc, '\t' )), ncol=5, byrow=T ))

    names( base.gc ) <- unlist( strsplit( basegcnames, '\t' ))
    for (j in 2:5) base.gc[,j] <- as.numeric( as.character( base.gc[,j]) )
    base.gc$`_Base` <- reorder( factor( gsub( "^", "b_", base.gc$`_Base` ) ), 1:67 )

  if (i == 1) {
  plot( as.numeric(base.gc$`_Base`), base.gc$G, ylim=c(10, 40), cex.axis=0.5, pch=16, type="n", xlab="position (bp)", xaxt='n', ylab="Q score", main="GC content" )
    axis( side=1, at=1:67, labels=base.gc$`_Base`, cex.axis=0.5 )
        points( base.gc$`_Base`, base.gc$G, type='b', cex=0.5, col="green", pch=16 )
        points( base.gc$`_Base`, base.gc$C, type='b', cex=0.5, col="black", pch=16 )
        points( base.gc$`_Base`, base.gc$T, type='b', cex=0.5, col="blue", pch=16 )
        points( base.gc$`_Base`, base.gc$A, type='b', cex=0.5, col="red", pch=16 )
      } else {
      points( base.gc$`_Base`, base.gc$G, type='b', cex=0.5, col="green", pch=16 )
      points( base.gc$`_Base`, base.gc$C, type='b', cex=0.5, col="black", pch=16 )
      points( base.gc$`_Base`, base.gc$T, type='b', cex=0.5, col="blue", pch=16 )
      points( base.gc$`_Base`, base.gc$A, type='b', cex=0.5, col="red", pch=16 )
      }
      print( i ) }

# close device
dev.off()
