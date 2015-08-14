# script to read in a folder of ..fastqc_data.txt files and plot their quality score by base pair

# these lines grab a list of the fastqc_data files in the directory and the name of the directory
files.list <- dir( pattern=".+fastqc\\_data\\.txt" )
nfiles <- length( files.list )
dirname <- unlist( strsplit( getwd(), "/" ))[ length( unlist( strsplit( getwd(), "/" ))) ]

# open the pdf device to plot
pdf( paste( dirname, "Qscores.pdf", sep='', coll='' ), width=14, height=7 )
par( mfrow=c( 1, 2 ), las=2 )

# the first loop draws the plot of *median* quality score by base pair
for ( i in 1:nfiles ){
  lines <- readLines( files.list[i] )
  headrows <- grep( ">>Per base sequence quality\tpass", lines )+1
  tailrows <- grep( ">>END_MODULE", lines)[grep( ">>END_MODULE", lines) > headrows][1]-1

  basequal <- gsub( "[- #]", "_", lines[ (headrows+1):tailrows ] )
  basequalnames <- gsub( "[- #]", "_", lines[ headrows ] )
  base.qual <- data.frame( matrix( unlist( strsplit( basequal, '\t' )), ncol=7, byrow=T ))

    names( base.qual ) <- unlist( strsplit( basequalnames, '\t' ))
    for (j in 2:7) base.qual[,j] <- as.numeric( as.character( base.qual[,j]) )
    base.qual$`_Base` <- reorder( factor( gsub( "^", "b_", base.qual$`_Base` ) ), 1:67 )

  if (i == 1) {
  plot( as.numeric(base.qual$`_Base`), base.qual$Mean, ylim=c(20, 40), cex.axis=0.5, pch=16, type="n", xlab="position (bp)", xaxt='n', ylab="Q score", main="Median" )
    axis( side=1, at=1:67, labels=base.qual$`_Base`, cex.axis=0.5 )
      abline( h=28, lty=2 )
        points( base.qual$`_Base`, base.qual$Median, type='b', cex=0.5, col=i )
      } else {
      points( base.qual$`_Base`, base.qual$Median, type='b', cex=0.5, col=i ) } }

# the second loop draws the plot of *mean* quality score by base pair
for (i in 1:nfiles ){
  lines <- readLines( files.list[i] )
  headrows <- grep( ">>Per base sequence quality\tpass", lines )+1
  tailrows <- grep( ">>END_MODULE", lines)[grep( ">>END_MODULE", lines) > headrows][1]-1

  basequal <- gsub( "[- #]", "_", lines[ (headrows+1):tailrows ] )
  basequalnames <- gsub( "[- #]", "_", lines[ headrows ] )
  base.qual <- data.frame( matrix( unlist( strsplit( basequal, '\t' )), ncol=7, byrow=T ))

    names( base.qual ) <- unlist( strsplit( basequalnames, '\t' ))
    for (j in 2:7) base.qual[,j] <- as.numeric( as.character( base.qual[,j]) )
    base.qual$`_Base` <- reorder( factor( gsub( "^", "b_", base.qual$`_Base` ) ), 1:67 )

  if (i == 1) {
  plot( as.numeric(base.qual$`_Base`), base.qual$Mean, ylim=c(20, 40), cex.axis=0.5, pch=16, type="n", xlab="position (bp)", xaxt='n', ylab="Q score", main="Mean" )
      axis( side=1, at=1:67, labels=base.qual$`_Base`, cex.axis=0.5 )
        abline( h=28, lty=2 )
          points( base.qual$`_Base`, base.qual$Mean, type='b', cex=0.5 )
      } else {
      points( base.qual$`_Base`, base.qual$mean, type='b', cex=0.5, col=i ) } }

# close device
dev.off()
