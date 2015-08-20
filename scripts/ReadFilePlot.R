setwd( "/Volumes/efish/2015_genomic_data/" )

system( "ls -latr 20150707_DNASeq_PE/*gz > 07_list.txt" )
system( "ls -latr 20150721_DNASeq_PE/*gz > 21_list.txt" )

dir()

head( read.table( "07_list.txt" ) )
reads07 <- read.table( "07_list.txt" )
reads21 <- read.table( "21_list.txt" )

names(reads21)
reads07$ID <- substring( as.character( reads07$V9 ), 1, 10 )
reads21$ID <- substring( as.character( reads21$V9 ), 1, 10 )
head(reads07)

readsboth <- merge( reads21, reads07, by="ID" )
head(readsboth)

str( readsboth )
readsboth$ID <- factor( readsboth$ID )

plot( readsboth$V5.x, readsboth$V5.y )

summary(c( readsboth$V5.x, readsboth$V5.y )/1e+9)

plot( readsboth$V5.x/1e+9, readsboth$V5.y/1e+9, xlim=c(0.15, 1.5), ylim=c(0.15, 1.5), ylab="20150707", xlab="20150721", col=c(1:64)[ readsboth$ID ], main="read file size (GB)" )
			text( 1.2, 0.6, "colour by barcode" )

system( "rm 07_list.txt 21_list.txt" )
