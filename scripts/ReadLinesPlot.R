setwd( "/Volumes/efish/2015_genomic_data/" )

system( "wc -l /mnt/scratch/pitchers/eFISH/20150721_DNASeq_PE/*.fastq > /mnt/research/efish/2015_genomic_data/21_list.txt" )
system( "wc -l /mnt/scratch/pitchers/eFISH/20150707_DNASeq_PE/*.fastq > /mnt/research/efish/2015_genomic_data/07_list.txt" )

dir()

head( read.table( "07_list.txt" ) )
reads07 <- read.table( "07_list.txt" )
reads21 <- read.table( "21_list.txt" )

names(reads21)
reads07$ID <- substring( as.character( reads07$V2 ), 1, 10 )
reads21$ID <- substring( as.character( reads21$V2 ), 1, 10 )
head(reads07)

readsboth <- merge( reads21, reads07, by="ID" )
head(readsboth)

str( readsboth )
readsboth$ID <- factor( readsboth$ID )

summary(c( readsboth$V1.x, readsboth$V1.y )/1e+06)

plot( readsboth$V1.x/1e+6, readsboth$V1.y/1e+6, xlim=c(6, 65), ylim=c(6, 65), ylab="20150707", xlab="20150721", col=c(1:64)[ readsboth$ID ], main="no. mega-rows in fastq" )
			text( 40, 30, "colour by barcode" )

system( "rm 07_list.txt 21_list.txt" )
