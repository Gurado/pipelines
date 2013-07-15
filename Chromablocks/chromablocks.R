library(GenomicFeatures)
library(Repitools)
library(chipseq)
library(BSgenome.Hsapiens.UCSC.hg19)
library(multicore)
library(rtracklayer)
options(cores=4)

# expects two bam files (signal, input) and run chromablocks 
# expects HG19 assembly

args <- commandArgs(trailingOnly = TRUE)

chip_arg<-args[1]
input_arg<-args[2]
output_arg<-args[3]

chromacall <- function(signal, input, output) 
{
	ip_name <- sub("[.][^.]*$", "", basename(signal))
	in_name <- sub("[.][^.]*$", "", basename(input))
	rs<-BAM2GRangesList(c("ip"=signal, "input"=input))

	chroma.sm <- ChromaBlocks(rs.ip=rs[1], rs.input=rs[2], Hsapiens, 1:24, preset="small")
        save(chroma.sm,file=file.path(output, paste(paste(ip_name,in_name,sep='-'),"sm","rda",sep='.')))	
	chroma.sm.regions <- regions(chroma.sm)
	export(chroma.sm.regions, con=file.path(output,paste(paste(ip_name,in_name,sep='-'),"sm","bed",sep='.')), format="bed")

	chroma.lg <- ChromaBlocks(rs.ip=rs[1], rs.input=rs[2], Hsapiens, 1:24, preset="large")
        save(chroma.lg,file=file.path(output, paste(paste(ip_name,in_name,sep='-'),"lg","rda",sep='.')))
	chroma.lg.regions <- regions(chroma.lg)
	export(chroma.lg.regions, con=file.path(output,paste(paste(ip_name,in_name,sep='-'),"lg","bed",sep='.')), format="bed")

}

chromacall(signal=chip_arg, input=input_arg, output=output_arg)
