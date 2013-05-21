library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(splines)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
chip<-args[1]
input<-args[2]
fname<-args[3]
path<-args[4]

ChIPQC <- function(rsChip, rsInput, name, path, windowSize=1000, dataPoints=1000)
{	
	hg19.windows <- genomeBlocks(Hsapiens, chrs=1:24, windowSize)
	ChIPcounts <- countOverlaps(hg19.windows, rsChIP)
	INPUTcounts <- countOverlaps(hg19.windows, rsINPUT)
	
	hg19.ordered = data.frame(cbind(ChIPcounts, INPUTcounts)[order(ChIPcounts),])
	
	hg19.ordered$ChIPcounts = cumsum(hg19.ordered$ChIPcounts)
	hg19.ordered$INPUTcounts = cumsum(hg19.ordered$INPUTcounts)
	
	
	hg19.ordered$ChIPcounts = hg19.ordered$ChIPcounts/hg19.ordered$ChIPcounts[length(hg19.windows)]
	hg19.ordered$INPUTcounts = hg19.ordered$INPUTcounts/hg19.ordered$INPUTcounts[length(hg19.windows)]	

	spaced <- round(seq(1,length(hg19.windows), by=length(hg19.windows)/dataPoints))
	df1 <- data.frame("bin"=c(1:dataPoints),"value"=hg19.ordered$ChIPcounts[spaced], "type"="ChIP")
	df2 <- data.frame("bin"=c(1:dataPoints),"value"=hg19.ordered$INPUTcounts[spaced], "type"="INPUT")
	df <- data.frame(rbind(df1,df2))
	maxdist <- which.max(abs(hg19.ordered$ChIPcounts[spaced]-hg19.ordered$INPUTcounts[spaced]))
	
	p1 <-ggplot(df, aes(x=bin, y=value, group=type)) + geom_line(aes(color=type)) 
	p1 <- p1 + geom_vline(xintercept = maxdist, colour="grey")
	p1 <- p1 + xlab("Percentage of bins") + ylab("Percentage of reads") + ggtitle(name)
	p1 <- p1 + theme(legend.position = c(0.1, 0.9), legend.background = element_rect(fill = "white", colour = NA))	
	
	
	# get derivative
	x1 = seq(1, dataPoints-2)
	x2 = seq(3, dataPoints)
	slope1 = 1/2*(hg19.ordered$ChIPcounts[spaced[x2]]-hg19.ordered$ChIPcounts[spaced[x1]])
	slope2 = 1/2*(hg19.ordered$INPUTcounts[spaced[x2]]-hg19.ordered$INPUTcounts[spaced[x1]])
	
	df3 <- data.frame(bin=c(1: length(slope1)), value=slope1, type="ChIP")
	df4 <- data.frame(bin=c(1: length(slope1)), value=slope2, type="INPUT")
	df5 <- data.frame(rbind(df3,df4))
	
	p2 <-ggplot(df5, aes(x=bin, y=value, group=type)) + geom_line(aes(color=type)) + ggtitle("Derivative") + xlab("Percentage of bins") + ylab("dy/dx")
	p2 <- p2 + geom_smooth(method = "lm",formula = y~bs(x, degree = 6),se = TRUE, alpha=0.5)
	p2 <- p2 + geom_vline(xintercept = maxdist, colour="grey")
	p2 <- p2 + theme(legend.position="none")
	
	# print plots
	pdf(paste(path,"/",name, ".pdf", sep=""), width = 10, height=5)
	grid.arrange(p1, p2 , ncol=2)
	dev.off()
}

rsChIP <- BAM2GRanges(chip)
rsINPUT <- BAM2GRanges(input)
#only count read starts, avoids reads being across multiple windows
rsChIP <- resize(rsChIP, 1, fix="start")
rsINPUT <- resize(rsINPUT, 1, fix="start")

ChIPQC(rsChip,rsInput, fname,path)
