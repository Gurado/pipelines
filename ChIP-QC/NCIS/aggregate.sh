#!/bin/sh -e

echo "ChIP	Control	Normalization	Binsize	SeqDepth	BackgroundRatio" > NCIS_aggregate.txt

for f in `ls results/*.txt`; do

	b=$(basename $f)
        t=(${b//./ })
	chip=(${t%-*})
	input=(${t##*-})

	tail -n 1 $f | awk -v i1=${chip} -v i2=${input} '{print i1"\t"i2"\t"$2"\t"$3"\t"$4"\t"$5}' >> NCIS_aggregate.txt
done

echo "library(gridExtra)" > NCIS_aggregate.R
echo "library(ggplot2)" >> NCIS_aggregate.R
echo "library(scales)" >> NCIS_aggregate.R
echo "data <- read.delim('NCIS_aggregate.txt', head=TRUE)" >> NCIS_aggregate.R
echo "ChIP_cell <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 2)" >> NCIS_aggregate.R
echo "Control_cell <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 2)" >> NCIS_aggregate.R
echo "ChIP_mark <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 3)" >> NCIS_aggregate.R
echo "Control_mark <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 3)" >> NCIS_aggregate.R
echo "data <- cbind(data, ChIP_cell, Control_cell, ChIP_mark, Control_mark)" >> NCIS_aggregate.R
echo "data <- data[with(data, order(ChIP_cell, ChIP_mark, Control_cell, Control_mark)),]"  >> NCIS_aggregate.R
echo "data[,1] <- factor(data[,1], levels=unique(data[,1]))"  >> NCIS_aggregate.R
echo "data[,2] <- factor(data[,2], levels=unique(data[,2]))"  >> NCIS_aggregate.R
echo "fmt <- '%.2f'" >> NCIS_aggregate.R
echo "p1 <- ggplot(data, aes(Control,ChIP)) + geom_tile(aes(fill = BackgroundRatio), colour = 'white') + scale_fill_gradient2(limits=c(0.25, 1.),midpoint=0.63, high=muted('red'), mid='yellow', low=muted('steelblue'), na.value='steelblue') + theme(axis.text.x = element_text(angle = 45, hjust=1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, BackgroundRatio)),size=2)"  >> NCIS_aggregate.R
echo "p2 <- ggplot(data, aes(Control,ChIP)) + geom_tile(aes(fill = Normalization), colour = 'white') + scale_fill_gradient2(low = 'steelblue', mid='yellow', high = 'red') + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, Normalization)),size=2)" >> NCIS_aggregate.R
echo "p3 <- ggplot(data, aes(x=Control,y=ChIP,z=Binsize), log10='z') + geom_tile(aes(fill = Binsize), colour = 'white') + scale_fill_gradient2(low = 'steelblue', mid='yellow', high = 'red') + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')" >> NCIS_aggregate.R
echo "pdf(file='NCIS_aggregate.pdf')" >> NCIS_aggregate.R
echo "p1" >> NCIS_aggregate.R
echo "p2" >> NCIS_aggregate.R
echo "p3" >> NCIS_aggregate.R
echo "dev.off()" >> NCIS_aggregate.R

/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/Rscript NCIS_aggregate.R
