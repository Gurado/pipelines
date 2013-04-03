#!/bin/sh -e

echo "ChIP	Control	Normalization	Enriched	Cumulative	FDR" > CHANCE_aggregate.txt

for f in `ls results/*.IPstrength`; do

	b=$(basename $f)
        t=(${b//./ })
	chip=(${t%-*})
	input=(${t##*-})

	FDR=$(cat $f | grep "^fdr," | awk -F\, '{print $2}')
	NORM=$(cat $f | grep "input_scaling_factor" | awk -F\, '{print $2}')
	ENRICHED=$(cat $f | grep "percent_genome_enriched" | awk -F\, '{print $2}')
	CUMMULATIVE=$(cat $f | grep "differential_percentage_enrichment" | awk -F\, '{print $2}')

        echo "${chip}	${input}	${NORM}	${ENRICHED}	${CUMMULATIVE}	${FDR}" >> CHANCE_aggregate.txt
done

echo "library(gridExtra)" > CHANCE_aggregate.R
echo "library(ggplot2)" >> CHANCE_aggregate.R
echo "library(scales)" >> CHANCE_aggregate.R
echo "data <- read.delim('CHANCE_aggregate.txt', head=TRUE)" >> CHANCE_aggregate.R
echo "data <- cbind(data, 'BackgroundRatio'=0)" >> CHANCE_aggregate.R
echo "data['BackgroundRatio'] <- (100 - data['Cumulative'])/100"  >> CHANCE_aggregate.R
echo "ChIP_cell <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 2)" >> CHANCE_aggregate.R
echo "Control_cell <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 2)" >> CHANCE_aggregate.R
echo "ChIP_mark <- sapply(strsplit(as.character(data[,1]), '_'), '[[', 3)" >> CHANCE_aggregate.R
echo "Control_mark <- sapply(strsplit(as.character(data[,2]), '_'), '[[', 3)" >> CHANCE_aggregate.R
echo "data <- cbind(data, ChIP_cell, Control_cell, ChIP_mark, Control_mark)" >> CHANCE_aggregate.R
echo "data <- data[with(data, order(ChIP_cell, ChIP_mark, Control_cell, Control_mark)),]"  >> CHANCE_aggregate.R
echo "data[,1] <- factor(data[,1], levels=unique(data[,1]))"  >> CHANCE_aggregate.R
echo "data[,2] <- factor(data[,2], levels=unique(data[,2]))"  >> CHANCE_aggregate.R
echo "data['FDR'] <- sapply(data['FDR'], log10)" >> CHANCE_aggregate.R
echo "data['FDR'] <- -data['FDR'] "  >> CHANCE_aggregate.R
echo "fmt <- '%.2f'" >> CHANCE_aggregate.R
echo "p1 <- ggplot(data, aes(Control,ChIP)) + geom_tile(aes(fill = BackgroundRatio), colour = 'white') + scale_fill_gradient2(limits=c(0.25, 1.),midpoint=0.63, high=muted('red'), mid='yellow', low=muted('steelblue'), na.value='steelblue') + theme(axis.text.x = element_text(angle = 45, hjust=1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, BackgroundRatio)),size=2)"  >> CHANCE_aggregate.R
echo "p2 <- ggplot(data, aes(Control,ChIP)) + geom_tile(aes(fill = Normalization), colour = 'white') + scale_fill_gradient2(name='input scaling factor', low = 'steelblue', mid='yellow', high = 'red') + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, Normalization)),size=2)" >> CHANCE_aggregate.R
echo "p3 <- ggplot(data, aes(x=Control,y=ChIP)) + geom_tile(aes(fill = Enriched), colour = 'white') + scale_fill_gradient2(name='percent genome enriched', midpoint=20, low = 'steelblue', mid='yellow', high = 'red') + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, Enriched)),size=2)" >> CHANCE_aggregate.R
echo "p4 <- ggplot(data, aes(Control,ChIP,FDR)) + geom_tile(aes(fill = FDR), colour = 'white') + scale_fill_gradient2(name='-log10(FDR)', limits=c(0,3), midpoint=1.5, high=muted('steelblue'), mid='yellow', low=muted('red'), na.value='steelblue') + theme(axis.text.x = element_text(angle = 45, hjust=1, size=8), axis.text.y=element_text(size=8)) + theme(legend.direction = 'horizontal', legend.position = 'top')+geom_text(aes(label=sprintf(fmt, FDR)),size=2)"  >> CHANCE_aggregate.R
echo "pdf(file='CHANCE_aggregate.pdf')" >> CHANCE_aggregate.R
echo "p1" >> CHANCE_aggregate.R
echo "p2" >> CHANCE_aggregate.R
echo "p3" >> CHANCE_aggregate.R
echo "p4" >> CHANCE_aggregate.R
echo "dev.off()" >> CHANCE_aggregate.R

/share/ClusterShare/software/contrib/Cancer-Epigenetics/tools/bin/Rscript CHANCE_aggregate.R
