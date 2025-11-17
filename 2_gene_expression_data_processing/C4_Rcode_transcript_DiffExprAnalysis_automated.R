library(DESeq2)
library(pheatmap)
library(RColorBrewer)	
library(PCAtools) 
library(genefilter)
library(viridis)
library(EnhancedVolcano)


#setwd('../Dropbox/IITKGP_stuff/Manuscripts/Mitochondria_CN_cancer/Results/myAnalysis_may2025/programs_GeneExprAnalysis') 
#setwd('../Dropbox/IITKGP_stuff/Manuscripts/Mitochondria_CN_cancer/Results/myAnalysis_may2025/programs_GeneExprAnalysis') 

sink("../results_GeneExprAnalysis/DiffExprAnalysis_transcript/Stats_Up_Down_reg.txt") 

canctypelist = c('Biliary-AdenoCA' , 
                 'Bladder-TCC',
	         'Breast-AdenoCA',
                 'CNS-GBM',
                 'CNS-Oligo',
                 'Cervix-SCC',
                 'ColoRect-AdenoCA',
                 'Head-SCC',
                 'Kidney-ChRCC',
                 'Kidney-RCC',
                 'Liver-HCC',
                 'Lung-AdenoCA',
                 'Lung-SCC',
                 'Lymph-BNHL',
                 'Lymph-CLL',
                 'Ovary-AdenoCA',
                 'Panc-AdenoCA',
                 'Prost-AdenoCA',
                 'SoftTissue-Leiomyo',
                 'SoftTissue-Liposarc',
                 'Stomach-AdenoCA',
                 'Thy-AdenoCA',
                 'Uterus-AdenoCA')

for (canctype in canctypelist) { 
	#print(canctype) 

        datafile = paste('../geneExpression_data/transcript/canctype_specific_data/',canctype,'_counts_data.txt',sep="") 
	sampfile = paste('../geneExpression_data/transcript/canctype_specific_data/Metadata/Metadata_',canctype,'.txt',sep="") 

	data=read.table(datafile,header=T)  ## Data
	samp=read.table(sampfile,header=T)  ## Metadata
        rm(datafile)
	rm(sampfile) 	

	datamatrix=as.data.frame(data[-1])
	rownames(datamatrix) = data[,1]
	#rownames(datamatrix) = make.unique(rownames(datamatrix)) ##
	sampmatrix=as.data.frame(samp[,c(1,3)])
	colnames(sampmatrix) = c('Samples','mtDNA_CN_class') 
	sampmatrix$mtDNA_CN_class = as.factor(sampmatrix$mtDNA_CN_class)
	rownames(sampmatrix) = gsub('-','.',sampmatrix$Samples)

	dds <- DESeqDataSetFromMatrix(countData = datamatrix, colData=sampmatrix, design=~mtDNA_CN_class)
	#dds

	vstdata <- vst(dds, blind = F) ##, nsub = 1000, fitType = "parametric")
	#head(assay(vstdata),3)

	
	#par(mfrow=c(1, 2))
	#mks <- estimateSizeFactors(dds)
	#lims <- c(-2, 20)
	#plot(log2(counts(mks, normalized=TRUE)[,1:2] + 1),pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)
	#plot(assay(vstdata)[,1:2], pch=16, cex=0.3, main="VST", xlim=lims, ylim=lims)
        #rm(lims) 
        #rm(mks)


	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/Distanceplots/',canctype,'_distmatrix.pdf',sep="") 

	sampleDists <- dist(t(assay(vstdata)))
	#sampleDists

	sampleDistMatrix <- as.matrix( sampleDists )
	rownames(sampleDistMatrix) <- vstdata$mtDNA_CN_class  #paste(vstdata$mtDNA_CN_class, vstdata$Samples, sep="-" )
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, "YlOrBr")) )(255)
	
	pdf(file = wrf, useDingbats=F) 
	pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
	dev.off() 

	rm(sampleDists)
        rm(sampleDistMatrix) 
        rm(colors) 
        rm(data)
        rm(samp) 
	rm(wrf)
		

## PCAtools 

	p <- pca(assay(vstdata), metadata = sampmatrix, removeVar = 0.1)
	#screeplot(p, axisLabSize = 12, titleLabSize = 16)

	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/PCAplots/',canctype,'_PCA.pdf',sep="")
	
	pdf(file = wrf, useDingbats=F) 
	bip = biplot(p, showLoadings = F, lab = NULL, pointSize = 1, sizeLoadingsNames = 5, colby='mtDNA_CN_class', legendPosition='right',ntopLoadings=15)
	print(bip) 
	dev.off() 
    
        rm(bip) 
	rm(wrf) 

	#plotloadings(p, rangeRetain = 0.01, labSize = 4.0, title = 'Loadings plot', 
    	#	subtitle = 'PC1, PC2, PC3, PC4, PC5', 
    	#	caption = 'Top 1% variables',
    	#	col = c('dodgerblue', 'grey', 'orange'),
    	#	drawConnectors = TRUE)

	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/PCAplots/',canctype,'_pairplots_PCA.pdf')
	pdf(file = wrf, width=8, height=6, useDingbats=F) 

	pplot = pairsplot(p, components = getComponents(p, c(1:5)),
    		triangle = TRUE, trianglelabSize = 6,
    		hline = 0, vline = 0,
    		pointSize = 0.2,
    		gridlines.major = FALSE, gridlines.minor = FALSE,colby='mtDNA_CN_class') 
	
	print(pplot) 
	dev.off()

	rm(pplot) 
	rm(wrf) 

	## biplot(p, x = 'PC4', y = 'PC3', colby = 'mtDNA_CN_class', hline = 0, vline = 0, legendPosition = 'right',showLoadings = F, lab=NULL, pointSize = 1, ntopLoadings=20) 
 	
	rm(p)



## ==== Differential expression ==== 

	prdds <- DESeq(dds)
	res <- results(prdds,contrast = c("mtDNA_CN_class", "HighCopy", "LowCopy"))
	#res 

	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/MAplots/',canctype,'_MAplot.pdf',sep="")
	pdf(file = wrf, useDingbats=F) 

	plotMA(res, ylim=c(-5,5))

	dev.off() 
	rm(wrf) 

	#mcols(res, use.names=TRUE)

	wfile = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/DESeq_results/',canctype,'_High_vs_Low_mtDNA.txt',sep='')
	#write.table(res,'../results_GeneExprAnalysis/DiffExprAnalysis_transcript/DESeq_results/Liver-HCC_Low_vs_High_mtDNA.txt',quote=F)
	write.table(res,wfile,quote=F)
	
	rm(wfile) 


	res_padjfilt = res[which((res$padj < 0.1) == 'TRUE'),]
	upreg = res[which((res$padj < 0.1) == 'TRUE' & (res$log2FoldChange > 0) == 'TRUE'),]
	dnreg = res[which((res$padj < 0.1) == 'TRUE' & (res$log2FoldChange < 0) == 'TRUE'),]

	wfile1 = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/DESeq_results/',canctype,'_Upreg_Low_vs_High_mtDNA.txt',sep='')
	write.table(upreg,wfile1,quote=F)
	rm(wfile1) 
	
	wfile2 = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/DESeq_results/',canctype,'_Dnreg_Low_vs_High_mtDNA.txt',sep='')
	write.table(dnreg,wfile2,quote=F)
	rm(wfile2) 

	prtext = paste(canctype,' Upreg:',nrow(upreg),' Dnreg:',nrow(dnreg))
	cat(c(prtext,"\n"))

	rm(prtext) 
	rm(res_padjfilt)
	rm(upreg)
	rm(dnreg) 


## ===== Volcano plot ======== 

	topDiffGenes <- head(order(res$padj,decreasing=F),5)
	topDiffGeneNames <- rownames(datamatrix[topDiffGenes,])

	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/Volcanoplots/',canctype,'_Volcanoplot.pdf',sep="")
	
	pdf(file = wrf, width=8, height=6, useDingbats=F) 

	evplot = EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj',
    		xlab = bquote(~Log[2]~ 'fold change'),
    		ylab = bquote(-~Log[10]~padj),
    		selectLab = topDiffGeneNames, 
    		pCutoff = 0.1,
    		FCcutoff = 1.0,
    		pointSize = 1.0,
    		labSize = 4.0,
    		labCol = 'black',
    		labFace = 'bold',
    		boxedLabels = TRUE,
    		colAlpha = 4/5,
    		legendLabels=c('Not sig.','Log (base 2) FC','p-value', 'p-value & Log (base 2) FC'),
    		legendPosition = 'right',
    		legendLabSize = 12,
    		legendIconSize = 4.0,
    		drawConnectors = T,
    		widthConnectors = 1.0,
    		colConnectors = 'black')

	print(evplot)
	dev.off() 
	rm(wrf) 

	rm(evplot) 
	rm(topDiffGenes) 
	rm(topDiffGeneNames) 


## ====== Topgene Normalized count plot ========= 

	topGene <- rownames(res)[which.min(res$padj)]
	## topGene <- rownames(res)[which.max(res$log2FoldChange)]

	geneCounts <- plotCounts(prdds, gene=topGene, intgroup= "mtDNA_CN_class", returnData=TRUE)

	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/TopGeneExprplot/',canctype,'_TopgeneExpr.pdf',sep="")

	pdf(file = wrf, useDingbats=F)
	ggp = ggplot(geneCounts, aes(x=mtDNA_CN_class, y=count, color= mtDNA_CN_class)) + scale_y_log10() + geom_point(position=position_jitter(width=0.2,height=0), size=3) + labs(title = topGene)
	print(ggp) 
	dev.off()

	rm(wrf) 	

	rm(ggp) 
	rm(geneCounts) 
	rm(topGene) 


## ========= Heatmap Top N gene expression ===========

	resfilt = res[which(res$padj < 0.1),]

	topVarGenes <- head(order(rowVars(assay(vstdata)),decreasing=TRUE),20)
	#topDiffGenes <- head(order(abs(resfilt$log2FoldChange),decreasing=F),40)
	topDiffGenes <- head(order(resfilt$padj,decreasing=F),40)

	mat <- assay(vstdata)[ topDiffGenes,]
	mat <- mat - rowMeans(mat)
	mat <- mat[,order(colnames(mat),decreasing=TRUE)]

	val <- colData(vstdata)
	val <- val[order(rownames(val),decreasing=T),]

	df <- as.data.frame(val)[,c("mtDNA_CN_class","Samples")]

	my_colors <- colorRampPalette(inferno(n = 256))(256)
	
	wrf = paste('../results_GeneExprAnalysis/DiffExprAnalysis_transcript/Top40heatmap/',canctype,'_Top40Heatmap.pdf',sep="")

	pdf(file = wrf, useDingbats=F)
	#pheatmap(mat, show_colnames = F, color=my_colors, cluster_cols=F, annotation_col=annotation_col)
	pheatmap(mat, show_colnames = F, color=my_colors, cluster_cols=F, annotation_col= df, annotation_legend = F)
	dev.off() 

	rm(prdds)
	rm(wrf) 
	rm(mat)
	rm(val)
	rm(df)
	rm(resfilt)
	rm(res)
	rm(vstdata) 
        rm(datamatrix)
        rm(sampmatrix) 
        rm(topVarGenes)
        rm(topDiffGenes) 
	rm(dds) 

}

sink()
