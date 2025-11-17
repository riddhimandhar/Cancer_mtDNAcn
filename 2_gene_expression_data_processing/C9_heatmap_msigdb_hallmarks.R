## dir: "C:/Users/Riddhiman Dhar/Dropbox/IITKGP_stuff/Manuscripts/Mitochondria_CN_cancer/programs_GeneExprAnalysis

library(gplots)
library(viridis) 

my_palette <- colorRampPalette(c("dodgerblue1", "white", "darkorange1"))(100)
# my_palette <- viridis(option='E', 50, direction = -1)


readfile = "../results_GeneExprAnalysis/PanCancer_PosNeg_genes/PanCancer_genes_logFC.txt"
data = read.table(readfile,header=T)

hmfile = "../../../datasets/hallmarks_msigdb/h.all.v2025.1.Hs.symbols.gmt"
sigline <- readLines(hmfile)

for (sig in sigline) { 
	elems = strsplit(sig,"\t") 
	num_elems = length(elems[[1]]) 
	#print(noquote(paste(elems[[1]][1],num_elems,sep=' ')))
	genelist = vector() 
	for (i in 3:num_elems) {
		genelist = append(genelist,noquote(elems[[1]][i]))  
	}

	print(elems[[1]][1])
	#print(genelist) 

	## heatmap across all cancer types 
	filt_data = data[which(rownames(data) %in% unique(genelist)),]  
			
	wrfile = paste('../results_GeneExprAnalysis/Heatmaps_EMT_signatures/hallmarks_msigdb/Pan_cancer',elems[[1]][1],'.pdf',sep='')
	pdf(file = wrfile, useDingbats = F, width = 16, height = 12) 
	heatmap.2(as.matrix(filt_data), scale = "none", dendrogram = "none", Rowv = T, Colv = F, trace = "none", 
	col=my_palette, breaks = seq(-1,1,length.out=length(my_palette)+1), keysize=1, key.xlab='log2FC',margins = c(12, 5))
	dev.off() 

	rm(num_elems)
	rm(elems) 
	rm(filt_data)
}


rm(sig)
rm(sigline)


ind_files <- list.files(path = "../../../Relevant_papers/Cancer_drivers_CSC_signatures/EMTome_signature/ind_Sigs/", pattern = "\\.txt$", full.names = F)

for (file in ind_files) {
	rdfile = paste('../../../Relevant_papers/Cancer_drivers_CSC_signatures/EMTome_signature/ind_Sigs/',file,sep='')
	tmp = read.table(rdfile,header=T,check.names=F) 
	cancers = strsplit(colnames(tmp), ",")
	print(cancers[[1]])
	num_canc = length(cancers[[1]])

	#for (canctype in cancers) {
		part = strsplit(file,'.txt') 
		
		if (cancers[[1]][1] == 'Pan_cancer') { 
			## heatmap across all cancer types 
			filt_data = data[which(rownames(data) %in% unique(tmp$Pan_cancer)),]  
			
			wrfile = paste('../results_GeneExprAnalysis/Heatmaps_EMT_signatures/EMTome_ind_Sigs/Pan_cancer',part,'.pdf',sep='')
			pdf(file = wrfile, useDingbats = F, width = 16, height = 12) 
			heatmap.2(as.matrix(filt_data), scale = "none", dendrogram = "none", Rowv = T, Colv = F, trace = "none", 
				  col=my_palette, breaks = seq(-1,1,length.out=length(my_palette)+1), keysize=1, key.xlab='log2FC',margins = c(12, 5))
			dev.off() 
		} 
		if (cancers[[1]][1] != 'Pan_cancer') { 
			## heatmap for the specific cancer types 

			filt_data = as.data.frame(data[which(rownames(data) %in% unique(tmp[,1])),])
			rn = rownames(filt_data) 
			filt_data = as.data.frame(filt_data[,which(colnames(filt_data) %in% cancers[[1]])])
			rownames(filt_data) = rn 		
	
			if (length(filt_data) != 0)  {
				
				filt_data = as.matrix(cbind(rep(0,dim(filt_data)[1]),filt_data)) 
				colnames(filt_data) = c('Dummy',cancers[[1]])

				wrfile = paste('../results_GeneExprAnalysis/Heatmaps_EMT_signatures/EMTome_ind_Sigs/',part,'.pdf',sep='')	
				pdf(file = wrfile, useDingbats=F, width = 6, height = 12) 
				heatmap.2(as.matrix(filt_data), scale = "none", dendrogram = "none", Rowv = T, Colv = F, trace = "none", 
					  col=my_palette, breaks = seq(-1,1, length.out=length(my_palette)+1), cexCol=1.2, keysize=1, key.xlab='log2FC',margins = c(12, 5))
				dev.off()
			} 
		}
		
		rm(filt_data) 
		rm(wrfile)
		rm(part)  
	#}
	rm(tmp) 
	rm(rdfile)
	rm(cancers)
}

rm(file)
rm(ind_files) 


rm(data) 	
