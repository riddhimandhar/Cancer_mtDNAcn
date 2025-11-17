## Gather statistics of genes from Pan_cancer signature lists 

### =======  EMTome signatures ==============


library(gplots)
library(viridis) 

readfile = "../results_GeneExprAnalysis/PanCancer_PosNeg_genes/PanCancer_genes_logFC.txt"
data = read.table(readfile,header=T)

num_rows = dim(data)[2]

ind_files <- list.files(path = "../../../Relevant_papers/Cancer_drivers_CSC_signatures/EMTome_signature/ind_Sigs/", pattern = "\\.txt$", full.names = F)

full_matrix = matrix(NA, nrow = num_rows, ncol = 0)

for (file in ind_files) {
	#file = 'MsigDB_v7.0.txt'
	rdfile = paste('../../../Relevant_papers/Cancer_drivers_CSC_signatures/EMTome_signature/ind_Sigs/',file,sep='')
	tmp = read.table(rdfile,header=T,check.names=F) 
	cancers = strsplit(colnames(tmp), ",")
	print(cancers[[1]])
	num_canc = length(cancers[[1]])

		
	if (cancers[[1]][1] == 'Pan_cancer') { 
		
		newcol = list()
		filt_data = data[which(rownames(data) %in% unique(tmp$Pan_cancer)),]  
		for (name in colnames(filt_data)) {
			print(name) 
			newname = gsub('\\.','-',name)
			indcancfile = paste('../results_GeneExprAnalysis/DiffExprAnalysis_gene/DESeq_results/',newname,'_High_vs_Low_mtDNA.txt',sep='') 
			ind_data = read.table(indcancfile,header=T)
			rm(newname) 

			posfrac = 0
		        negfrac = 0 
			vals = as.numeric(filt_data[,name])
			pvs = ind_data[rownames(filt_data),]$padj

		 	for (k in 1:length(vals)) {
				if (vals[k]>0 & is.na(pvs[k]) == 0 & pvs[k] <= 0.1) { posfrac = posfrac+1 } 
				#if (vals[k] >= 0.2) { posfrac = posfrac+1 } 
				if (vals[k]<0 & is.na(pvs[k]) == 0 & pvs[k] <= 0.1) { negfrac = negfrac+1 }
				#if (vals[k] <= -0.2) { negfrac = negfrac+1 }  
			}
			posfrac = round(posfrac/length(vals),3) 
			negfrac = round(negfrac/length(vals),3) 	
			negenrich = round((0.01+negfrac)/(0.01+posfrac),3) 
			newcol = append(newcol,negenrich)
			#print(name)
			rm(vals) 
			rm(pvs)
		}
		rm(name)
		rm(filt_data) 	

		full_matrix = cbind(full_matrix, file = newcol)
		colnames(full_matrix)[(colnames(full_matrix) == 'file')] = file 
		rm(newcol)

	} 

	rm(tmp) 
	rm(rdfile)
	rm(cancers)
}

rm(file)
rm(ind_files) 

rownames(full_matrix) = colnames(data)
print(full_matrix)

write.table(t(full_matrix),'../results_GeneExprAnalysis/Heatmaps_EMT_signatures/EMTome_LomtDNA-CN_negEnrich_PanCancer_Sigs.txt',quote=F)

rm(full_matrix)
rm(data) 	


tab = read.table('../results_GeneExprAnalysis/Heatmaps_EMT_signatures/EMTome_LomtDNA-CN_negEnrich_PanCancer_Sigs.txt')
boxplot(tab,las=2,ylim=c(0,10),col='dodgerblue3')
abline(h=1,lwd=2,lty=2,col='red') 
dev.copy2pdf(file='../results_GeneExprAnalysis/Heatmaps_EMT_signatures/EMTome_LomtDNA-CN_negEnrich_PanCancer_padj0.1.pdf',width=9,height=6)
dev.off()


