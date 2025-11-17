## ========== using expresion of mtDNA encoded gene (MT-CO2) to predict hallmark gene expression ==== 

library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)


data = read.csv('../../geneExpression_data/gene/Processed3_pcawg.rnaseq.transcript.expr.counts.csv',header=T)
donid = read.csv('../../geneExpression_data/gene/DonorID_canctype_mtDNA-CN.csv',header=T)

dim(data)
rownames(data) = data$Gene 
data = data[-1]
 
length(donid$mtDNA.CN)


hmfile = "../../../../datasets/hallmarks_msigdb/h.all.v2025.1.Hs.symbols.gmt"
sigline <- readLines(hmfile)


sink('results_mtGeneExpr_pred_HallmarkExpr/Hallmark_wise_cor_pval.txt') 
cat(paste('Hallmark','Mean_cor','Sd_cor','Num_sig_pval','Tot_genes','Perc','\n',sep=' '))
sink()

CON <- file("results_mtGeneExpr_pred_HallmarkExpr/Hallmark_wise_cor_pval.txt", "a")


### sigline = c('MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB')

sig_data = data[which(rownames(data) == 'MT-CO2'),] 

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
	#print(sig)

	subdir = elems[[1]][1] 
	dir.create(file.path('results_mtGeneExpr_pred_HallmarkExpr',subdir ), showWarnings = FALSE)
	sink(paste('results_mtGeneExpr_pred_HallmarkExpr/',subdir,'/','0_Prediction_correlation_results.txt',sep=''))

	filt_data = data[which(rownames(data) %in% unique(genelist)),]  

	corlist = vector()
	pvlist = 0

	for (k in 1:dim(filt_data)[1]) { 
		####  
		inputData = data.frame(gene = (t(filt_data[k,])),MT_CO2 = t(sig_data)) 
		colnames(inputData) = c('gene','MT_CO2')
		name = colnames(t(filt_data[k,]))	
	
		trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData)) # indices for 80% training data
		trainingData <- inputData[trainingIndex,]  # training data
		testData <- inputData[-trainingIndex,]  # test data
	
		lmmod <- lm(gene ~ MT_CO2, data=trainingData)
		summary(lmmod) ## Performance on training data 

		preds <- predict(lmmod, testData, na.action=na.pass)  # predict on test data

		actual <- testData$gene
		full <- as.data.frame(cbind(actual,preds))


		wrf = paste('results_mtGeneExpr_pred_HallmarkExpr/',subdir,'/',name,'_expr_pred_MT-CO2.pdf',sep='')
          
		pdf(file=wrf,useDingbats=F)
		ggp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0.5,0.8,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
		theme_minimal() + labs(x="Actual expression",y="Predicted expression",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') 
		print(ggp)
		dev.off() 


		actual=full$actual[which(full$preds!="NA")]
		preds=full$preds[which(full$preds!="NA")]

		pv=cor.test(actual,preds,method='spearman')
		cat(paste(name,pv$estimate,pv$p.value,'\n',sep='\t'))
		## print(cor.test(actual,preds))

		corlist = append(corlist,pv$estimate) 
		if(pv$p.value<0.05) {
			pvlist = pvlist + 1 
		}
 
		rm(wrf) 
		rm(inputData)
		rm(lmmod)
		rm(name)
		rm(actual)
		rm(preds)
		rm(ggp) 
		rm(pv)
		rm(full)
		rm(trainingData)
		rm(testData)		
	}
	sink()
	avgcor = round(mean(corlist),2)
	sdcor = round(sd(corlist),2)

	perc = round((pvlist/dim(filt_data)[1]*100),2)
	towrite = paste(elems[[1]][1],' ',avgcor,sdcor,pvlist,dim(filt_data)[1],perc,'\n',sep=' ')
	print(towrite)
	### cat(towrite, file = "results_mtGeneExpr_pred_HallmarkExpr/Hallmark_wise_cor_pval.txt", append = TRUE)
	cat(towrite, file = CON)

	rm(perc)
	rm(towrite)
	rm(sdcor)
	rm(avgcor)
	rm(elems)
	rm(filt_data)
	rm(subdir)
	rm(corlist)
	rm(pvlist) 
}


close(CON) 

rm(hmfile)

rm(data)
rm(sig)
rm(sig_data)
rm(sigline)
rm(donid) 
rm(CON) 