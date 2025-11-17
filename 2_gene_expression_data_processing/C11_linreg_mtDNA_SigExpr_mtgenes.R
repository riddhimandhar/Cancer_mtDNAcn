## ======== using mtDNA to predict mtDNA encoded gene expression ==== 

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


genelist = c('MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB')

filt_data = data[which(rownames(data) %in% unique(genelist)),]  


sink('results_mtDNApred_mtDNAgeneExpr/Prediction_correlation_results.txt')

	for (k in 1:dim(filt_data)[1]) { 
		inputData = data.frame(gene = (t(filt_data[k,])),mtCN = donid$mtDNA.CN) 
		colnames(inputData) = c('gene','mtCN')
		name = colnames(t(filt_data[k,]))
	
		trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData)) # indices for 80% training data
		trainingData <- inputData[trainingIndex,]  # training data
		testData <- inputData[-trainingIndex,]  # test data
	
		lmmod <- lm(gene ~ mtCN, data=trainingData)
		summary(lmmod) ## Performance on training data 

		preds <- predict(lmmod, testData, na.action=na.pass)  # predict on test data

		actual <- testData$gene
		full <- as.data.frame(cbind(actual,preds))

		wrf = paste('results_mtDNApred_mtDNAgeneExpr/',name,'_expr_pred_mtCN.pdf',sep='')
          
		pdf(file=wrf,useDingbats=F)
		ggp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0.5,0.8,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
		theme_minimal() + labs(x="Actual expression",y="Predicted expression",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') 
		print(ggp)
		dev.off() 


		actual=full$actual[which(full$preds!="NA")]
		preds=full$preds[which(full$preds!="NA")]

		pv = cor.test(actual,preds,method='spearman')
		cat(paste(name,pv$estimate,pv$p.value,'\n',sep='\t'))

		# print(cor.test(actual,preds))

		rm(name)
		rm(actual)
		rm(preds)
		rm(pv)
		rm(full)
		rm(ggp)
		rm(wrf)
		rm(lmmod) 
		rm(testData)
		rm(trainingData)
		rm(trainingIndex) 
	}

sink()


rm(donid)
rm(filt_data)
rm(data)
rm(genelist) 



