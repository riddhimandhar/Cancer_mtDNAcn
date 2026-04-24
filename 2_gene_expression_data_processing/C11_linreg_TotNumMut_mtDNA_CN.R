library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)


data = read.table('../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt',header=T)

filt_data = data[,c(9,12,15,16,20)] 


inputData <- filt_data[,which(apply(filt_data, 2, var, na.rm=TRUE) != 0)]


encoded_gender <- model.matrix(~donor_sex-1, data=data)
inputData = cbind(inputData,encoded_gender)


lmData = inputData[which((is.na(inputData$MT_CopyNumber) == 'FALSE') & (is.na(inputData$Total_Num_Mutations) == 'FALSE') & (is.na(inputData$donor_age_at_diagnosis_years) == 'FALSE')),]

lmData = lmData[,which(colSums(abs(na.omit(lmData)))!=0)]
lmData = lmData[,which(colSums(is.na(lmData))<0.9*nrow(lmData))]
lmData = lmData[which(rowSums(is.na(lmData))<0.2*ncol(lmData)),]

lmData = lmData[,which(colSums(na.omit(lmData)==0)<0.8*nrow(lmData))]


 full_linmod <- lm(Total_Num_Mutations ~ ., data = lmData)

 mod_sum=summary(full_linmod)
 coefs = mod_sum$coefficients[,1]
 pval = mod_sum$coefficients[,4]
 key_variables=names(which(mod_sum$coefficients[,4]<0.05))
 print(key_variables) 

 selected = as.data.frame(coefs[names(coefs) %in% key_variables])
 selected = cbind(key_variables,selected)
 colnames(selected) = c('Feature', 'coefficient')
 selected = selected[order(abs(selected$coefficient),decreasing=T),]

 write.table (selected,"../results_June2025/lin_reg_models_rand_forest/Linmod_key_variables_TotNumMut.txt",quote=F,row.names=F)
 rm(selected) 

 cor_p = vector()
 cor_s = vector()
 pval_p = vector()
 pval_s = vector() 
 fl = 0 

 for (ct in 1:100) { 

  trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData)) # indices for 80% training data
  trainingData <- inputData[trainingIndex,]  # training data
  testData <- inputData[-trainingIndex,]  # test data

  lmmod <- lm(Total_Num_Mutations ~ . , data=trainingData)
  
  ## summary(lmmod) ## Performance on training data 

  expv=summary(lmmod)[9][[1]]
  preds <- predict(lmmod, testData, na.action=na.pass)  # predict on test data

  actual <- testData$Total_Num_Mutations
  full <- as.data.frame(cbind(actual,preds))

  actual=full$actual[which(full$preds!="NA")]
  preds=full$preds[which(full$preds!="NA")]

  rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))
  tss <- sum(na.omit(actual - mean(na.omit(actual))) ^ 2/length(na.omit(actual-mean(na.omit(actual)))))
  rsq <- 1 - rss/tss

  #rsq  ## Rsq on test data 
  cs = cor.test(actual,preds,method='spearman')
  cor_s = append(cor_s, cs$estimate)
  pval_s = append(pval_s, cs$p.value) 
 
  cp = cor.test(actual,preds)
  cor_p = append(cor_p, cp$estimate)
  pval_p = append(pval_p, cp$p.value) 

  rm(cs)
  rm(cp) 

  
  if(fl==0) {
     tt = paste(cor_s,pval_s)
     fl = 1
     pdf(file="../results_June2025/lin_reg_models_rand_forest/Linmod_actual_predicted_mutational_load.pdf",useDingbats=F) 
     gp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
          theme_minimal() + labs(x="Actual mutational load",y="Predicted mutational load",title= tt) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  
     print(gp)
     dev.off() 
     rm(gp) 
  }
}

print(paste(mean(cor_s),sd(cor_s)))
print(paste(mean(pval_s),sd(pval_s)))
print(paste(mean(cor_p),sd(cor_p)))
print(paste(mean(pval_p),sd(pval_p)))

rm(cor_s)
rm(pval_s)
rm(cor_p)
rm(pval_p)
