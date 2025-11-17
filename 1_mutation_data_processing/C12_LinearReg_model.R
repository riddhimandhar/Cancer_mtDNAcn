library(pls)
library(dplyr)
library(tidyr)


data = read.table('../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt',header=T)

filt_data = data[,c(9,12:18,20:65,71:77,80:114)] 


inputData <- filt_data[,which(apply(filt_data, 2, var, na.rm=TRUE) != 0)]

lmData = inputData[which(is.na(inputData$MT_CopyNumber) == 'FALSE'),] 
lmData = lmData[,which(colSums(abs(na.omit(inputData)))!=0)]
lmData = lmData[,which(colSums(is.na(lmData))<0.9*nrow(lmData))]
lmData = lmData[,which(colSums(lmData==0)<0.9*nrow(lmData))]
lmData = lmData[which(rowSums(is.na(lmData))<0.2*ncol(lmData)),]




### =========== Predictive modeling =======================

## Ridge and Lasso regression 


library(ridge)
library(glmnet)


linRidgeMod <- linearRidge(MT_CopyNumber ~ ., data = lmData)
mod_sum=summary(linRidgeMod)

names(mod_sum)

npc=mod_sum$chosen.nPCs

head(mod_sum$summaries$summary1$coefficients)

coefs = mod_sum$summaries$summary1$coefficients[,1]
pval=mod_sum$summaries$summary1$coefficients[,5]

key_variables=names(which(mod_sum$summaries$summary1$coefficients[,5]<0.05))
print(key_variables) 

selected = as.data.frame(coefs[names(coefs) %in% key_variables])
selected = cbind(key_variables,selected)
colnames(selected) = c('Feature', 'coefficient')
selected = selected[order(abs(selected$coefficient),decreasing=T),]

write.table(selected,"../results_June2025/lin_reg_models_rand_forest/Ridge_key_variables.txt",quote=F,row.names=F)
rm(selected) 


ridge_inputData=inputData[,which(colnames(inputData) %in% key_variables)]

ridge_inputData=cbind(inputData$MT_CopyNumber,ridge_inputData)

ridge_inputData=rename(ridge_inputData, 'MT_CopyNumber'='inputData$MT_CopyNumber')


  ridge_trainingIndex <- sample(1:nrow(ridge_inputData), 0.80*nrow(ridge_inputData)) # indices for 80% training data
  ridge_trainingData <- ridge_inputData[ridge_trainingIndex,]  # training data
  ridge_testData <- ridge_inputData[-ridge_trainingIndex,]  # test data

  lmmod_ridge <- lm(MT_CopyNumber ~ . , data=ridge_trainingData)

  ## summary(lmmod_ridge) ## Performance on training data 

  expv=summary(lmmod_ridge)[9][[1]]
  preds <- predict(lmmod_ridge, ridge_testData, na.action=na.pass)  # predict on test data

  actual <- ridge_testData$MT_CopyNumber
  full <- as.data.frame(cbind(actual,preds))

  actual=full$actual[which(full$preds!="NA")]
  preds=full$preds[which(full$preds!="NA")]

  rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))
  tss <- sum(na.omit(actual - mean(na.omit(actual))) ^ 2/length(na.omit(actual-mean(na.omit(actual)))))
  rsq <- 1 - rss/tss

  rsq  ## Rsq on test data 

  ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
  theme_minimal() + labs(x="Actual mtDNA CN",y="Predicted mtDNA CN",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

  dev.copy2pdf(file="../results_June2025/lin_reg_models_rand_forest/Ridge_actual_predicted_mtDNA-CN.pdf",useDingbats=F) 
  dev.off() 

  cor.test(actual,preds,method='spearman')
  cor.test(actual,preds)



### ==== Lasso ==========

len=length(names(lmData))
x = model.matrix(MT_CopyNumber~., lmData)[,-1]
y = lmData %>% select(MT_CopyNumber) %>% unlist() %>% as.numeric()
grid = 10^seq(10, -2, length = 100)

cv.out = cv.glmnet(x, y, alpha = 1)
bestlam = cv.out$lambda.min

out = glmnet(x, y, alpha = 1, lambda = grid)

lasso_coef = predict(out, type = "coefficients", s = bestlam)[1:len,] 


## New dataset with selected features 

impvar=lasso_coef[lasso_coef != 0]
print(impvar)
impvarlist=names(impvar[-1])
implen=length(names(impvar))

selected = as.data.frame(impvar[names(impvar) %in% impvarlist])
selected = cbind(impvarlist,selected)
colnames(selected) = c('Feature', 'coefficient')
selected = selected[order(abs(selected$coefficient),decreasing=T),]

write.table(selected,"../results_June2025/lin_reg_models_rand_forest/Lasso_key_variables.txt",quote=F,row.names=F)
rm(selected) 


lasso_inputData=inputData[,which(colnames(inputData) %in% impvarlist)]

lasso_inputData=cbind(inputData$MT_CopyNumber,lasso_inputData)

lasso_inputData=rename(lasso_inputData, 'MT_CopyNumber'='inputData$MT_CopyNumber')



## Training and testing 

  lasso_trainingIndex <- sample(1:nrow(lasso_inputData), 0.80*nrow(lasso_inputData)) # indices for 80% training data
  lasso_trainingData <- lasso_inputData[lasso_trainingIndex,] # training data
  lasso_testData <- lasso_inputData[-lasso_trainingIndex,] # test data

  lmmod_lasso <- lm(MT_CopyNumber ~ . , data=lasso_trainingData)
  expv=summary(lmmod_lasso)[9][[1]]
  expv
  
  preds <- predict(lmmod_lasso, lasso_testData, na.action=na.pass)  # predict on test data
  actual <- lasso_testData$MT_CopyNumber
  full <- as.data.frame(cbind(actual,preds))

  actual=full$actual[which(full$preds!="NA")]
  preds=full$preds[which(full$preds!="NA")]
  rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))
  tss <- sum(na.omit(actual - mean(na.omit(actual))) ^ 2/length(na.omit(actual-mean(na.omit(actual)))))
  rsq <- 1 - rss/tss

  rsq


  ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
  theme_minimal() + labs(x="Actual mtDNA CN",y="Predicted mtDNA CN",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

  dev.copy2pdf(file="../results_June2025/lin_reg_models_rand_forest/Lasso_actual_predicted_mtDNA-CN.pdf",useDingbats=F) 
  dev.off() 

  cor.test(actual,preds,method='spearman')
  cor.test(actual,preds)




## === randomForest ==========

library(randomForest) 


trainingIndex <- sample(1:nrow(lmData), 0.80*nrow(lmData)) # indices for 80% training data
trainingData <- lmData[trainingIndex,]  # training data
testData <- lmData[-trainingIndex,]  # test data


full_model = randomForest(MT_CopyNumber~.,data=lmData,na.action=na.pass)
full_model$importance 
rfimpvar = rev(sort(full_model$importance[,1]))

print(rfimpvar) 


model_rf = randomForest(MT_CopyNumber~.,data=trainingData,na.action=na.pass)
preds <- predict(model_rf, testData, na.action=na.pass)

actual = testData$MT_CopyNumber

full <- as.data.frame(cbind(actual,preds))


  ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
  theme_minimal() + labs(x="Actual mtDNA CN",y="Predicted mtDNA CN",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

  dev.copy2pdf(file="../results_June2025/lin_reg_models_rand_forest/RF_actual_predicted_mtDNA-CN.pdf",useDingbats=F) 
  dev.off() 

  cor.test(actual,preds,method='spearman')
  cor.test(actual,preds)



## ============= Prepare data with rfimpvar and plot a heatmap 

library(dplyr)

df = as.data.frame(rfimpvar)

impvar_list = vector()

i = 1
for (val in rownames(df)) {
	impvar_list[i] = val
        i = i+1 
} 

impvar_list[i] = 'MT_CopyNumber'
impvar_list[i+1] = 'DonorID'

selected_data = data %>% select(all_of(impvar_list))
write.table(selected_data,"../results_June2025/lin_reg_models_rand_forest/RF_selected_data.txt", quote=F,row.names=F)

top25features = impvar_list[c(1:25)]
selected_data_top25 = data %>% select(all_of(top25features))
write.table(selected_data_top25,"../results_June2025/lin_reg_models_rand_forest/RFtop25_selected_data.txt", quote=F,row.names=F)

final_top25 = as.data.frame(cbind(top25features,df[c(1:25),]))
colnames(final_top25) = c('Feature','Importance_score_RF')
write.table(final_top25,"../results_June2025/lin_reg_models_rand_forest/RF_top25_features.txt",quote=F,row.names=F)


rm(df)
rm(selected_data)
rm(impvar_list) 
rm(list=ls())
