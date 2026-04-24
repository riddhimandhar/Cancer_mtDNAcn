library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)


data = read.table('../results_June2025/Collated_PurFilt_DonorWise_data_with_norms_ploidyAdj.txt',header=T)

filt_data = data[,c(9,13:18,21:65,71:77,80:110,125:126)] 

inputData <- filt_data[,which(apply(filt_data, 2, var, na.rm=TRUE) != 0)]


encoded_gender <- model.matrix(~donor_sex-1, data=data)
inputData = cbind(inputData,encoded_gender)


lmData = inputData[which((is.na(inputData$adj_mtcn) == 'FALSE') & (is.na(inputData$adj_mut) == 'FALSE') & (is.na(inputData$donor_age_at_diagnosis_years) == 'FALSE')),]

lmData = inputData[which(is.na(inputData$adj_mtcn) == 'FALSE'),] 
lmData = lmData[,which(colSums(abs(na.omit(inputData)))!=0)]
lmData = lmData[,which(colSums(is.na(lmData))<0.9*nrow(lmData))]

lmData = lmData[which(rowSums(is.na(lmData))<0.2*ncol(lmData)),]

lmData = lmData[,which(colSums(lmData==0)<0.8*nrow(lmData))]


### =========== Predictive modeling =======================

## Ridge and Lasso regression 


library(ridge)
library(glmnet)


linRidgeMod <- linearRidge(adj_mtcn~., data = lmData)
mod_sum=summary(linRidgeMod)

names(mod_sum)

npc=mod_sum$chosen.nPCs

head(mod_sum$summaries[[npc]]$coefficients)

coefs = mod_sum$summaries[[npc]]$coefficients[,1]
pval= mod_sum$summaries[[npc]]$coefficients[,5]

key_variables=names(which(mod_sum$summaries[[npc]]$coefficients[,5]<0.05))
print(key_variables) 

selected = as.data.frame(coefs[names(coefs) %in% key_variables])
selected = cbind(key_variables,selected)
colnames(selected) = c('Feature', 'coefficient')
selected = selected[order(abs(selected$coefficient),decreasing=T),]

write.table(selected,"../results_June2025/lin_reg_models_rand_forest_ploidyAdj/Ridge_key_variables.txt",quote=F,row.names=F)
rm(selected) 


ridge_inputData=inputData[,which(colnames(inputData) %in% key_variables)]
ridge_inputData=cbind(inputData$adj_mtcn,ridge_inputData)
ridge_inputData=rename(ridge_inputData, 'adj_mtcn'='inputData$adj_mtcn')


cor_p = vector()
cor_s = vector()
pval_p = vector()
pval_s = vector() 
fl = 0 

for (ct in 1:100) { 

  ridge_trainingIndex <- sample(1:nrow(ridge_inputData), 0.80*nrow(ridge_inputData)) # indices for 80% training data
  ridge_trainingData <- ridge_inputData[ridge_trainingIndex,]  # training data
  ridge_testData <- ridge_inputData[-ridge_trainingIndex,]  # test data

  lmmod_ridge <- lm(adj_mtcn ~ . , data=ridge_trainingData)

  ## summary(lmmod_ridge) ## Performance on training data 

  expv=summary(lmmod_ridge)[9][[1]]
  preds <- predict(lmmod_ridge, ridge_testData, na.action=na.pass)  # predict on test data

  actual <- ridge_testData$adj_mtcn
  full <- as.data.frame(cbind(actual,preds))

  actual=full$actual[which(full$preds!="NA")]
  preds=full$preds[which(full$preds!="NA")]

  rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))
  tss <- sum(na.omit(actual - mean(na.omit(actual))) ^ 2/length(na.omit(actual-mean(na.omit(actual)))))
  rsq <- 1 - rss/tss

  # rsq  ## Rsq on test data 

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
     pdf(file="../results_June2025/lin_reg_models_rand_forest_ploidyAdj/Ridge_actual_predicted_adj_mtDNA-CN.pdf",useDingbats=F) 
     gp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
          theme_minimal() + labs(x="Actual mtDNA CN",y="Predicted mtDNA CN",title= tt) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  
     print(gp) 
     dev.off()
     rm(gp) 
     rm(tt)
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



### ==== Lasso ==========


len=length(names(lmData))
x = model.matrix(adj_mtcn~., lmData)[,-1]
y = lmData %>% select(adj_mtcn) %>% unlist() %>% as.numeric()
grid = 10^seq(10, -2, length = 100)

cv.out = cv.glmnet(x, y, alpha = 1)
bestlam = cv.out$lambda.min

out = glmnet(x, y, alpha = 1, lambda = grid)

lasso_coef = coef(out, s = bestlam)[1:len,] 


## New dataset with selected features 

impvar=lasso_coef[lasso_coef != 0]
print(impvar)
impvarlist=names(impvar[-1])
implen=length(names(impvar))

selected = as.data.frame(impvar[names(impvar) %in% impvarlist])
selected = cbind(impvarlist,selected)
colnames(selected) = c('Feature', 'coefficient')
selected = selected[order(abs(selected$coefficient),decreasing=T),]

write.table(selected,"../results_June2025/lin_reg_models_rand_forest_ploidyAdj/Lasso_key_variables.txt",quote=F,row.names=F)
rm(selected) 

lasso_inputData=inputData[,which(colnames(inputData) %in% impvarlist)]
lasso_inputData=cbind(inputData$adj_mtcn,lasso_inputData)
lasso_inputData=rename(lasso_inputData, 'adj_mtcn'='inputData$adj_mtcn')


## Training and testing 

cor_p = vector()
cor_s = vector()
pval_p = vector()
pval_s = vector() 
fl = 0 

for (ct in 1:100) { 

  lasso_trainingIndex <- sample(1:nrow(lasso_inputData), 0.80*nrow(lasso_inputData)) # indices for 80% training data
  lasso_trainingData <- lasso_inputData[lasso_trainingIndex,] # training data
  lasso_testData <- lasso_inputData[-lasso_trainingIndex,] # test data

  lmmod_lasso <- lm(adj_mtcn ~ . , data=lasso_trainingData)
  expv=summary(lmmod_lasso)[9][[1]]
  ## expv
  
  preds <- predict(lmmod_lasso, lasso_testData, na.action=na.pass)  # predict on test data
  actual <- lasso_testData$adj_mtcn
  full <- as.data.frame(cbind(actual,preds))

  actual=full$actual[which(full$preds!="NA")]
  preds=full$preds[which(full$preds!="NA")]
  rss <- sum(na.omit(preds - actual) ^ 2)/length(na.omit(preds-actual))
  tss <- sum(na.omit(actual - mean(na.omit(actual))) ^ 2/length(na.omit(actual-mean(na.omit(actual)))))
  rsq <- 1 - rss/tss

  ## rsq
 
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
    pdf(file="../results_June2025/lin_reg_models_rand_forest_ploidyAdj/Lasso_actual_predicted_adj_mtDNA-CN.pdf",useDingbats=F) 
    gp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
         theme_minimal() + labs(x="Actual mtDNA CN",y="Predicted mtDNA CN",title= tt) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  
    print(gp)
    dev.off()
    rm(gp)
    rm(tt)
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
 


## === randomForest ==========


library(randomForest) 

full_model = randomForest(adj_mtcn~.,data=lmData, na.action=na.pass, importance=T)
full_model$importance 
rfimpvar = rev(sort(full_model$importance[,2]))
print(rfimpvar) 
print(full_model$importanceSD)

cor_p = vector()
cor_s = vector()
pval_p = vector()
pval_s = vector() 
fl = 0 

for (ct in 1:100) { 

  trainingIndex <- sample(1:nrow(lmData), 0.80*nrow(lmData)) # indices for 80% training data
  trainingData <- lmData[trainingIndex,]  # training data
  testData <- lmData[-trainingIndex,]  # test data

  model_rf = randomForest(adj_mtcn~.,data=trainingData,na.action=na.pass)
  preds <- predict(model_rf, testData, na.action=na.pass)

  actual = testData$adj_mtcn

  full <- as.data.frame(cbind(actual,preds))

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
     pdf(file="../results_June2025/lin_reg_models_rand_forest_ploidyAdj/RF_actual_predicted_adj_mtDNA-CN.pdf",useDingbats=F)  
     gp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
     theme_minimal() + labs(x="Actual adj mtDNA CN",y="Predicted adj mtDNA CN",title= tt) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  
     print(gp) 
     dev.off() 
     rm(gp) 
     rm(tt)
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
 


## ============= Prepare data with rfimpvar and plot a heatmap 


library(dplyr)

df = as.data.frame(rfimpvar)

impvar_list = rownames(df)
impvar_list = append(impvar_list,'adj_mtcn')  # 'MT_CopyNumber'
#impvar_list = append(impvar_list,'DonorID') 

selected_data = inputData %>% select(all_of(impvar_list))
write.table(selected_data,"../results_June2025/lin_reg_models_rand_forest_ploidyAdj/RF_selected_data.txt", quote=F,row.names=F)

top25features = impvar_list[c(1:25)]
selected_data_top25 = data %>% select(all_of(top25features))
write.table(selected_data_top25,"../results_June2025/lin_reg_models_rand_forest_ploidyAdj/RFtop25_selected_data.txt", quote=F,row.names=F)

final_top25 = as.data.frame(cbind(top25features,df[c(1:25),]))
colnames(final_top25) = c('Feature','Importance_score_RF')
write.table(final_top25,"../results_June2025/lin_reg_models_rand_forest_ploidyAdj/RF_top25_features.txt",quote=F,row.names=F)


rm(df)
rm(selected_data)
rm(impvar_list) 
rm(list=ls())
