library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)


data = read.table('../../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt',header=T)

inputData <- data

trainingIndex <- sample(1:nrow(inputData), 0.80*nrow(inputData)) # indices for 80% training data
trainingData <- inputData[trainingIndex,]  # training data
testData <- inputData[-trainingIndex,]  # test data


lmmod <- lm(Total_Num_Mutations ~ MT_CopyNumber, data=trainingData)
summary(lmmod) ## Performance on training data 

preds <- predict(lmmod, testData, na.action=na.pass)  # predict on test data

actual <- testData$EGFR
full <- as.data.frame(cbind(actual,preds))

actual=full$actual[which(full$preds!="NA")]
preds=full$preds[which(full$preds!="NA")]

ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0,1,0.7), size = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
theme_minimal() + labs(x="Actual num mutations",y="Predicted num mutations",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  
#geom_abline(intercept = 0, slope = 1, color = "brown", linetype = "dashed")

dev.copy2pdf(file="Actual_predicted_TotNumMut_mtDNAcn.pdf",useDingbats=F) 
dev.off() 

cor.test(actual,preds,method='spearman')



