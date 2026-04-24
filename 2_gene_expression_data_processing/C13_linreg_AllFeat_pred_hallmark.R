## ========== using all factors to predict EMT hallmark gene expression Linear Regression ==== 

library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Dict)
library(ridge)


data = read.csv('../../results_joint_MutExpr/Collated_PurFilt_DonorWise_data_norms_stage_encoding_withExpr.csv',header=T)
dim(data)
data$MT_SNV[is.na(data$MT_SNV)] <- 0
data$MT_INDEL[is.na(data$MT_INDEL)] <- 0


hmfile = "../../../../datasets/hallmarks_msigdb/h.all.v2025.1.Hs.symbols.gmt"
sigline <- readLines(hmfile)
#sig = sigline[14] ## EMT hallmark 

sink('results_AllFeat_pred_HallmarkExpr_LR/1_Hallmark_wise_cor_pval.txt') 
cat(paste('Hallmark','Mean_cor','Sd_cor','Num_sig_pval','Tot_genes','Perc','\n',sep=' '))
sink()

CON <- file("results_AllFeat_pred_HallmarkExpr_LR/1_Hallmark_wise_cor_pval.txt", "a")



for (sig in sigline) {

 elems = strsplit(sig,"\t") 
 num_elems = length(elems[[1]]) 
 #print(noquote(paste(elems[[1]][1],num_elems,sep=' ')))
 genelist = vector() 
 for (i in 3:num_elems) {
 	genelist = append(genelist,noquote(elems[[1]][i]))  
 }

 print(elems[[1]][1])
 print(genelist) 


 oth_feature_list=c('DonorID','donor_age_at_diagnosis_years','MT_CopyNumber','Tumor_Purity','Tumor_Ploidy',
                'GenesCdsPromsAff_MtLoc','GenesStrMut_MtLoc','GenesWeakMut_MtLoc','GenesTranscrEff_MtLoc',
		'NumGenesCdsPromsAff','NumGenesStrMut','NumGenesWeakMut','NumGenesTranscrEff',
                'PercSNP','PercINS','PercDEL','MT_SNV','MT_INDEL','NUC_SNV','NUC_INDEL','donor_sex','donor_sexmale')

 mtdna_genes = c('MT.ND1','MT.ATP6','MT.CO1','MT.CO2','MT.CYB')
 full_var_list = c(oth_feature_list,genelist,mtdna_genes) 

 new_data =data[which(names(data) %in% full_var_list)] 
 colnames(new_data) 
 
 encoded_gender <- model.matrix(~donor_sex-1, data=new_data)
 new_data = new_data[-c(1:2)]
 new_data = cbind(encoded_gender,new_data)

 rm(full_var_list)
   

 subdir = elems[[1]][1] 
 print(subdir)

 dir.create(file.path('results_AllFeat_pred_HallmarkExpr_LR',subdir ), showWarnings = FALSE)
 sink(paste('results_AllFeat_pred_HallmarkExpr_LR/',subdir,'/','0_Prediction_correlation_results.txt',sep=''))
  
 corlist = vector()
 pvlist = 0 
 
   for (k in 1:length(genelist)) { 
   
        ## if(!genelist[k] %in% colnames(new_data)) { 
        ##     print(genelist[k])
        ## }

	if(genelist[k] %in% colnames(new_data)) { 
           inputData = data.frame(cbind(new_data[names(new_data) %in% c(oth_feature_list,mtdna_genes)],new_data[genelist[k]]))
	   col_id = dim(inputData)[2]	
           name = colnames(inputData)[col_id]
	   colnames(inputData)[col_id] = 'gene'	

         
           lmmod <- lm(gene ~ ., data=inputData)
           key_variables = vector() ## for 'lm'
           for (ct in 2:col_id) { 
                  if(summary(lmmod)[[4]][col_id*3+ct] < 0.05) {
                  key_variables = append(key_variables,rownames(summary(lmmod)[[4]])[ct])
               }
           }
           if(length(key_variables) == 0) { 
              key_variables = rownames(summary(lmmod)[[4]])[2]
           }
                      

           select_inputData=as.data.frame(inputData[,which(colnames(inputData) %in% key_variables)])
           select_inputData=cbind(inputData$gene,select_inputData)
           select_inputData=rename(select_inputData, 'gene'='inputData$gene')
	
           ilist1 = vector() 
           ilist2 = vector() 
  
           for(iter in 1:10)  {
    
             trainingIndex <- sample(1:nrow(select_inputData), 0.80*nrow(select_inputData)) # indices for 80% training data
    	     trainingData <- select_inputData[trainingIndex,]  # training data
	     testData <- select_inputData[-trainingIndex,]  # test data

	     train_select <- lm(gene ~ . , data=trainingData)
             preds <- predict(train_select, testData, na.action=na.pass)
             actual <- testData$gene
	     full <- as.data.frame(cbind(actual,preds))

 	     actual=full$actual[which(full$preds!="NA")]
	     preds=full$preds[which(full$preds!="NA")]
             ## print(cor.test(actual,preds))

             pv=cor.test(actual,preds,method='spearman')
             ilist1 = append(ilist1,pv$estimate)
             ilist2 = append(ilist2,pv$p.value)

	     if(iter == 1)  { 
                wrf = paste('results_AllFeat_pred_HallmarkExpr_LR/',subdir,'/',name,'_expr_pred_AllFeat.pdf',sep='')
                tt = paste(round(pv$estimate,2),pv$p.value)
	        pdf(file=wrf,useDingbats=F)
	        ggp = ggplot(full, aes(x=actual, y=preds)) + geom_point(color=rgb(0,0.5,0.8,0.7), linewidth = 1.5) + geom_smooth(method=lm, color="darkgrey", size=1) +
	        theme_minimal() + labs(x="Actual expression",y="Predicted expression",title= tt) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') 
	        print(ggp)
	        dev.off() 
             }
       	     rm(pv)  
	   }

           cat(paste(name,mean(ilist1),mean(ilist2),'\n',sep='\t'))
	
           corlist = append(corlist,mean(ilist1)) 
	   if(mean(ilist2) < 0.05) {
		pvlist = pvlist + 1 
       	   }
           

    	   rm(wrf) 
           rm(select_inputData)
	   rm(inputData)
	   rm(name)
	   rm(actual)
	   rm(preds)
	   rm(ggp) 
	   rm(pv)
	   rm(full)
	   rm(trainingData)
	   rm(testData)	
           rm(lmmod)	
           rm(col_id) 
           #rm(train_ridge)
           #rm(ridgemod)
        }

    }
    sink()
    avgcor = round(mean(corlist),2)
    sdcor = round(sd(corlist),2)

    perc = round((pvlist/length(genelist)*100),2)
    towrite = paste(elems[[1]][1],' ',avgcor,sdcor,pvlist,length(genelist),perc,'\n',sep=' ')
    print(towrite)
    ## cat(towrite, file = "results_AllFeat_pred_HallmarkExpr_LR/1_Hallmark_wise_cor_pval.txt", append = TRUE)
    cat(towrite, file = CON)

    rm(perc)
    rm(towrite)
    rm(sdcor)
    rm(avgcor)
    
    rm(elems)
    rm(subdir)
    rm(corlist)
    rm(pvlist) 

}

close(CON) 


rm(hmfile)

rm(data)
rm(sig)
rm(sigline)
rm(CON) 

rm(list=ls())

### 
