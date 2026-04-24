## ========== using all factors to predict EMT hallmark gene expression Linear Regression ==== 

library(pls)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Dict)
library(ridge)


data = read.csv('../../results_joint_MutExpr/Collated_PurFilt_DonorWise_data_norms_stage_encoding_withExpr.csv',header=T)
dim(data)


hmfile = "../../../../datasets/hallmarks_msigdb/h.all.v2025.1.Hs.symbols.gmt"
sigline <- readLines(hmfile)
#sig = sigline[14] ## EMT hallmark 

newlist1 = list()
newlist2 = list()
newfeatlist = vector()

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


 oth_feature_list=c('DonorID','donor_age_at_diagnosis_years','donor_sex','MT_CopyNumber','Tumor_Purity','Tumor_Ploidy',
                'GenesCdsPromsAff_MtLoc','GenesStrMut_MtLoc','GenesWeakMut_MtLoc','GenesTranscrEff_MtLoc',
		'NumGenesCdsPromsAff','NumGenesStrMut','NumGenesWeakMut','NumGenesTranscrEff',
                'PercSNP','PercINS','PercDEL','MT_SNV','MT_INDEL','NUC_SNV','NUC_INDEL')

 mtdna_genes = c('MT.ND1','MT.ATP6','MT.CO1','MT.CO2','MT.CYB')
 full_var_list = c(oth_feature_list,genelist,mtdna_genes) 

 new_data =data[which(names(data) %in% full_var_list)] 
 colnames(new_data) 
 
 new_data = new_data[-1]

 rm(full_var_list)
   

 subdir = elems[[1]][1] 
 print(subdir)
 ## dir.create(file.path('results_AllFeat_pred_HallmarkExpr_LR',subdir ), showWarnings = FALSE)

 featlist = list()
 featlist2 = list() 
 
   for (k in 1:length(genelist)) { 
   
        if(!genelist[k] %in% colnames(new_data)) { 
           print(genelist[k])
        }
	if(genelist[k] %in% colnames(new_data)) { 
           inputData = data.frame(cbind(new_data[names(new_data) %in% c(oth_feature_list,mtdna_genes)],new_data[genelist[k]]))
	   col_id = dim(inputData)[2]	
	   colnames(inputData)[col_id] = 'gene'	

           lmmod <- lm(gene ~ ., data=inputData)  ## lm


           for (ct in 2:col_id) { 
               
               feat = rownames(summary(lmmod)[[4]])[ct]  ## for 'lm'
   	       tval = summary(lmmod)[[4]][col_id*2+ct]   ## for 'lm'
	       pval = summary(lmmod)[[4]][col_id*3+ct]   ## for 'lm'
            
	       # print(c(feat,tval,pval))

               val1 = vector()
               val2 = vector()
               val1 = featlist[[feat]]
               val2 = featlist2[[feat]] 
               featlist[[feat]] = c(val1,abs(tval))
               featlist2[[feat]] = c(val2,pval)
               rm(val1)
               rm(val2)
               rm(feat)
               rm(tval)
               rm(pval) 
           }
           rm(inputData)
           ## rm(ridgemod)
           rm(lmmod)
           rm(col_id) 
        }
    }



 sink(paste('results_AllFeat_pred_HallmarkExpr_LR/',subdir,'_tval_Features.txt',sep=''))
 
 for (key in names(featlist)) {
    value <- featlist[[key]]
    crt = vector()
    crt = newlist1[[elems[[1]][1]]]
    newlist1[[elems[[1]][1]]] = c(crt,median(value))

    if(subdir == 'HALLMARK_ADIPOGENESIS') {  
      newfeatlist = append(newfeatlist,key) 
    }

    cat(paste(key,"\t"))
    cat(value)
    cat("\n")
    rm(value)
    rm(crt) 
 }
 sink() 
 rm(key)

 pdf(file = paste('results_AllFeat_pred_HallmarkExpr_LR/plots/',subdir,'_tval_Features.pdf',sep=''),useDingbats=F) 
 boxplot(featlist[[1]],featlist[[2]],featlist[[3]],featlist[[4]],featlist[[5]],featlist[[6]],featlist[[7]],featlist[[8]],
        featlist[[9]],featlist[[10]],featlist[[11]],featlist[[12]],featlist[[13]],featlist[[14]],featlist[[15]],featlist[[16]],
        featlist[[17]],featlist[[18]],featlist[[19]],featlist[[20]],featlist[[21]],featlist[[22]],featlist[[23]],featlist[[24]],
        featlist[[25]],col=rgb(0,0.5,0.7,1),pch=20)
 
 dev.off() 

 sink(paste('results_AllFeat_pred_HallmarkExpr_LR/',subdir,'_pval_Features.txt',sep=''))
 
 for (key in names(featlist2)) {
    value <- featlist2[[key]]
    crt = vector()
    crt = newlist2[[elems[[1]][1]]]
    newlist2[[elems[[1]][1]]] = c(crt,median(value))

    cat(paste(key,"\t"))
    cat(value)
    cat("\n")
    rm(crt) 
    rm(value)
 }
 sink() 
 rm(key)
 

 pdf(file = paste('results_AllFeat_pred_HallmarkExpr_LR/plots/',subdir,'_pval_Features.pdf',sep=''),useDingbats=F) 
 boxplot(featlist2[[1]],featlist2[[2]],featlist2[[3]],featlist2[[4]],featlist2[[5]],featlist2[[6]],featlist2[[7]],featlist2[[8]],
        featlist2[[9]],featlist2[[10]],featlist2[[11]],featlist2[[12]],featlist2[[13]],featlist2[[14]],featlist2[[15]],featlist2[[16]],
        featlist2[[17]],featlist2[[18]],featlist2[[19]],featlist2[[20]],featlist2[[21]],featlist2[[22]],featlist2[[23]],featlist2[[24]],
        featlist2[[25]], log='y',col=rgb(0,0.5,0.7,1),pch=20)
 
 dev.off() 
 
 rm(featlist)
 rm(featlist2) 
 rm(subdir)

 #boxplot(featlist2[[1]],featlist2[[2]],featlist2[[3]],featlist2[[4]],featlist2[[5]],featlist2[[13]],featlist2[[23]],log='y')
 #boxplot(featlist[[1]],featlist[[2]],featlist[[3]],featlist[[4]],featlist[[5]],featlist[[13]],featlist[[23]])	

}

### print newfeatlist, newlist1, and newlist2 


sink(paste('results_AllFeat_pred_HallmarkExpr_LR/All_Sig_Heatmap_tval_data.txt',sep=''))

cat('Signature\t')
cat(newfeatlist)
cat('\n')

for (key in names(newlist1)) { 
   cat(paste(key,'\t'))
   value = newlist1[[key]]
   cat(value)
   cat('\n')
   rm(value)
}
sink()
rm(key)


sink(paste('results_AllFeat_pred_HallmarkExpr_LR/All_Sig_Heatmap_pval_data.txt',sep=''))

cat('Signature\t')
cat(newfeatlist)
cat('\n')


for (key in names(newlist2)) { 
   cat(paste(key,'\t'))
   value = newlist2[[key]]
   cat(value)
   cat('\n')
   rm(value) 
}

sink()
rm(key)

rm(elems)
rm(new_data)
rm(hmfile)
rm(data)
rm(sig)
rm(sigline)
rm(newfeatlist)
rm(newlist1)
rm(newlist2) 

