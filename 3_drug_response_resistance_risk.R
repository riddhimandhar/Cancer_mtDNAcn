data = read.table("../results_June2025/Collated_DonorWise_data_with_norms.txt",header=T)

table(data$first_therapy_type)

chemo = data[which(data$first_therapy_type == 'chemotherapy' & data$first_therapy_response != 'unknown'),]


## ========================== Chemo ========================================= 

chemo_cmpresp = chemo[which(chemo$first_therapy_response == 'complete_response'),] 
chemo_stabled = chemo[which(chemo$first_therapy_response == 'stable_disease'),]
chemo_disprog = chemo[which(chemo$first_therapy_response == 'disease_progression'),]


boxplot(chemo_cmpresp$norm_mtCN,chemo_stabled$norm_mtCN,chemo_disprog$norm_mtCN,log='y',boxwex=0.6,pch=20,lwd=2)
hist(chemo$norm_mtCN)


ncmp = dim(chemo_cmpresp)[1]
nstd = dim(chemo_stabled)[1]
ndpg = dim(chemo_disprog)[1]


m=range(na.omit(chemo$norm_mtCN))
lolim=0
uplim=max(na.omit(chemo$norm_mtCN))

int=(uplim-lolim)/8

ourbreaks=c(m[1],seq(from=lolim,to=uplim-int,by=int),m[2])

sort.chemo=chemo[order(chemo$norm_mtCN),]


bins.chemo <- cut(sort.chemo$norm_mtCN,breaks=ourbreaks,include.lowest=T)
level_bins1 <- levels(bins.chemo)

for(i in 1:length(level_bins1)) {   
  assign(paste0("chemo.bins_", i),
         sort.chemo[bins.chemo == levels(bins.chemo)[i], ])
  l1 = sort.chemo[bins.chemo == levels(bins.chemo)[i], ] 
  cE =  length(which(l1$first_therapy_response == 'disease_progression') == TRUE)
  c1 =  length(which(l1$first_therapy_response == 'stable_disease') == TRUE)
  c2 =  length(which(l1$first_therapy_response == 'complete_response') == TRUE)
  
  res_risk = (0.1+cE/ndpg)/(0.1+(c1+c2)/(nstd+ncmp))

  print(levels(bins.chemo)[i])
  print(res_risk) 
  print(cE)
  print(c1)
  print(c2)

  rm(l1)
  rm(c1)
  rm(c2)
  rm(cE)
  rm(res_risk)
}



rm(list=ls())

