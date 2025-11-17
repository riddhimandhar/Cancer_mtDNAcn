library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr) 

theme_set(
  theme_minimal()
  )
	
my_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 16),
  legend.title = element_text(size=14),
  plot.title = element_text(size=16, hjust=0.5)
)


data = read.table("../results_June2025/PCAWG_MitoData_collected.txt",header=T)

my_col=rgb(0,0.5,0.8,0.6) 


## =========== Normal tissue =================

filt_data = data[is.na(data$Normal_Tissue) == 0,]

head(filt_data) 




canc = which(table(filt_data$histology_abbreviation)> 2)

filt_data = filt_data[which(filt_data$histology_abbreviation %in% names(canc)),]

# unique(filt_data$histology_abbreviation)) 

table(filt_data$histology_abbreviation)

# Biliary-AdenoCA  Breast-AdenoCA     Eso-AdenoCA       Liver-HCC    Panc-AdenoCA 
#             15              10              10               8             158 
# Panc-Endocrine Stomach-AdenoCA 
#             80              37 



ggplot(filt_data,aes(x=histology_abbreviation,y=Ratio_MTCopyNumber_CancerToNormal)) + geom_boxplot(fill = 'dodgerblue1') + 
      scale_fill_viridis(discrete = T, direction = -1, option = "H") +  theme_minimal() + scale_y_continuous(trans='log10') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      geom_hline(yintercept=1,col='blue',linetype='dashed')
	
dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Boxplot_tissue_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()



#plot(filt_data$MT_CopyNumber,filt_data$Ratio_MTCopyNumber_CancerToNormal,pch=20,log='xy',col=rgb(0,0.5,0.8,0.7))

x = filt_data$MT_CopyNumber
y = filt_data$Ratio_MTCopyNumber_CancerToNormal
comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="MT_CopyNumber",y="Ratio_MTCopyNumber_CancerToNormal",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Corr_tissue_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()

cor.test(filt_data$MT_CopyNumber,filt_data$Ratio_MTCopyNumber_CancerToNormal,method='spearman')

## rs = 0.37, p = 1.33e-09
## rp = 0.33, p = 7.81e-08

rm(filt_data) 
rm(canc)
rm(x)
rm(y)
rm(comb) 



## ====================== Normal blood =================

filt_data2 = data[is.na(data$Normal_Blood) == 0,]

head(filt_data2) 

canc = which(table(filt_data2$histology_abbreviation)> 2)

filt_data2 = filt_data2[which(filt_data2$histology_abbreviation %in% names(canc)),]

# unique(filt_data2$histology_abbreviation)) 

table(filt_data2$histology_abbreviation)

# Biliary-AdenoCA      Bone-Benign       Bone-Epith   Bone-Osteosarc 
#              17               14                9               36 
#  Breast-AdenoCA      Breast-LobularCA      CNS-Medullo 
#             102                7              146 
#   CNS-PiloAstro      Eso-AdenoCA         Head-SCC       Kidney-RCC 
#              89               88               13               95 
#       Liver-HCC       Lymph-BNHL        Lymph-CLL      Myeloid-AML 
#             256               83              100                6 
#     Myeloid-MPN    Ovary-AdenoCA     Panc-AdenoCA 
#              22               38               51 
#  Panc-Endocrine    Prost-AdenoCA    Skin-Melanoma 
#               5              182               70 


ggplot(filt_data2,aes(x=histology_abbreviation,y=Ratio_MTCopyNumber_CancerToNormal)) + geom_boxplot(fill = 'dodgerblue1') + 
      scale_fill_viridis(discrete = T, direction = -1, option = "H") +  theme_minimal() + scale_y_continuous(trans='log10') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      geom_hline(yintercept=1,col='blue',linetype='dashed')
	
dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Boxplot_blood_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()



##plot(filt_data2$MT_CopyNumber,filt_data2$Ratio_MTCopyNumber_CancerToNormal,pch=20,log='xy',col=rgb(0,0.5,0.8,0.7))

x = filt_data2$MT_CopyNumber
y = filt_data2$Ratio_MTCopyNumber_CancerToNormal
comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="MT_CopyNumber",y="Ratio_MTCopyNumber_CancerToNormal",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Corr_blood_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()

cor.test(filt_data2$MT_CopyNumber,filt_data2$Ratio_MTCopyNumber_CancerToNormal,method='spearman')

## rs = 0.74, p < 10-15
## rp = 0.63, p < 10-15 


rm(filt_data2) 
rm(canc)
rm(x)
rm(y) 
rm(comb)


## ====================== Normal other =================

filt_data3 = data[is.na(data$Normal_Other) == 0,]

head(filt_data3) 


table(filt_data3$histology_abbreviation)

# Breast-AdenoCA     Lymph-BNHL    Myeloid-AML  Ovary-AdenoCA   Panc-AdenoCA 
#             1             17              9              9             22 
# Prost-AdenoCA 
#            10 


canc = which(table(filt_data3$histology_abbreviation)> 2)

filt_data3 = filt_data3[which(filt_data3$histology_abbreviation %in% names(canc)),]

# unique(filt_data3$histology_abbreviation)) 



ggplot(filt_data3,aes(x=histology_abbreviation,y=Ratio_MTCopyNumber_CancerToNormal)) + geom_boxplot(fill = 'dodgerblue1') + 
      scale_fill_viridis(discrete = T, direction = -1, option = "H") +  theme_minimal() + scale_y_continuous(trans='log10') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      geom_hline(yintercept=1,col='blue',linetype='dashed')
	
dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Boxplot_other_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()



#plot(filt_data3$MT_CopyNumber,filt_data3$Ratio_MTCopyNumber_CancerToNormal,pch=20,log='xy',col=rgb(0,0.5,0.8,0.7))

x = filt_data3$MT_CopyNumber
y = filt_data3$Ratio_MTCopyNumber_CancerToNormal
comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="MT_CopyNumber",y="Ratio_MTCopyNumber_CancerToNormal",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

dev.copy2pdf(file='../results_June2025/cancer_normal_mtCN_plots/Corr_other_mtCN_ratio.pdf',useDingbats=F,width=8,height=6) 
dev.off()

cor.test(filt_data3$MT_CopyNumber,filt_data3$Ratio_MTCopyNumber_CancerToNormal,method='spearman')

## rs = 0.19, p = 0.18
## rp = 0.38, p = 0.005 

rm(filt_data3) 
rm(canc)
rm(x)
rm(y)
rm(comb) 


## --------------

rm(data) 


