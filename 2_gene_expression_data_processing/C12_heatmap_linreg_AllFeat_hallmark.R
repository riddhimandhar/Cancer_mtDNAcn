library(gplots)
library(viridis) 

## ======== p-val ========


data = read.table("results_AllFeat_pred_HallmarkExpr_LR/All_Sig_Heatmap_pval_data.txt",header=T)

row_names = data$Signature
data = data[-1]
new_data = t(data)
colnames(new_data) = row_names 

pdf(file='Hallmark_AllFeat_pred_pval.pdf',width=20, height=30, useDingbats=F)
heatmap.2(as.matrix(new_data), scale = "none", dendrogram = "none", Rowv = F, Colv = F, trace = "none", col=cividis(100,direction=-1),breaks = seq(0,0.2,0.002))
dev.off()

rm(data) 
rm(new_data)
rm(row_names)


data = read.table("results_AllFeat_pred_HallmarkExpr_LR/All_Sig_Heatmap_tval_data.txt",header=T)

row_names = data$Signature
data = data[-1]
new_data = t(data)
colnames(new_data) = row_names 

pdf(file='Hallmark_AllFeat_pred_tval.pdf', width=20, height=30, useDingbats=F)
heatmap.2(as.matrix(new_data), scale = "none", dendrogram = "none", Rowv = F, Colv = F, trace = "none", col=plasma(100,direction=1))
dev.off()

rm(data) 
rm(new_data)
rm(row_names)
