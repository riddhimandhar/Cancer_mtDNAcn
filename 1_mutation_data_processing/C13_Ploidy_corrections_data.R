## Adjusted mtDNA CN and adjusted mutational load caclulation - correction for correlation with tumor ploidy


data = read.table('../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt',header=T)


## ================= Adjusted total number of mutations w.r.t. Ploidy =================


library(ggplot2) 
library(ggpubr) 
library(dplyr)

theme_set(
  theme_minimal()
  )

theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
		
my_Theme = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 16),
  legend.title = element_text(size=14),
  plot.title = element_text(size=16, hjust=0.5)
)


my_col=rgb(0,0,0,0.3) 

x = data$Tumor_Ploidy 
y = data$Total_Num_Mutations


comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="Tumor ploidy",y="Total no. of mutations",title="Pan-Cancer") + scale_y_continuous(trans = "log10")

dev.copy2pdf(file="../results_June2025/plots/TumorPloidy_vs_TotNumMut.pdf",useDingbats=F) 
dev.off()




## === Calculation of adjusted total number of mutations === 

comb = as.data.frame(cbind(x,y,data$DonorID))

library(dplyr)
comb =  comb %>% rename(DonorID = V3)


filt_comb = comb [which(is.na(x) == 'FALSE' & is.na(y) == 'FALSE'),]


model <- lm(y ~ x) 

intcp=summary(model)[[4]][[1]]
cf=summary(model)[[4]][[2]]

num = nrow(filt_comb)

adj_mut = c()

for (i in 1:num) {
	pred = intcp + cf * as.numeric(filt_comb$x[i])
        adj_mut[i] = as.numeric(filt_comb$y[i])-pred 
}

min_adj_mut = min(adj_mut) 


for (i in 1:num) {
        adj_mut[i] = adj_mut[i]-min_adj_mut 
}


new_filtcomb = cbind(filt_comb,adj_mut) 

new_filtcomb = new_filtcomb %>% rename (Tumor_Ploidy = x) 
new_filtcomb = new_filtcomb %>% rename ( Total_Num_Mutations= y) 


rm(x)
rm(y)


### Further check 

my_col=rgb(0,0,0,0.3) 

x = as.numeric(new_filtcomb$Tumor_Ploidy)
y = as.numeric(new_filtcomb$adj_mut)

comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="Tumor ploidy",y="Adj total no. of mutations",title="Pan-Cancer") + scale_y_continuous(trans = "log10")
 

dev.copy2pdf(file="../results_June2025/plots/TumorPloidy_vs_adjusted_TotNumMut.pdf",useDingbats=F) 
dev.off()


# Right join using dplyr

new_data1 <- right_join(data, new_filtcomb[,3:4], by = "DonorID")

rm(comb)
rm(x)
rm(y)
rm(new_filtcomb) 



## ## ================= Calculation adjusted MT_CopyNumber w.r.t. Ploidy =================

my_col=rgb(0,0,0,0.3) 

x = data$Tumor_Ploidy 
y = data$MT_CopyNumber


comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="Tumor ploidy",y="MT_CopyNumber",title="Pan-Cancer") + scale_y_continuous(trans = "log10")

dev.copy2pdf(file="../results_June2025/plots/TumorPloidy_vs_MT_CopyNumber.pdf",useDingbats=F) 
dev.off()




## ============= Calculation of adjusted MT_CopyNumber ============

comb = as.data.frame(cbind(x,y,data$DonorID))

library(dplyr)
comb =  comb %>% rename(DonorID = V3)


filt_comb = comb [which(is.na(x) == 'FALSE' & is.na(y) == 'FALSE'),]


model <- lm(y ~ x) 

intcp=summary(model)[[4]][[1]]
cf=summary(model)[[4]][[2]]

num = nrow(filt_comb)

adj_mtcn = c()

for (i in 1:num) {
	pred = intcp + cf * as.numeric(filt_comb$x[i])
        adj_mtcn[i] = as.numeric(filt_comb$y[i])-pred 
}

min_adj_mtcn = min(adj_mtcn) 


for (i in 1:num) {
        adj_mtcn[i] = adj_mtcn[i]-min_adj_mtcn 
}


new_filtcomb = cbind(filt_comb,adj_mtcn) 

new_filtcomb = new_filtcomb %>% rename (Tumor_Ploidy = x) 
new_filtcomb = new_filtcomb %>% rename ( MT_CopyNumber = y) 

rm(x)
rm(y)


### Further check 

my_col=rgb(0,0,0,0.3) 

x = as.numeric(new_filtcomb$Tumor_Ploidy)
y = as.numeric(new_filtcomb$adj_mtcn)

comb = cbind(x,y)

ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="Tumor ploidy",y="Adj MT_CopyNumber",title="Pan-Cancer") + scale_y_continuous(trans = "log10")
 

dev.copy2pdf(file="../results_June2025/plots/TumorPloidy_vs_adjusted_mtCN.pdf",useDingbats=F) 
dev.off()


new_data2 <- right_join(new_data1, new_filtcomb[,3:4], by = "DonorID")


write.table(new_data2, '../results_June2025/Collated_PurFilt_DonorWise_data_with_norms_ploidyAdj.txt',quote=F,row.names=F)


rm(comb)
rm(x)
rm(y)
rm(new_filtcomb) 
rm(data)
rm(new_data1)
rm(new_data2)
