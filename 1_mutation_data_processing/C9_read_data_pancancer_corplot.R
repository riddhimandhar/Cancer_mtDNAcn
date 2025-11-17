### 

data = read.table("../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt",header=T)

table(data$Cancer)


### Correlation plots 


library(ggplot2) 
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

my_col=rgb(0,0.5,0.8,0.3) 


## == MT-CN vs TotNumMutations ==== 

x = data$Total_Num_Mutations
y = data$MT_CopyNumber

comb = cbind(x,y)


ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 2) + geom_smooth(method=lm, color="blue", size=2) +
  my_Theme + labs(x="No. of mutations",y="mtDNA CN",title="Pan-cancer") + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')  

dev.copy2pdf(file="../results_June2025/plots/panCancer_corrs/Pancancer_mtDNA-CN_vs_numMutations.pdf",useDingbats=F) 
dev.off() 

cor.test(x,y,method='spearman')


