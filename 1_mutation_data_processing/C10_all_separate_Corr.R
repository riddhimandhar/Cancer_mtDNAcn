library(ggplot2) 
library(ggpubr)
library(gridExtra)
 

theme_set(
  theme_minimal()
  )

	
my_Theme = theme(
  axis.title.x = element_text(size = 9),
  axis.text.x = element_text(size = 7),
  axis.text.y = element_text(size = 7),
  axis.title.y = element_text(size = 9),
  legend.title = element_text(size=9),
  plot.title = element_text(size=9, hjust=0.8)
)


load_data = read.table("../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt",header=T)
list = names(which(table(load_data$Cancer)>=3))
data = load_data[load_data$Cancer %in% list,]

newdata = data[,c(1,2,5:128)]  


### === MT-CN vs TotNumMutations ============= 

my_col= '#CC3333'   ## '#009933'

fullplot = list()
i = 1

for (canctype in sort(unique(data$Cancer))){

	tmp =  newdata[which(newdata$Cancer == canctype),]
        x = tmp$Total_Num_Mutations
        y = tmp$MT_CopyNumber
 	n1 = length(na.omit(x))
        n2 = length(na.omit(y)) 
       
        if (n1>=3 & n2 >=3) {
         rv=cor.test(x,y,method="spearman")
         rves=rv$estimate
         rvp=rv$p.value
         rvm = paste(canctype,paste('rs',round(rves,2),sep='= '),paste('p',round(rvp,2),sep='='),sep=', ')

        comb = cbind(x,y)

        tmp_plot = ggplot(comb, aes(x=x, y=y)) + geom_point(color=my_col, size = 1) + geom_smooth(method=lm, color="#660033", linewidth=1) +  ## "#003300"
        my_Theme + labs(x="Total no. of mutations",y="mtDNA CN",title= rvm) + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10')
	## theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))     

        fullplot[[i]] = tmp_plot
        i = i+1 

        rm(rv)
	rm(comb)
        rm(rves)
        rm(rvp)
        rm(rvm)
       }
       rm(n1)
       rm(n2)
       rm(x)
       rm(y)
       rm(tmp) 

}


ml = gridExtra::marrangeGrob(fullplot,layout_matrix = matrix(seq_len(3*3),nrow=3,ncol=3,byrow=T))

ggsave("../results_June2025/individual_corr_plots/Corr_MT_CN_vs_TotNumMut.pdf",ml)

dev.off() 


