
#### =============== PurityFilt =================== 


data = read.table("../results_June2025/Collated_DonorWise_data.txt",header=T)

data = data[which(data$Tumor_Purity >= 0.5),]

newdata = data[,c(1,2,5:126)]  



keys = list()
values = list()
values2 = list()

i = 1

for (canctype in unique(data$Cancer)){

	tmp =  newdata[which(newdata$Cancer == canctype),]
        x = tmp$MT_CopyNumber
        x2 = tmp$Total_Num_Mutations 
 	n = length(na.omit(x))
        n2 = length(na.omit(x))  
       
        if (n>=3) {
	  med_mtcn = median(na.omit(x))
        }
	if(n<3) {
 	  med_mtcn = mean(na.omit(x))
        }

	if (n2>=3) {
	  med_nummut = median(na.omit(x2)) 
        }
	if(n2<3) {
 	  med_nummut = mean(na.omit(x))
        }

        print(canctype) 
        print(med_mtcn)
        print(med_nummut)
	keys[[i]] = canctype
        values[[i]] = med_mtcn
        values2[[i]] = med_nummut 
        i = i+1

        rm(n)
        rm(x)
 	rm(x2)
	rm(n2)
       
        rm(tmp) 
}

rm(i) 


newdict = list()
newdict2 = list()

for (j in seq_along(keys)) {
        #print(keys[[j]])
        #print(values[j])
        newdict[keys[[j]]] = values[j]
        newdict2[keys[[j]]] = values2[j]
      }
 
for (j in seq(1,length(newdict2))){
  key <- newdict[j]
  #print(paste(names(key), key))
  #print(names(key))
  print(key)
}

rm(keys)
rm(values)
rm(values2) 
rm(j)


norm_mtCN = vector()
norm_nummut = vector() 

i = 1

for (i in seq(1,nrow(data))) {

           val = data$MT_CopyNumber[i]/as.numeric(newdict[data$Cancer[i]])
           val2 = data$Total_Num_Mutations[i]/as.numeric(newdict2[data$Cancer[i]])
           if(!is.na(val)) { 
             val = round(val,3)
           }
           if(!is.na(val2)) { 
             val2 = round(val2,3)
           }
 
           norm_mtCN[i] = val 
           norm_nummut[i] = val2

           rm(val)
           rm(val2) 
           #print(norm_mtCN[i])
}

data$norm_mtCN = norm_mtCN
data$norm_TotNumMut = norm_nummut 


write.table(data,file='../results_June2025/Collated_PurFilt_DonorWise_data_with_norms.txt',quote=F,row.names=F)



rm(list=ls())
