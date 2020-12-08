.libPaths("C:/Users/LCOIN/R-4.0.2/library")

library('jsonlite');
library('ggplot2')
library(tidyr)
library(reshape2)
RHOME="../../R"
setwd("../data/new_samples")
source(paste(RHOME, "plotPlasmaFuncs.R", sep="/"));

names = c("1524", "1249", "1494", "1084","065", "098");

x = 2:1999
P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x)
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x) 
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x)

.cumsum1<-function(x) sum(x) -cumsum(x)

P0_c = data.frame(cbind(apply(P0,2,cumsum), apply(P0,2,.cumsum1)))
shared_c = data.frame(cbind(apply(shared,2,cumsum), apply(shared,2,.cumsum1)))
unique_c = data.frame(cbind(apply(unique,2,cumsum), apply(unique,2,.cumsum1)))
ncol = dim(P0_c)[2]
ratios = matrix(nrow = dim(P0_c)[1],ncol = ncol)
for(i in 1:ncol){
  
  shared_c[,i] = shared_c[,i]/P0_c[,i]
  unique_c[,i] = unique_c[,i]/P0_c[,i]
  ratios[,i] = shared_c[,i]/unique_c[,i]
}
ratios = data.frame(ratios)
names(ratios) = c(paste(names(unique),"LTE",sep="."),paste(names(unique),"GT",sep="."))
probs =c(0.25,.5,.75)

ratios_ = cbind(x,.addAvg(ratios, "LTE", probs  =probs),.addAvg(ratios, "GT", probs  =probs))
#ratios = cbind(x,ratios)
df = pivot_longer(ratios_, names_to = "ID", values_to = "proportion", cols=names(ratios)[-1])  %>% 
  separate(ID, c('ID', 'type'), sep='\\.', remove = T)
ratios0 = ratios_[,c(1,grep("q_25",names(ratios_)))]
ratios1 = ratios_[,c(1,grep("q_50",names(ratios_)))]
ratios2 = ratios_[,c(1,grep("q_75",names(ratios_)))]

df0 = pivot_longer(ratios0, names_to="ID",values_to="q25",cols = names(ratios0)[-1])
df1 = pivot_longer(ratios1, names_to="ID",values_to="q50",cols = names(ratios1)[-1])
df2 = pivot_longer(ratios2, names_to="ID",values_to="q75",cols = names(ratios2)[-1])
df_all = cbind(df1,df0,df2)[,c(1,2,3,6,9)]

textsize=20
ggp<-ggplot(df_all, aes(x,q50, linetype=ID))+
  ggtitle("Ratio of shared to unique")+
  ylab("Ratio")+
  geom_line()+
  geom_ribbon(aes(ymin=q25, ymax=q75, fill = ID), alpha=0.1)+
  scale_x_continuous(trans='log10',
                     breaks = c(2,5,10,20,50,100,200,500,1000,2000))+
  ylim(0,0.4)+
  theme_bw()+
  theme(text = element_text(size=textsize))
ggp

ggsave("ratio.png", plot=ggp, width = 30, height = 30, units = "cm")

ggp1 = ggplot(df, aes(x,proportion, fill=ID,color=ID, linetype=type))+
  geom_line()+
  scale_x_continuous(trans='log10')+
  ggtitle("Ratio of shared to unique")+
  theme_bw()+
  theme(text = element_text(size=textsize))
ggp1

ggsave("ratio1.png", plot=ggp1, width = 30, height = 30, units = "cm")

ncol1 = dim(P0)[2]
for(i in 1:ncol1){
  P0[,i] = P0[,i]/sum(P0[,i])
}

#density  = data.frame(P0)
#density = cbind(x,density[,!(names(density) %in% "combined")])
###DENSITY PLOTS

x = 2:1999
P0 = readDataAll(names[1:4], '_P0_read_length_count.json', x,"Tumour_all")
shared = readDataAll(names[1:4], '_P0.shared.somatic_reads_read_length_count.json', x,"Tumour_shared") 
unique = readDataAll(names[1:4], '_P0.unique.somatic_reads_read_length_count.json', x, "Tumour_unique")
P01 = readDataAll(names[5:6], '_P0_read_length_count.json', x, "Benign_all")
merged = readDataAll(names[5:6], '_P0.somatic_reads_read_length_count.json', x, "Benign_filtered")
mergedT = readDataAll(names[1:4], '_P0.somatic_reads_read_length_count.json', x,"Tumour_filtered")

P0_d = data.frame(apply(P0, 2, makeDensity,smooth=10))
shared_d = data.frame(apply(shared, 2, makeDensity,smooth=10))
unique_d =data.frame( apply(unique, 2, makeDensity,smooth=10))
P01_d = data.frame(apply(P01, 2, makeDensity,smooth=10))
merged_d = data.frame(apply(merged,2,makeDensity,smooth=10))
merged_Td = data.frame(apply(mergedT,2,makeDensity,smooth=10))

.addAvg(P0_d, "", probs  =probs)

l = list(Benign_all_frag = P01_d, 
         Benign_somatic_frag = merged_d,
         Tumour_all_frag = P0_d,
         Tumour_somatic_frag = merged_Td,
         Tumour_somatic_shared_frag = shared_d,
         Tumour_somatic_uniq_frag = unique_d)

l = lapply(l, .addAvg,"",probs=probs)
.getDens<-function(x, l, nme,id_nme="ID"){
  res = list(x=x)
  for(i in 1:length(l)){
    res[[i+1]] = l[[i]][,which(names(l[[i]]) %in% nme)]
  }
  names(res) = c("x",names(l))
  density=data.frame(res)
  df2 = pivot_longer(density, names_to = id_nme, values_to = nme, cols=names(density)[-1])
  df2
}
density = .getDens(x,l,"X_q_50")
q_25 = .getDens(x,l,"X_q_25", "ID")
q_50 = .getDens(x,l,"X_q_50","ID1")
q_75 = .getDens(x,l,"X_q_75","ID2")

df2 = cbind(q_25, q_50, q_75)[,c(1,2,3,6,9)]

ggp2<-ggplot(df2, aes(x, X_q_50,  fill=ID, color=ID))+
  scale_y_continuous(trans='log10')+
  scale_color_manual(values = c("blue", "black", "green4", "turquoise", "red", "gold4"))+
  ggtitle("Density plot(20bp sliding window)")+
  geom_line()+
  theme_bw()+
  theme(text = element_text(size=textsize))+
  xlab("Fragment length(bp)")+
  ylab("Density") #+
  #geom_ribbon(aes(ymin=X_q_25, ymax=X_q_75, fill = ID),  alpha=0.1)
ggp2

ggsave("density1.png", plot=ggp2, width = 30, height = 30, units = "cm")

