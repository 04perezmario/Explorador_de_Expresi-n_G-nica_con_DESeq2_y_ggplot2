library(readr)
library(ggpubr)
library(tidyr)
counts<-read.table("/home/labs/lslab/studentlslab/out/petits/STAR_mapped/readcounts.txt",sep="\t")
DESeq.out<-read.table("/home/labs/lslab/studentlslab/out/petits/STAR_mapped/DSEQ_Results.txt",sep="\t")

counts.long<-pivot_longer(counts["MYC",],cols=c(colnames(counts["MYC",])),names_to=c("sample"),values_to=c("raw"))
View(counts.long)
counts.long$group<-counts.long$sample

for (i in 1:length(counts.long$sample)){
counts.long$group[i]<-unlist(strsplit(counts.long$sample,"_")[[i]])[1]
}
counts.long$group<-as.factor(counts.long$group)

p<-ggboxplot(as.data.frame(counts.long),x="group",y="raw")
p

stats.df<-data.frame(.y.=  
                     group1=
                     group2        
                     p   
                     p.adj 
                     p.format 
                     p.signif 
                     method)
p + stat_pvalue_manual(, label = "p = {p.adj}")
