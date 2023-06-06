library(ggpubr)

#Read the exons table

exons.tab=read.table('exons_tab.txt',sep='\t',header = T)

#Remove rows that correspond to exons of type I/II that were classified as 0 in other tissues

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Remove rows that repeat the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Compute exons lengths and create data frame

exons.tab$Length=exons.tab$End-exons.tab$Start+1

df=data.frame(Length=c(exons.tab$Length[exons.tab$Type==2],
                       exons.tab$Length[exons.tab$Type==1],
                       exons.tab$Length[exons.tab$Type==0]),
              Type=c(rep(2,sum(exons.tab$Type==2)),
                     rep(1,sum(exons.tab$Type==1)),
                     rep(0,sum(exons.tab$Type==0))))


pdf("ExonsLengths.pdf",width=7,height=5)

ggboxplot(df,x='Type',y='Length',fill='Type',palette = 'npg',ylab='Exon Length (bp)',outlier.shape = NA,ylim=c(0,2200))+
  theme(legend.position = 'none')

dev.off()

