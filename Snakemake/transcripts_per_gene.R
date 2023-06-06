library(ggpubr)

#Read the exons table

exons.tab=read.table('exons_tab.txt',sep='\t',header = T)

#Remove rows that correspond to exons of type I/II that were classified as 0 in other tissues

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Remove rows that repeat the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

pdf("TranscriptsPerGene.pdf",width=7,height=5)

ggboxplot(exons.tab,x='Type',y='TotalTranscripts',fill='Type',palette = 'npg',ylab='#Transcripts',outlier.shape = NA,ylim=c(0,50))+
  theme(legend.position = 'none')

dev.off()