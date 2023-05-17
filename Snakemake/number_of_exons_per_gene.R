library(ggpubr)

#Read the exons table

exons.tab=read.table('exons_tab.txt',sep='\t',header = T)

#Remove rows that correspond to exons of type I/II that were classified as 0 in other tissues

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Remove rows that repeat the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Create arrays with the gene IDs of exons that contain type I, type II and type 0 exons

type.1.genes=unique(exons.tab$GeneId[exons.tab$Type==1])

type.2.genes=unique(exons.tab$GeneId[exons.tab$Type==2])

type.0.genes=unique(exons.tab$GeneId[exons.tab$Type==0])

#create a data frame with the number of exons per gene for each type of genes

df=data.frame(Type=c(rep(0,sum((exons.tab$GeneId %in% type.0.genes) & !(duplicated(exons.tab$GeneId)))),
                     rep(1,sum((exons.tab$GeneId %in% type.1.genes) & !(duplicated(exons.tab$GeneId)))),
                     rep(2,sum((exons.tab$GeneId %in% type.2.genes) & !(duplicated(exons.tab$GeneId))))),
              NumExons=c(exons.tab$NumExons[(exons.tab$GeneId %in% type.0.genes) & !(duplicated(exons.tab$GeneId))],
                         exons.tab$NumExons[(exons.tab$GeneId %in% type.1.genes) & !(duplicated(exons.tab$GeneId))],
                         exons.tab$NumExons[(exons.tab$GeneId %in% type.2.genes) & !(duplicated(exons.tab$GeneId))]))

pdf("ExonsPerGene.pdf",width=7,height=5)

ggboxplot(df,x='Type',y='NumExons',fill='Type',palette = 'npg',ylab='#Exons',outlier.shape = NA,ylim=c(0,30))+theme(legend.position = 'none')

dev.off()


