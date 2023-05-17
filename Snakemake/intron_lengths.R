library(ggpubr)

#Read the exons table

exons.tab=read.table('exons_tab.txt',sep='\t',header = T)

#Remove rows that correspond to exons of type I/II that were classified as 0 in other tissues

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Remove rows that repeat the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Create an array of exon directionality:

dir= unlist(lapply(exons.tab$Name,function(x)strsplit(x,split='\\|')[[1]][4]))

#Compute intron lengths:

exons.tab$UpIntronLength=ifelse(dir=='+',exons.tab$UpIntronEnd-exons.tab$UpIntronStart+1,exons.tab$UpIntronStart-exons.tab$UpIntronEnd+1)

exons.tab$DownIntronLength=ifelse(dir=='+',exons.tab$DownIntronEnd-exons.tab$DownIntronStart+1,exons.tab$DownIntronEnd-exons.tab$DownIntronStart+1)

#Create a data frame and plot downstream intron lengths:

pdf("DownstreamIntronLengths.pdf",width=7,height=5)

df=data.frame(Length=c(exons.tab$DownIntronLength[exons.tab$Type==2 & exons.tab$DownIntronLength>0],
                       exons.tab$DownIntronLength[exons.tab$Type==1 & exons.tab$DownIntronLength>0],
                       exons.tab$DownIntronLength[exons.tab$Type==0 & exons.tab$DownIntronLength>0]),
              Type=c(rep(2,sum(exons.tab$Type==2 & exons.tab$DownIntronLength>0)),
                     rep(1,sum(exons.tab$Type==1 & exons.tab$DownIntronLength>0)),
                     rep(0,sum(exons.tab$Type==0 & exons.tab$DownIntronLength>0))))

ggboxplot(df,x='Type',y='Length',fill='Type',palette = 'npg',ylab='Downstream Intron Length (bp)',outlier.shape = NA,
          ylim=c(0,7000))+
  theme(legend.position = 'none')+geom_hline(yintercept=median(df$Length[df$Type==0]), linetype='dotted', col = 'red')

dev.off()

pdf("UpstreamIntronLengths.pdf",width=7,height=5)

df=data.frame(Length=c(exons.tab$UpIntronLength[exons.tab$Type==2 & exons.tab$UpIntronLength>0],
                       exons.tab$UpIntronLength[exons.tab$Type==1 & exons.tab$UpIntronLength>0],
                       exons.tab$UpIntronLength[exons.tab$Type==0 & exons.tab$UpIntronLength>0]),
              Type=c(rep(2,sum(exons.tab$Type==2 & exons.tab$UpIntronLength>0)),
                     rep(1,sum(exons.tab$Type==1 & exons.tab$UpIntronLength>0)),
                     rep(0,sum(exons.tab$Type==0 & exons.tab$UpIntronLength>0))))

ggboxplot(df,x='Type',y='Length',fill='Type',palette = 'npg',ylab='Upstream Intron Length (bp)',outlier.shape = NA,
          ylim=c(0,6500))+
  theme(legend.position = 'none')+geom_hline(yintercept=median(df$Length[df$Type==0]), linetype='dotted', col = 'red')

dev.off()


