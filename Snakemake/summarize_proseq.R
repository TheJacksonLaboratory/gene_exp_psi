library(ggpubr)

library(rstatix)

#Read the bed file containing the coordinates of the downstream introns of every exon 

exons.coord=read.table('exons_coord.bed',header=F,sep='\t')

#Read the bed files with the containing the intersections of intron coordinates and read 3' ends

bedfiles=list.files('.',pattern = 'exons_proseq_*')

file2=read.table(bedfiles[2])

file2=file2[file2$V6==file2$V13,]  #make sure that read and intron directionalities match

file1=read.table(bedfiles[1])

file1=file1[file1$V6==file1$V13,] #make sure that read and intron directionalities match

#For each gene, summarize the reads that fall into the downstream introns of each of its exon types, including the two biological samples

res.tab=do.call(rbind,lapply(unique(exons.coord$V7),function(x)
{

  intersect.file=file1
    
  intersect.file=intersect.file[intersect.file$V7==x,]  #Get the reads that intersect with gene x

  #Sum the lengths (in base pairs) of downstream introns of type 0,I, and II exons that belong to gene x
  
  len.type.0=sum(exons.coord$V3[exons.coord$V7==x & grepl('\\_0$',exons.coord$V4)]-exons.coord$V2[exons.coord$V7==x & grepl('\\_0$',exons.coord$V4)]+1)
  
  len.type.1=sum(exons.coord$V3[exons.coord$V7==x & grepl('\\_1$',exons.coord$V4)]-exons.coord$V2[exons.coord$V7==x & grepl('\\_1$',exons.coord$V4)]+1)
  
  len.type.2=sum(exons.coord$V3[exons.coord$V7==x & grepl('\\_2$',exons.coord$V4)]-exons.coord$V2[exons.coord$V7==x & grepl('\\_2$',exons.coord$V4)]+1)
  
  len.exp=sum(exons.coord$V5[exons.coord$V7==x]) #This is the length of all the introns together
  
  #Count the number of reads that fall into downstream introns of the 3 exon types, in the control sample. Divide this number by the bp length. 
  #If the gene doesn't have that type (e.g. it only has type 0 and I but not type II) set the value to NA
  
  num.type.2.control=ifelse(len.type.2>0,sum(grepl('\\_2$',intersect.file$V4))/len.type.2,NA)
  
  num.type.0.control=ifelse(len.type.0>0,sum(grepl('\\_0$',intersect.file$V4))/len.type.0,NA)
  
  num.type.1.control=ifelse(len.type.1>0,sum(grepl('\\_1$',intersect.file$V4))/len.type.1,NA)
  
  exp.control=length(unique(intersect.file$V11))/len.exp  #The number of reads that intersect with the gene in total
  
  #Compute the same values for the second biological sample (case)
  
  intersect.file=file2
  
  intersect.file=intersect.file[intersect.file$V7==x,]
  
  num.type.2.case=ifelse(len.type.2>0,sum(grepl('\\_2$',intersect.file$V4))/len.type.2,NA)
  
  num.type.0.case=ifelse(len.type.0>0,sum(grepl('\\_0$',intersect.file$V4))/len.type.0,NA)
  
  num.type.1.case=ifelse(len.type.1>0,sum(grepl('\\_1$',intersect.file$V4))/len.type.1,NA)
  
  exp.case=length(unique(intersect.file$V11))/len.exp

  #Return the values for gene x, which will be appended to the result
  
  res=matrix(c(num.type.1.case,exp.case,'1',
               num.type.2.case,exp.case,'2',
               num.type.0.case,exp.case,'0',
               num.type.1.control,exp.control,'1',
               num.type.2.control,exp.control,'2',
               num.type.0.control,exp.control,'0'),ncol=3,byrow = T)
      
  return(res)
  
}))


#Create logical variables that select non-zero values for each type of exon (because if an exon is not expressed we cannot estimate its speed of transcription)

samp.1=unlist(lapply(1:nrow(res.tab),function(i)res.tab[i,3]=='1' & ifelse(!is.na(res.tab[i,1]),as.numeric(res.tab[i,1])>0,F)))

samp.2=unlist(lapply(1:nrow(res.tab),function(i)res.tab[i,3]=='2' & ifelse(!is.na(res.tab[i,1]),as.numeric(res.tab[i,1])>0,F)))

samp.0=unlist(lapply(1:nrow(res.tab),function(i)res.tab[i,3]=='0' & ifelse(!is.na(res.tab[i,1]),as.numeric(res.tab[i,1])>0,F)))

#Create a data frame with the data, compute p-values and plot

df=data.frame(type=c(rep('1',sum(samp.1)),rep('2',sum(samp.2)),rep('0',sum(samp.0))),
              per.bp.downstream.intron=as.numeric(c(res.tab[samp.1,1],res.tab[samp.2,1],res.tab[samp.0,1]))/as.numeric(c(res.tab[samp.1,2],res.tab[samp.2,2],res.tab[samp.0,2])))
    
stat.test <-( df %>%
                group_by('type') %>%
                wilcox_test(per.bp.downstream.intron ~ type,paired = F) %>%
                add_significance())
stat.test

pdf("WissinkOverlap.pdf",width=7,height=5)

ggboxplot(df,x='type',y='per.bp.downstream.intron',fill='type',xlab='Type',ylab='Reads per BP (intron/gene)',
        palette = 'npg',outlier.shape = NA,ylim=c(0,3))+ theme(legend.position = 'none')+ stat_pvalue_manual(
          stat.test,  label = "p", tip.length = 0.001,step.increase =rep(0.01,3) ,y.position = 2,label.size = 4,bracket.shorten = 0.1
        )

dev.off()

wilcox.test(df$per.bp.downstream.intron[df$type=='2'],df$per.bp.downstream.intron[df$type=='0'],paired = F)

wilcox.test(df$per.bp.downstream.intron[df$type=='1'],df$per.bp.downstream.intron[df$type=='0'],paired = F)



