library(ggpubr)

#read the exon coordinates and intersection with ChIP-Seq files

exons.coord=read.table('exons_coord_chip.bed',header=F,sep='\t')

bedfiles=list.files('./',pattern = 'exons_chip_intersection_*')

#For each gene, sum the number of binding in the data and divide by the length of the respective exons

res.tab=do.call(rbind,lapply(unique(exons.coord$V7),function(x)
{
  res=matrix(ncol=4,nrow=0)
  
  #get exon lengths
  
  len.type.0=sum(exons.coord$V5[exons.coord$V7==x & grepl('\\_0$',exons.coord$V4)])
  
  len.type.1=sum(exons.coord$V5[exons.coord$V7==x & grepl('\\_1$',exons.coord$V4)])
  
  len.type.2=sum(exons.coord$V5[exons.coord$V7==x & grepl('\\_2$',exons.coord$V4)])
  
  if ((len.type.1==0 && len.type.2==0) || len.type.0==0)
    
    return(NULL)
  
  for (bedfile in bedfiles)
  {
    if (file.size(bedfile)==0)
      
      next
    
    intersect.file=read.table(bedfile)
    
    intersect.file=intersect.file[intersect.file$V7==x,]
    
    if (nrow(intersect.file)==0)
      
      next
    
    #Get number of bindings
    
    num.type.2=sum(grepl('\\_2$',intersect.file$V4))
    
    num.type.0=sum(grepl('\\_0$',intersect.file$V4))
    
    num.type.1=sum(grepl('\\_1$',intersect.file$V4))
    
    #add the next binding count to the results
    
    res=rbind(res,c(x,num.type.0/max(len.type.0,1),num.type.1/max(len.type.1,1),num.type.2/max(len.type.2,1)))
  }    
  
  if (nrow(res)==0)
    
    return(NULL)
  
  ret=res
  
  #Sum the total number of bindings for this gene
  
  if (nrow(res)>1)
    
    ret=c(res[1,1],apply(res[,2:4],2,function(v)sum(as.numeric(v))))
  
  #Compute the relative type I/II vs type 0 bindings
  
  if (len.type.1==0)
  {
    ret[3]=NA
  }else{
    ret[3]=as.numeric(ret[3])/as.numeric(ret[2])
  }
  
  if (len.type.2==0)
  {
    ret[4]=NA
  }else{
    
    ret[4]=as.numeric(ret[4])/as.numeric(ret[2])
  }
  
  return(ret[c(1,3,4)])  #return the values for this gene
  
}))

#write to file and plot the results

write.table(res.tab,'chip_intersection_summary.txt',col.names = T,row.names = F,quote = F)  

res.tab=res.tab[!(res.tab[,2]=='Inf') & !(res.tab[,3]=='Inf'),]

df=data.frame(Type=c(rep('1 vs. 0',nrow(res.tab)),rep('2 vs. 0',nrow(res.tab))),Bindings=c(as.numeric(res.tab[,2]),
                                                                                           as.numeric(res.tab[,3])))

pdf("ChIPSeqPlot.pdf",width=7,height=5)

ggboxplot(df,x='Type',y='Bindings',fill='Type',palette = 'npg',ylab='Bindings per BP Ratio',ylim=c(0,7),outlier.shape = NA)+
  theme(legend.position = 'none')+geom_hline(yintercept=1, linetype='dotted', col = 'red')

dev.off()

