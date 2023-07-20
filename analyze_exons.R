#combine the results from the different tissues into one table of T1/T2/T0 exons

library(intervals)

library(seqinr)

#create a matrix that will contain the information about the exons:

exons.tab=matrix(nrow=0,ncol=18) #***

colnames(exons.tab)=c('Name','Type','Chr','Start','End','UpIntronStart','UpIntronEnd','DownIntronStart','DownIntronEnd','Coef','Slope','ExonNumber',
                      'NumExons','TranscriptsContain','TotalTranscripts','GeneId','Tissue','CD1cor') #***

tissues=c('Spleen','Thyroid','Brain - Cortex','Adrenal Gland','Breast - Mammary Tissue','Heart - Left Ventricle',
          'Liver','Pituitary','Pancreas') #Lung

load('after_exon_sig_next.RData')

for (tissue in tissues)
{
  load(paste0('types_',tissue,'.RData'))  #load the results for the relevant tissue

  cyclin.d1=colSums(gtex.counts[gtex.counts$gene_id=='ENSG00000110092',-c(1,2)]) #***
  
  rm(gtex.counts)  #free memory 

  rm(transcript.to.exon)

  adj.p=p.adjust(cor.vals[,2],method='BH')  #transform the p-values into BH scores
  
  for (exon.num in which(adj.p<=0.05 & cor.vals[,3]>=0.5 & rowMeans(exon.gene.counts)>=20))  #consider exons with low BH, high R^2 and above a threshold of gene expression
  {
    exon.name=unlist(strsplit(rownames(exon.gene.counts)[exon.num],split=' '))[[2]]
    
    exon.signature=exon.signatures[exon.signatures[,1]==exon.name,2]
    
    equivalent.exons=exon.signatures[exon.signatures[,2]==exon.signature,1]
    
    for (next.exon in equivalent.exons)  #The same will apply to all the exons with the same signature
    {
      
      exon.type=ifelse(cor.vals[exon.num,1]>0,1,2)
       
      exon.length=unlist(strsplit(next.exon,'\\|'))  #Get the location by splitting the exon name
      
      exon.start=as.integer(exon.length[2])
      
      exon.end=as.integer(exon.length[3])
     
      exon.chr=exon.length[1]
      
      exon.dir=exon.length[4]
      
      exon.length=abs(exon.end-exon.start+1)
      
      exon.gene=unlist(strsplit(rownames(exon.gene.counts)[exon.num],split=' '))[1]
      
      exon.tr=transcript.ids[exon.names==next.exon][1]
      
      gtf.part=gtf.file[gene.ids==exon.gene,]  #relevant GTF rows of the exon's gene
      
      merged.exons=interval_union(Intervals(gtf.part[,4:5]))   #To compute introns, merge all overlapping exons
    
      exon.row=as.integer(interval_overlap(Intervals(c(exon.start,exon.end)),merged.exons))  #Find the location of the exon in the merged GTF
      
      #get the introm coordinates according to the directionality of the gene
      
      if (exon.dir=='+')
      {
        
        down.intron.start=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row,2]+1)
        
        down.intron.end=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row+1,1]-1)
        
        up.intron.start=ifelse(exon.row==1,0,merged.exons[exon.row-1,2]+1)
        
        up.intron.end=ifelse(exon.row==1,0,merged.exons[exon.row,1]-1)
        
        exon.loc=exon.row
        
      }else{
        
        up.intron.start=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row+1,1]-1)
        
        up.intron.end=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row,2]+1)
        
        down.intron.start=ifelse(exon.row==1,0,merged.exons[exon.row-1,2]+1)
        
        down.intron.end=ifelse(exon.row==1,0,merged.exons[exon.row,1]-1)
        
        exon.loc=(nrow(gtf.part)-exon.row+1)
        
      }
      
      #fit a linear model to get the relevant stats for the exon
      fit.mod=summary(lm(exon.gene.counts[exon.num,]~exon.proportions[exon.num,],na.action='na.omit'))
      
      cd1.cor=cor(cyclin.d1,exon.proportions[exon.num,]) #***
      
      #add all the exon info
      exons.tab=rbind(exons.tab,
                      c(next.exon,exon.type,exon.chr,exon.start,exon.end,up.intron.start,up.intron.end,down.intron.start,down.intron.end,
                        fit.mod$coefficients[1,1],fit.mod$coefficients[2,1],exon.loc,nrow(merged.exons),
                        paste(unique(transcript.ids[exon.names==next.exon & gene.ids==exon.gene]),collapse = ','), length(unique(transcript.ids[gene.ids==exon.gene]))
                        ,exon.gene,tissue,cd1.cor))#***
      
    }
    
  }
  
  #Retrieve similar information for type 0 exons, i.e. those with a high BH score, above minimal gene expression and any R^2
  
  for (exon.num in which(adj.p>0.5 & rowMeans(exon.gene.counts)>=20 & !is.na(cor.vals[,3])))
  {
    
    exon.name=unlist(strsplit(rownames(exon.gene.counts)[exon.num],split=' '))[[2]]
    
    exon.signature=exon.signatures[exon.signatures[,1]==exon.name,2]
    
    equivalent.exons=exon.signatures[exon.signatures[,2]==exon.signature,1]
    
    for (next.exon in equivalent.exons)
    {
      
        exon.type=0
        
        exon.length=unlist(strsplit(next.exon,'\\|'))
        
        exon.start=as.integer(exon.length[2])
        
        exon.end=as.integer(exon.length[3])
        
        exon.chr=exon.length[1]
        
        exon.dir=exon.length[4]
        
        exon.length=abs(exon.end-exon.start+1)
        
        exon.tr=transcript.ids[exon.names==next.exon][1]
        
        exon.gene=unlist(strsplit(rownames(exon.gene.counts)[exon.num],split=' '))[1]
        
        gtf.part=gtf.file[gene.ids==exon.gene,]
        
        merged.exons=interval_union(Intervals(gtf.part[,4:5]))
        
        exon.row=as.integer(interval_overlap(Intervals(c(exon.start,exon.end)),merged.exons))
        
        if (exon.dir=='+')
        {
          
          down.intron.start=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row,2]+1)
          
          down.intron.end=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row+1,1]-1)
          
          up.intron.start=ifelse(exon.row==1,0,merged.exons[exon.row-1,2]+1)
          
          up.intron.end=ifelse(exon.row==1,0,merged.exons[exon.row,1]-1)
          
          exon.loc=exon.row
          
        }else{
          
          up.intron.start=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row+1,1]-1)
          
          up.intron.end=ifelse(exon.row==nrow(merged.exons),0,merged.exons[exon.row,2]+1)
          
          down.intron.start=ifelse(exon.row==1,0,merged.exons[exon.row-1,2]+1)
          
          down.intron.end=ifelse(exon.row==1,0,merged.exons[exon.row,1]-1)
          
          exon.loc=(nrow(gtf.part)-exon.row+1)
          
        }
        
       
        
        exons.tab=rbind(exons.tab,
                        c(next.exon,exon.type,exon.chr,exon.start,exon.end,up.intron.start,up.intron.end,down.intron.start,down.intron.end,
                          0,0,exon.loc,nrow(merged.exons),
                          paste(unique(transcript.ids[exon.names==next.exon & gene.ids==exon.gene]),collapse = ','), length(unique(transcript.ids[gene.ids==exon.gene]))
                          ,exon.gene,tissue))
    }
  }
}


write.table(exons.tab,'exons_tab.txt',sep='\t',col.names=T,row.names=F,quote=F)


