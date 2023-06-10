#read all gtex
#convert from transcript counts to exon counts
#create another matrix with proportions
#compute correlation between rows of the two matrices

library(data.table)

library(parallel)

num.cores=30

tissues=c('Lung','Spleen','Thyroid','Brain - Cortex','Adrenal Gland','Breast - Mammary Tissue','Heart - Left Ventricle',
          'Liver','Pituitary','Pancreas')

load('after_exon_sig_next.RData')

rm(exon.signatures)

rm(num.gene.transcripts)

rm(num.exon.transcripts)

args=commandArgs(trailingOnly = TRUE) 

tissue.index=as.integer(args[1])

gtf.file=read.table('Homo_sapiens.GRCh38.91.gtf',sep='\t')

gtf.file=gtf.file[gtf.file$V3=='exon',]

exon.names=paste(gtf.file$V1,gtf.file$V4,gtf.file$V5,gtf.file$V7,sep='|')

transcript.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',6)))

transcript.ids=gsub("\"",'',transcript.ids)

gene.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',2)))

gene.ids=gsub("\"",'',gene.ids)

gtex.counts=fread('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz',data.table = F)

meta.data=read.csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',sep='\t',header=T)

meta.data=meta.data[which(meta.data$SMTSD==tissues[tissue.index]),]

gtex.counts=gtex.counts[,c(1,2,which(colnames(gtex.counts) %in% meta.data$SAMPID))]

meta.data=meta.data[meta.data$SAMPID %in% colnames(gtex.counts),]

meta.data=meta.data[match(meta.data$SAMPID,colnames(gtex.counts)[-c(1,2)]),]

sum(meta.data$SAMPID!=colnames(gtex.counts)[-c(1,2)])

gtex.counts$transcript_id=gsub('\\.[0-9]*','',gtex.counts$transcript_id)

gtex.counts$gene_id=gsub('\\.[0-9]*','',gtex.counts$gene_id)


print(tissues[tissue.index])


exon.proportions=do.call(rbind,lapply(unique(transcript.to.exon[,2]),function(x){
  
  tte.row=which(transcript.to.exon[,2]==x)

  genes=transcript.to.exon[tte.row,3]
  
  res=do.call(rbind,lapply(unique(genes),function(x){
    
    in.transcripts=transcript.to.exon[tte.row[genes==x],1]
    
    all.gene.transcripts=gtex.counts$transcript_id[gtex.counts$gene_id==x]
    
    colSums(gtex.counts[gtex.counts$transcript_id %in% in.transcripts,-c(1,2)])/
     colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)])
  }))
  
  res[is.nan(res) | is.infinite(res)]=NA
  
  rownames(res)=paste(unique(genes),x)
  
  res
  
}))

print('exon.proportions created')



exon.gene.counts=do.call(rbind,lapply(unique(transcript.to.exon[,2]),function(x){
  
  tte.row=which(transcript.to.exon[,2]==x)
  
  genes=transcript.to.exon[tte.row,3]
  
  res=do.call(rbind,lapply(unique(genes),function(x){
    
    all.gene.transcripts=gtex.counts$transcript_id[gtex.counts$gene_id==x]
    
    colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)])
  }))
  
  rownames(res)=paste(unique(genes),x)
  
  res[is.nan(res) | is.infinite(res)]=NA
  
  res
  
}))



print('exon.gene.counts created')

cor.vals=do.call(rbind,lapply(1:nrow(exon.proportions),function(i){
  
  if(sum(is.na(exon.proportions[i,]))>=length(exon.proportions[i,])/2)
    return(c(NA,NA,NA))
  
  if (quantile(exon.gene.counts[i,],0.95)/quantile(exon.gene.counts[i,],0.05)<2)
    return(c(NA,NA,NA))
  
  if (sd(exon.proportions[i,],na.rm=T)==0 || sd(exon.gene.counts[i,],na.rm=T)==0)
    return(c(NA,NA,NA))
  
  res=summary(lm(exon.gene.counts[i,]~exon.proportions[i,],na.action='na.omit'))
  
  c(res$coefficients[2,1],res$coefficients[2,4],res$r.squared)
}))

colnames(cor.vals)=c('coef','pval','r.squared')

save.image(paste0('types_',tissues[tissue.index],'.RData'))




