library(data.table)

library(biomaRt)

#read isoform-level annotations 

isopret.tab=read.table('isoform_function_list_CC.txt',sep='\t',header=T)

isopret.tab=rbind(isopret.tab,read.table('isoform_function_list.txt',sep='\t',header=T))

isopret.tab=rbind(isopret.tab,read.table('isoform_function_list_BP.txt',sep='\t',header=T))

isopret.tab=unique(isopret.tab)

#read the exons table

exons.tab=fread('exons_tab.txt',sep='\t',header=T)

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#create sets of UHP/DHP type 0 transcripts

type.1.transcripts=unique(unlist(lapply(
  exons.tab$TranscriptsContain[exons.tab$Type==1],strsplit,split=',')))

type.2.transcripts=unique(unlist(lapply(
  exons.tab$TranscriptsContain[exons.tab$Type==2],strsplit,split=',')))

type.0.transcripts=unique(setdiff(unlist(lapply(
  exons.tab$TranscriptsContain,strsplit,split=',')),c(type.1.transcripts,type.2.transcripts)))

type.1.transcripts=type.1.transcripts[type.1.transcripts %in% isopret.tab$Ensembl.ID]

type.2.transcripts=type.2.transcripts[type.2.transcripts %in% isopret.tab$Ensembl.ID]

type.0.transcripts=type.0.transcripts[type.0.transcripts %in% isopret.tab$Ensembl.ID]

isopret.tab=isopret.tab[isopret.tab$Ensembl.ID %in% c(type.1.transcripts,type.2.transcripts,type.0.transcripts),]

#conduct hypergeometric test for each GO term
enrich=do.call(rbind,lapply(unique(isopret.tab$Go.Terms),function(term){
  
  m=sum(type.1.transcripts %in% isopret.tab$Ensembl.ID)
  
  k = sum(isopret.tab$Go.Terms==term)
  
  n = sum(!unique(isopret.tab$Ensembl.ID) %in% type.1.transcripts)
  
  q=sum((isopret.tab$Go.Terms==term) & (isopret.tab$Ensembl.ID %in% type.1.transcripts))
  
  if (q<5)
    
    return(NULL)
  
  pval=phyper(m=m,
              k=k ,
              n=n ,
              q=q-1,lower.tail=F)
  
  
  return(c(term,pval,q/k))
  
}))

colnames(enrich)=c('term','p.value','Coverage')
#Bonferroni correction
enrich=cbind(enrich,p.adjust(enrich[,2]))

colnames(enrich)[ncol(enrich)]='adj.p'

enrich=enrich[order(as.numeric(enrich[,4])),]
#download GO term names from Biomart
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="https://www.ensembl.org")

enrich=enrich[as.numeric(enrich[,4])<=0.05,]

enrich=cbind(enrich,rep('',nrow(enrich)))

for (i in 1:nrow(enrich))
{
  term=getBM(attributes = c("go_id","name_1006"),
        filters    = "go",
        values     = enrich[i,1], 
        mart       = mart)
  
  enrich[i,5]=term$name_1006[term$go_id==enrich[i,1]]

  
}
