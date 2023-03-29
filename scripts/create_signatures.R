#!/usr/bin/env Rscript

#create data structure to store information about exons and the transcripts they belong to
library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two arguments must be supplied (input file).n", call.=FALSE)
} 

gtf_file=args[1]
gtex_file=args[2]


#Get exon names, transcript and gene IDs from the GTF file:
# gtf_file should be the path to Homo_sapiens.GRCh38.91.gtf file
gtf.file=read.table(gtf_file,sep='\t')

gtf.file=gtf.file[gtf.file$V3=='exon',]

exon.names=paste(gtf.file$V1,gtf.file$V4,gtf.file$V5,gtf.file$V7,sep='|')

transcript.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',6)))

transcript.ids=gsub("\"",'',transcript.ids)

gene.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',2)))

gene.ids=gsub("\"",'',gene.ids)

#gtex_file should be the path to GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz
gtex.counts=fread(gtex_file, data.table=F)

#compute a unique signature for exons that belong to the same set of transcripts
exon.signatures=cbind(unique(exon.names),unlist(lapply(unique(exon.names),function(x)paste(unique(transcript.ids[exon.names==x]),collapse='-'))))

#remove exons that only belong to transcript that are not in GTEx
exon.signatures=exon.signatures[unlist(lapply(exon.signatures[,2],function(x)sum(unlist(strsplit(x,split='-')) %in% transcript.ids)>0)),]

#remove exons of genes that only have one transcript
exon.signatures=exon.signatures[unlist(lapply(exon.signatures[,1],function(x)length(unique(transcript.ids[gene.ids %in% gene.ids[exon.names==x]]))>1)),]

#create a table that maps transcripts to exons
transcript.to.exon=cbind(transcript.ids,exon.names,gene.ids)

#count the number of transcripts of each gene
num.gene.transcripts=unlist(lapply(gene.ids,function(x)length(unique(transcript.ids[gene.ids==x]))))

#count the number of transcripts that each exon belongs to
num.exon.transcripts=unlist(lapply(exon.names,function(x)length(unique(transcript.ids[exon.names==x]))))

#count the number of exons of each transcript
num.exons.transcript=unlist(lapply(transcript.ids,function(x)length(unique(exon.names[transcript.ids==x]))))

#remove exons that are in all the gene's transcripts, transcripts that contain one exon and transcripts of genes that only have one transcript
transcript.to.exon=transcript.to.exon[(num.exon.transcripts<num.gene.transcripts) & num.gene.transcripts>1  & num.exons.transcript>1,]

#keep only one exon for each signature
transcript.to.exon=transcript.to.exon[transcript.to.exon[,2] %in% exon.signatures[!duplicated(exon.signatures[,2]),1],]

save.image('after_exon_sig_next.RData')

