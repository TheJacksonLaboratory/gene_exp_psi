#create data structure to store information about exons and the transcripts they belong to

library(data.table)

library(parallel)

num.cores=30

#Get exon names, transcript and gene IDs from the GTF file:

gtf.file=read.table('Homo_sapiens.GRCh38.91.gtf',sep='\t')

gtf.file=gtf.file[gtf.file$V3=='exon',]

exon.names=paste(gtf.file$V1,gtf.file$V4,gtf.file$V5,gtf.file$V7,sep='|')

transcript.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',6)))

transcript.ids=gsub("\"",'',transcript.ids)

gene.ids=gsub(';','',unlist(lapply(strsplit(as.character(gtf.file[,9]),split=' '),'[[',2)))

gene.ids=gsub("\"",'',gene.ids)

gtex.counts=fread('GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz',data.table = F)

#compute a unique signature for exons that belong to the same set of transcripts

exon.signatures=cbind(unique(exon.names),unlist(mclapply(unique(exon.names),function(x)paste(unique(transcript.ids[exon.names==x]),collapse='-'),mc.cores = num.cores)))

#remove exons of genes that only have one transcript

exon.signatures=exon.signatures[unlist(mclapply(exon.signatures[,1],function(x)length(unique(transcript.ids[gene.ids %in% gene.ids[exon.names==x]]))>1,mc.cores = num.cores)),]

#create a table that maps transcripts to exons

transcript.to.exon=cbind(transcript.ids,exon.names,gene.ids)

#count the number of transcripts of each gene

num.gene.transcripts=unlist(mclapply(gene.ids,function(x)length(unique(transcript.ids[gene.ids==x])),mc.cores = num.cores))

#count the number of transcripts that each exon belongs to

num.exon.transcripts=unlist(mclapply(exon.names,function(x)length(unique(transcript.ids[exon.names==x])),mc.cores = num.cores))

#count the number of exons of each transcript

num.exons.transcript=unlist(mclapply(transcript.ids,function(x)length(unique(exon.names[transcript.ids==x])),mc.cores = num.cores))

#remove exons that are in all the gene's transcripts, transcripts that contain one exon and transcripts of genes that only have one transcript

transcript.to.exon=transcript.to.exon[(num.exon.transcripts<num.gene.transcripts) & num.gene.transcripts>1  & num.exons.transcript>1,]

#keep only one exon for each signature

transcript.to.exon=transcript.to.exon[transcript.to.exon[,2] %in% exon.signatures[!duplicated(exon.signatures[,2]),1],]

save.image('after_exon_sig_next.RData')

