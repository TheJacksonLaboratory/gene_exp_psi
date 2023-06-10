library(data.table)

library(ggplot2)

path.to.gtf='Homo_sapiens.GRCh38.91.gtf'

path.to.gtex.counts='GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz'

path.to.gtex.metadata='GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'

exons.tab=read.table('exons_tab.txt',header = T,sep='\t')

chr=exons.tab$Chr[exons.tab$Type==1 & exons.tab$Tissue=='Brain - Cortex'][1] #exon chromosome

start=exons.tab$Start[exons.tab$Type==1  & exons.tab$Tissue=='Brain - Cortex'][1]   #exon start

end=exons.tab$End[exons.tab$Type==1  & exons.tab$Tissue=='Brain - Cortex'][1]  #exon end

#read the transcripts and exons in order to map the exon to the gene and transcripts that contain it

gtf.file=fread(path.to.gtf,sep='\t',data.table = FALSE)

gtf.file=gtf.file[gtf.file$V3=='exon',]

sub.gtf.file=gtf.file[paste(gtf.file[,1],gtf.file[,4],gtf.file[,5],sep='|')==paste(chr,start,end,sep='|'),]

transcript.with.exon=unique(gsub(';','',unlist(lapply(strsplit(as.character(sub.gtf.file[,9]),split=' '),'[[',6))))

transcript.with.exon=gsub("\"",'',transcript.with.exon)

gene.with.exon=unique(gsub(';','',unlist(lapply(strsplit(as.character(sub.gtf.file[,9]),split=' '),'[[',2))))[1]

gene.with.exon=gsub("\"",'',gene.with.exon)

all.gene.transcripts=gtf.file[grepl(gene.with.exon,gtf.file[,9]),]

all.gene.transcripts=unique(gsub(';','',unlist(lapply(strsplit(as.character(all.gene.transcripts[,9]),split=' '),'[[',6))))

all.gene.transcripts=gsub("\"",'',all.gene.transcripts)

#Next, obtain the GTEx TPM counts for the relevant tissue

gtex.counts=fread(path.to.gtex.counts,data.table = F,nrows=2,skip = 1)

meta.data=read.csv(path.to.gtex.metadata,sep='\t',header=T)

meta.data=meta.data[which(meta.data$SMTSD=='Brain - Cortex'),]

read.cols=c(1,2,which(colnames(gtex.counts) %in% meta.data$SAMPID))

gtex.counts=fread(path.to.gtex.counts,data.table = F,select=read.cols,sep='\t')

meta.data=meta.data[meta.data$SAMPID %in% colnames(gtex.counts),]

meta.data=meta.data[match(meta.data$SAMPID,colnames(gtex.counts)[-c(1,2)]),]

sum(meta.data$SAMPID!=colnames(gtex.counts)[-c(1,2)])  #This should return 0 if the GTEx input is correct.

gtex.counts$transcript_id=gsub('\\.[0-9]*','',gtex.counts$transcript_id)

gtex.counts$gene_id=gsub('\\.[0-9]*','',gtex.counts$gene_id)

#compute the PSI


psi=as.numeric(colSums(gtex.counts[gtex.counts$transcript_id %in% transcript.with.exon,-c(1,2)])/
                 colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)]))     #compute PSI

#compute gene expression
expression=colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)])

#plot expression/PSI
plot(psi,expression,main=paste(gene.with.exon,'-',chr,start,end))

df = data.frame(unlist(psi),unlist(expression))
colnames(df) <- c("PSI", "Expression")
z <- lm(expression ~ psi)

pdf("TypeIexample.pdf",width=7,height=5)

ggplot(df, aes(x=PSI, y=Expression)) +
  geom_point(alpha=0.5, colour="#3C5488", size=2) +
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))  + geom_line(aes(PSI,predict(z)), color = "blue") 

dev.off()



chr=exons.tab$Chr[exons.tab$Type==2  & exons.tab$Tissue=='Brain - Cortex'][1] #exon chromosome

start=exons.tab$Start[exons.tab$Type==2  & exons.tab$Tissue=='Brain - Cortex'][1]   #exon start

end=exons.tab$End[exons.tab$Type==2  & exons.tab$Tissue=='Brain - Cortex'][1]  #exon end


sub.gtf.file=gtf.file[paste(gtf.file[,1],gtf.file[,4],gtf.file[,5],sep='|')==paste(chr,start,end,sep='|'),]

transcript.with.exon=unique(gsub(';','',unlist(lapply(strsplit(as.character(sub.gtf.file[,9]),split=' '),'[[',6))))

transcript.with.exon=gsub("\"",'',transcript.with.exon)

gene.with.exon=unique(gsub(';','',unlist(lapply(strsplit(as.character(sub.gtf.file[,9]),split=' '),'[[',2))))[1]

gene.with.exon=gsub("\"",'',gene.with.exon)

all.gene.transcripts=gtf.file[grepl(gene.with.exon,gtf.file[,9]),]

all.gene.transcripts=unique(gsub(';','',unlist(lapply(strsplit(as.character(all.gene.transcripts[,9]),split=' '),'[[',6))))

all.gene.transcripts=gsub("\"",'',all.gene.transcripts)

#Next, obtain the GTEx TPM counts for the relevant tissue

gtex.counts=fread(path.to.gtex.counts,data.table = F,nrows=2,skip = 1)

meta.data=read.csv(path.to.gtex.metadata,sep='\t',header=T)

meta.data=meta.data[which(meta.data$SMTSD=='Brain - Cortex'),]

read.cols=c(1,2,which(colnames(gtex.counts) %in% meta.data$SAMPID))

gtex.counts=fread(path.to.gtex.counts,data.table = F,select=read.cols,sep='\t')

meta.data=meta.data[meta.data$SAMPID %in% colnames(gtex.counts),]

meta.data=meta.data[match(meta.data$SAMPID,colnames(gtex.counts)[-c(1,2)]),]

sum(meta.data$SAMPID!=colnames(gtex.counts)[-c(1,2)])  #This should return 0 if the GTEx input is correct.

gtex.counts$transcript_id=gsub('\\.[0-9]*','',gtex.counts$transcript_id)

gtex.counts$gene_id=gsub('\\.[0-9]*','',gtex.counts$gene_id)

#compute the PSI


psi=as.numeric(colSums(gtex.counts[gtex.counts$transcript_id %in% transcript.with.exon,-c(1,2)])/
                 colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)]))     #compute PSI

#compute gene expression
expression=colSums(gtex.counts[gtex.counts$transcript_id %in% all.gene.transcripts,-c(1,2)])

#plot expression/PSI
plot(psi,expression,main=paste(gene.with.exon,'-',chr,start,end))

df = data.frame(unlist(psi),unlist(expression))
colnames(df) <- c("PSI", "Expression")
z <- lm(expression ~ psi)

pdf("TypeIIexample.pdf",width=7,height=5)

ggplot(df, aes(x=PSI, y=Expression)) +
  geom_point(alpha=0.5, colour="#3C5488", size=2) +
  theme_bw() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=14))  + geom_line(aes(PSI,predict(z)), color = "blue") 

dev.off()
