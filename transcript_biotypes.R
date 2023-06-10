library(biomaRt)

library(ggpubr)

#Read the exons table

exons.tab=read.table('exons_tab.txt',sep='\t',header = T)

#Remove rows that correspond to exons of type I/II that were classified as 0 in other tissues

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Remove rows that repeat the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Get all the transcripts in the exons table

exon.transcripts=lapply(1:nrow(exons.tab),function(i)unlist(strsplit(exons.tab$TranscriptsContain[i],split=',')))

#Get transcript biotypes from biomart

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host="https://www.ensembl.org")

transcript.types=getBM(attributes = c("ensembl_transcript_id","transcript_biotype"),
                       filters    = "ensembl_transcript_id",
                       values     = unique(unlist(exon.transcripts)), 
                       mart       = mart)

#Prepare arrays with transcripts that contain each type of exons

type.1.transcripts=unlist(lapply(exons.tab$TranscriptsContain[exons.tab$Type==1],function(x)unlist(strsplit(x,','))))

type.2.transcripts=unlist(lapply(exons.tab$TranscriptsContain[exons.tab$Type==2],function(x)unlist(strsplit(x,','))))

type.0.transcripts=unlist(lapply(exons.tab$TranscriptsContain[exons.tab$Type==0],function(x)unlist(strsplit(x,','))))

#Prepare arrays with the biotypes of the transcripts in each array:

tr.1.types=transcript.types$transcript_biotype[transcript.types$ensembl_transcript_id %in% type.1.transcripts]

tr.2.types=transcript.types$transcript_biotype[transcript.types$ensembl_transcript_id %in% type.2.transcripts]

tr.0.types=transcript.types$transcript_biotype[transcript.types$ensembl_transcript_id %in% type.0.transcripts]

#Create a table of biotypes for each type of transcripts, divided by the total number of biotypes for that transcript type

type.counts=c(table(tr.0.types)/length(tr.0.types),table(tr.1.types)/length(tr.1.types),table(tr.2.types)/length(tr.2.types))

#Create a data frame and plot

df=data.frame(TranscriptType=names(type.counts),FracTranscripts=type.counts,
              ExonType=c(rep(0,length(unique(tr.0.types))),rep(1,length(unique(tr.1.types))),rep(2,length(unique(tr.2.types)))))

pdf("TranscriptBiotypes.pdf",width=7,height=5)

ggbarplot(df,x='ExonType',y='FracTranscripts',fill = 'TranscriptType')+ theme(legend.title = element_blank()) 

dev.off()

