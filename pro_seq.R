#This files reads the PRO-Seq bam files that were downloaded from ENCODE (Wissink et al.), and intersects them with the downstream introns of the exons
#in the exons_tab.txt table.  Wissink et al. has two biological samples, cases and controls.

library(BRGenomics)

library(data.table)

#import and combine the bam files that correspond to the first biological sample
ps <- import_bam(file=c('ENCFF098OEZ.bam','ENCFF672CXH.bam','ENCFF251WTK.bam','ENCFF857RTX.bam'),
                 mapq = 20, 
                 revcomp = TRUE,
                 shift = -1,
                trim.to = "3p",
                 paired_end = FALSE)

ps=do.call("c", ps)

#Write the coordinates of the 3' end of each read to a bed file
write.table(cbind(as.character(ps@seqnames),data.frame(ps@ranges)[,1:2],paste0('pro_seq_',1:length(ps)),data.frame(ps@ranges)[,3],strand(ps)),
            paste0('proseq_1.bed'),sep='\t',row.names = F,col.names = F,quote = F)

#import and combine the bam files that correspond to the second biological sample
ps <- import_bam(file=c('ENCFF944ENF.bam','ENCFF480XBL.bam','ENCFF710QMM.bam','ENCFF933PEI.bam'),
                 mapq = 20, 
                 revcomp = TRUE,
                 shift = -1,
                 trim.to = "3p",
                 paired_end = FALSE)

ps=do.call("c", ps) 

#Write the coordinates of the 3' end of each read to a bed file
write.table(cbind(as.character(ps@seqnames),data.frame(ps@ranges)[,1:2],paste0('pro_seq_',1:length(ps)),data.frame(ps@ranges)[,3],strand(ps)),
            paste0('proseq_2.bed'),sep='\t',row.names = F,col.names = F,quote = F)

#read the exons tab
exons.tab=fread('../exons_tab.txt',sep='\t',header=T)

#Remove rows that correspond to tissues in which type I/II exons were classified as type 0
exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Keep only genes that have at least one type I or type II exon
exons.tab=exons.tab[exons.tab$GeneId %in% exons.tab$GeneId[exons.tab$Type!=0],]

#Remove rows that correspond to the same exon in different tissues
exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Remove exons that don't have a downstream intron (last exons)
exons.tab=exons.tab[exons.tab$DownIntronStart>0,]

#In the last GTF there was a single type 0 exon of length 1bp exactly
exons.tab=exons.tab[exons.tab$End>exons.tab$Start,]

#Write the downstream intron coordinates of all the exons into a bed file
write.table(cbind(
  
  paste0('chr',exons.tab$Chr),
  exons.tab$DownIntronStart,
  exons.tab$DownIntronEnd,
  paste0('exon',1:nrow(exons.tab),'_type_',exons.tab$Type),
  exons.tab$End-exons.tab$Start+1,
  unlist(lapply(exons.tab$Name,function(x)strsplit(x,split='\\|')[[1]][4])),
  exons.tab$GeneId
)

,'exons_coord.bed',sep='\t',row.names = F,col.names = F,quote = F)

#Intersect the bed files that was created for each of the two biological samples with the bed file containing the downsteam intron coordinates 
for (idx in 1:2)
  system(paste0('bedtools intersect -wa -wb -a exons_coord.bed -b proseq_',idx,'.bed > exons_proseq_',idx,'_intersection.bed'))





