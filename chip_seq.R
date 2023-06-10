#read the exons tab

exons.tab=read.table('exons_tab.txt',sep='\t',header=T)

#Remove rows that correspond to tissues in which type I/II exons were classified as type 0

exons.tab=exons.tab[(exons.tab$Type!=0) | (!exons.tab$Name %in% exons.tab$Name[exons.tab$Type!=0]),]

#Keep only genes that have at least one type I or type II exon

exons.tab=exons.tab[exons.tab$GeneId %in% exons.tab$GeneId[exons.tab$Type!=0],]

#Remove rows that correspond to the same exon in different tissues

exons.tab=exons.tab[!duplicated(exons.tab$Name),]

#Write the exon coordinates of all the exons into a bed file

write.table(cbind(
  
  paste0('chr',exons.tab$Chr),
  exons.tab$Start,
  exons.tab$End,
  paste0('exon',1:nrow(exons.tab),'_type_',exons.tab$Type),
  exons.tab$End-exons.tab$Start+1,
  unlist(lapply(exons.tab$Name,function(x)strsplit(x,split='\\|')[[1]][4])),
  exons.tab$GeneId
)

,'exons_coord_chip.bed',sep='\t',row.names = F,col.names = F,quote = F)

#Intersect the bed file of each sample with the bed file containing the exon coordinates 

for (bedfile in c('ENCFF836KEF.bed','ENCFF616GPO.bed','ENCFF750UAQ.bed','ENCFF061DKX.bed','ENCFF502RVO.bed','ENCFF886PSD.bed','ENCFF890XQW.bed','ENCFF680DFX.bed','ENCFF538YRS.bed','ENCFF843PJA.bed','ENCFF681FKP.bed','ENCFF159PYD.bed','ENCFF216PEY.bed','ENCFF473YGW.bed','ENCFF848IHI.bed','ENCFF039TJY.bed','ENCFF834UYS.bed','ENCFF885ZWB.bed','ENCFF211PWC.bed','ENCFF634JRD.bed','ENCFF730QNU.bed','ENCFF258NNN.bed','ENCFF235UTX.bed','ENCFF519DDH.bed','ENCFF745RUH.bed','ENCFF870KJE.bed','ENCFF218GFN.bed','ENCFF456FTV.bed','ENCFF322DAE.bed','ENCFF394NVG.bed','ENCFF043FHJ.bed','ENCFF535TAL.bed','ENCFF341YAZ.bed','ENCFF410CPY.bed','ENCFF658XFZ.bed','ENCFF292PBX.bed','ENCFF129APS.bed','ENCFF195PBD.bed','ENCFF909CTD.bed','ENCFF166EUJ.bed','ENCFF245CAL.bed','ENCFF918UHB.bed','ENCFF946CSJ.bed','ENCFF107SJD.bed','ENCFF607OBG.bed','ENCFF712JXT.bed','ENCFF921FKB.bed','ENCFF898ZLY.bed','ENCFF448ZOJ.bed','ENCFF880CLF.bed','ENCFF618BDP.bed','ENCFF451XZC.bed','ENCFF952WTH.bed','ENCFF396ZQL.bed','ENCFF672HZG.bed','ENCFF698HVQ.bed','ENCFF064ASM.bed','ENCFF226BZE.bed','ENCFF493EPF.bed','ENCFF105LDY.bed','ENCFF082DKP.bed','ENCFF683YQZ.bed','ENCFF779EDF.bed','ENCFF869HZM.bed','ENCFF567ELA.bed','ENCFF132MQB.bed','ENCFF838LKZ.bed','ENCFF536KQX.bed','ENCFF389ZBU.bed','ENCFF151LHR.bed','ENCFF923UMY.bed','ENCFF471ZTS.bed','ENCFF175PDD.bed','ENCFF673HLW.bed','ENCFF558UJR.bed','ENCFF615EAT.bed','ENCFF160KVW.bed','ENCFF227ELA.bed','ENCFF371CYY.bed','ENCFF574UKD.bed','ENCFF569MJS.bed','ENCFF376CIT.bed','ENCFF673VSN.bed','ENCFF657NRR.bed','ENCFF653SLQ.bed','ENCFF229KAY.bed','ENCFF065MQK.bed','ENCFF608MSL.bed','ENCFF540RGI.bed','ENCFF411PMA.bed','ENCFF770ENO.bed','ENCFF872EBX.bed','ENCFF254MJA.bed','ENCFF993GPP.bed','ENCFF627GBZ.bed','ENCFF482BDW.bed','ENCFF012VQA.bed','ENCFF355MNE.bed','ENCFF842JME.bed','ENCFF681MRC.bed','ENCFF555UBC.bed','ENCFF469KDX.bed','ENCFF900JDD.bed','ENCFF967TFR.bed','ENCFF354VWZ.bed'))
  system(paste0('bedtools intersect -wa -wb -a exons_coord_chip.bed -b ',bedfile,' > exons_chip_intersection_',bedfile))

