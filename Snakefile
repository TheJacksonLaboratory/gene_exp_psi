import os
import tempfile
from snakemake.shell import shell

rule all:
    input:
       'WissinkOverlap.pdf',
       'TranscriptsPerGene.pdf',
       'ExonsLengths.pdf',
       'UpstreamIntronLengths.pdf',
       'DownstreamIntronLengths.pdf',
       'TranscriptBiotypes.pdf',
       'TypeIexample.pdf',
       'TypeIIexample.pdf',
       'ChIPSeqPlot.pdf'

rule get_gtex_and_gtf:
    output:
        gtex_counts="GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
        gtex_meta="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        gtf="Homo_sapiens.GRCh38.91.gtf"
    threads: 1
    run:
        shell("wget https://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz")
        shell("gunzip Homo_sapiens.GRCh38.91.gtf.gz")
        shell("wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz")
        shell("wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

rule create_signatures:
      input:
        gtex_counts="GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
        gtex_meta="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        gtf="Homo_sapiens.GRCh38.91.gtf"
      output:
        sigs='after_exon_sig_next.RData'
      threads:
        30
      run:
        shell("Rscript create_signatures.R")
        
rule run_gtex_tissues:
      input:
        sigs='after_exon_sig_next.RData'
      output:
        expand("types_{tissues}.RData",tissues=['Lung','Spleen','Thyroid','Brain - Cortex','Adrenal Gland','Breast - Mammary Tissue','Heart - Left Ventricle',
          'Liver','Pituitary','Pancreas'])
      threads:
        30
      run:
        run_string=''
        for i in [1,2,3,4,5,6,7,8,9,10]:
           run_string=run_string+'Rscript ge_frac_cor.R '+str(i)+' &'
        shell(run_string+'\nwait')
        

rule analyze_exons:
     input:
       expand("types_{tissues}.RData",tissues=['Lung','Spleen','Thyroid','Brain - Cortex','Adrenal Gland','Breast - Mammary Tissue','Heart - Left Ventricle',
          'Liver','Pituitary','Pancreas'])
     output:
       'exons_tab.txt'
     threads:
      1
     run:
       shell("Rscript analyze_exons.R")


rule run_wissink_analysis:
    input:
        i='exons_tab.txt'
    output:
        o='WissinkOverlap.pdf'
    threads:
        1
    run:
        shell("wget https://www.encodeproject.org/files/ENCFF098OEZ/@@download/ENCFF098OEZ.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF672CXH/@@download/ENCFF672CXH.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF251WTK/@@download/ENCFF251WTK.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF857RTX/@@download/ENCFF857RTX.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF944ENF/@@download/ENCFF944ENF.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF480XBL/@@download/ENCFF480XBL.bam")
	shell("wget https://www.encodeproject.org/files/ENCFF710QMM/@@download/ENCFF710QMM.bam")
        shell("wget https://www.encodeproject.org/files/ENCFF933PEI/@@download/ENCFF933PEI.bam")
        shell("Rscript pro_seq.R")
	shell("Rscript summarize_proseq.R") 

   	
rule compare_num_exons:
    input:
       i='exons_tab.txt'
    output:
       o='ExonsPerGene.pdf'
    threads:
       1
    run:
       shell("Rscript number_of_exons_per_gene.R")


rule compare_transcripts_per_gene:
    input:
       i='exons_tab.txt'
    output:
       o='TranscriptsPerGene.pdf'
    threads:
       1
    run:
       shell("Rscript transcripts_per_gene.R")


rule run_exon_lengths:
     input:
       i='exons_tab.txt'
     output:
       o='ExonsLengths.pdf'
     threads:
       1
     run:
       shell("Rscript exon_lengths.R") 


rule run_intron_lengths:
     input:
       i='exons_tab.txt'
     output:
       o1='UpstreamIntronLengths.pdf',
       o2='DownstreamIntronLengths.pdf'
     threads:
       1
     run:
       shell("Rscript intron_lengths.R")


rule transcript_biotypes:
     input:
       i='exons_tab.txt'
     output:
       o='TranscriptBiotypes.pdf'
     threads:
       1
     run:
       shell("Rscript transcript_biotypes.R")

rule plot_examples:
     input:
        gtex_counts="GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
        gtex_meta="GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
        gtf="Homo_sapiens.GRCh38.91.gtf"
     output:
        "TypeIexample.pdf",
        "TypeIIexample.pdf"
     threads:
         1
     run:
         shell("Rscript exp_psi_plots.R")

rule chip_seq:
     input:
       i='exons_tab.txt'
     output:
       o='ChIPSeqPlot.pdf'
     threads:
       1  
     run:
          shell("wget https://www.encodeproject.org/files/ENCFF836KEF/@@download/ENCFF836KEF.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF616GPO/@@download/ENCFF616GPO.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF750UAQ/@@download/ENCFF750UAQ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF061DKX/@@download/ENCFF061DKX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF502RVO/@@download/ENCFF502RVO.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF886PSD/@@download/ENCFF886PSD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF890XQW/@@download/ENCFF890XQW.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF680DFX/@@download/ENCFF680DFX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF538YRS/@@download/ENCFF538YRS.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF843PJA/@@download/ENCFF843PJA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF681FKP/@@download/ENCFF681FKP.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF159PYD/@@download/ENCFF159PYD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF216PEY/@@download/ENCFF216PEY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF473YGW/@@download/ENCFF473YGW.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF848IHI/@@download/ENCFF848IHI.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF039TJY/@@download/ENCFF039TJY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF834UYS/@@download/ENCFF834UYS.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF885ZWB/@@download/ENCFF885ZWB.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF211PWC/@@download/ENCFF211PWC.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF634JRD/@@download/ENCFF634JRD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF730QNU/@@download/ENCFF730QNU.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF258NNN/@@download/ENCFF258NNN.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF235UTX/@@download/ENCFF235UTX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF519DDH/@@download/ENCFF519DDH.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF745RUH/@@download/ENCFF745RUH.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF870KJE/@@download/ENCFF870KJE.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF218GFN/@@download/ENCFF218GFN.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF456FTV/@@download/ENCFF456FTV.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF322DAE/@@download/ENCFF322DAE.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF394NVG/@@download/ENCFF394NVG.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF043FHJ/@@download/ENCFF043FHJ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF535TAL/@@download/ENCFF535TAL.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF341YAZ/@@download/ENCFF341YAZ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF410CPY/@@download/ENCFF410CPY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF658XFZ/@@download/ENCFF658XFZ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF292PBX/@@download/ENCFF292PBX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF129APS/@@download/ENCFF129APS.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF195PBD/@@download/ENCFF195PBD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF909CTD/@@download/ENCFF909CTD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF166EUJ/@@download/ENCFF166EUJ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF245CAL/@@download/ENCFF245CAL.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF918UHB/@@download/ENCFF918UHB.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF946CSJ/@@download/ENCFF946CSJ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF107SJD/@@download/ENCFF107SJD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF607OBG/@@download/ENCFF607OBG.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF712JXT/@@download/ENCFF712JXT.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF921FKB/@@download/ENCFF921FKB.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF898ZLY/@@download/ENCFF898ZLY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF448ZOJ/@@download/ENCFF448ZOJ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF880CLF/@@download/ENCFF880CLF.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF618BDP/@@download/ENCFF618BDP.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF451XZC/@@download/ENCFF451XZC.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF952WTH/@@download/ENCFF952WTH.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF396ZQL/@@download/ENCFF396ZQL.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF672HZG/@@download/ENCFF672HZG.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF698HVQ/@@download/ENCFF698HVQ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF064ASM/@@download/ENCFF064ASM.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF226BZE/@@download/ENCFF226BZE.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF493EPF/@@download/ENCFF493EPF.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF105LDY/@@download/ENCFF105LDY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF082DKP/@@download/ENCFF082DKP.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF683YQZ/@@download/ENCFF683YQZ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF779EDF/@@download/ENCFF779EDF.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF869HZM/@@download/ENCFF869HZM.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF567ELA/@@download/ENCFF567ELA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF132MQB/@@download/ENCFF132MQB.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF838LKZ/@@download/ENCFF838LKZ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF536KQX/@@download/ENCFF536KQX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF389ZBU/@@download/ENCFF389ZBU.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF151LHR/@@download/ENCFF151LHR.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF923UMY/@@download/ENCFF923UMY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF471ZTS/@@download/ENCFF471ZTS.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF175PDD/@@download/ENCFF175PDD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF673HLW/@@download/ENCFF673HLW.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF558UJR/@@download/ENCFF558UJR.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF615EAT/@@download/ENCFF615EAT.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF160KVW/@@download/ENCFF160KVW.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF227ELA/@@download/ENCFF227ELA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF371CYY/@@download/ENCFF371CYY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF574UKD/@@download/ENCFF574UKD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF569MJS/@@download/ENCFF569MJS.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF376CIT/@@download/ENCFF376CIT.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF673VSN/@@download/ENCFF673VSN.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF657NRR/@@download/ENCFF657NRR.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF653SLQ/@@download/ENCFF653SLQ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF229KAY/@@download/ENCFF229KAY.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF065MQK/@@download/ENCFF065MQK.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF608MSL/@@download/ENCFF608MSL.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF540RGI/@@download/ENCFF540RGI.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF411PMA/@@download/ENCFF411PMA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF770ENO/@@download/ENCFF770ENO.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF872EBX/@@download/ENCFF872EBX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF254MJA/@@download/ENCFF254MJA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF993GPP/@@download/ENCFF993GPP.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF627GBZ/@@download/ENCFF627GBZ.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF482BDW/@@download/ENCFF482BDW.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF012VQA/@@download/ENCFF012VQA.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF355MNE/@@download/ENCFF355MNE.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF842JME/@@download/ENCFF842JME.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF681MRC/@@download/ENCFF681MRC.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF555UBC/@@download/ENCFF555UBC.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF469KDX/@@download/ENCFF469KDX.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF900JDD/@@download/ENCFF900JDD.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF967TFR/@@download/ENCFF967TFR.bed.gz")
          shell("wget https://www.encodeproject.org/files/ENCFF354VWZ/@@download/ENCFF354VWZ.bed.gz")
          shell("gunzip ENCFF836KEF.bed.gz")
          shell("gunzip ENCFF616GPO.bed.gz")
          shell("gunzip ENCFF750UAQ.bed.gz")
          shell("gunzip ENCFF061DKX.bed.gz")
          shell("gunzip ENCFF502RVO.bed.gz")
          shell("gunzip ENCFF886PSD.bed.gz")
          shell("gunzip ENCFF890XQW.bed.gz")
          shell("gunzip ENCFF680DFX.bed.gz")
          shell("gunzip ENCFF538YRS.bed.gz")
          shell("gunzip ENCFF843PJA.bed.gz")
          shell("gunzip ENCFF681FKP.bed.gz")
          shell("gunzip ENCFF159PYD.bed.gz")
          shell("gunzip ENCFF216PEY.bed.gz")
          shell("gunzip ENCFF473YGW.bed.gz")
          shell("gunzip ENCFF848IHI.bed.gz")
          shell("gunzip ENCFF039TJY.bed.gz")
          shell("gunzip ENCFF834UYS.bed.gz")
          shell("gunzip ENCFF885ZWB.bed.gz")
          shell("gunzip ENCFF211PWC.bed.gz")
          shell("gunzip ENCFF634JRD.bed.gz")
          shell("gunzip ENCFF730QNU.bed.gz")
          shell("gunzip ENCFF258NNN.bed.gz")
          shell("gunzip ENCFF235UTX.bed.gz")
          shell("gunzip ENCFF519DDH.bed.gz")
          shell("gunzip ENCFF745RUH.bed.gz")
          shell("gunzip ENCFF870KJE.bed.gz")
          shell("gunzip ENCFF218GFN.bed.gz")
          shell("gunzip ENCFF456FTV.bed.gz")
          shell("gunzip ENCFF322DAE.bed.gz")
          shell("gunzip ENCFF394NVG.bed.gz")
          shell("gunzip ENCFF043FHJ.bed.gz")
          shell("gunzip ENCFF535TAL.bed.gz")
          shell("gunzip ENCFF341YAZ.bed.gz")
          shell("gunzip ENCFF410CPY.bed.gz")
          shell("gunzip ENCFF658XFZ.bed.gz")
          shell("gunzip ENCFF292PBX.bed.gz")
          shell("gunzip ENCFF129APS.bed.gz")
          shell("gunzip ENCFF195PBD.bed.gz")
          shell("gunzip ENCFF909CTD.bed.gz")
          shell("gunzip ENCFF166EUJ.bed.gz")
          shell("gunzip ENCFF245CAL.bed.gz")
          shell("gunzip ENCFF918UHB.bed.gz")
          shell("gunzip ENCFF946CSJ.bed.gz")
          shell("gunzip ENCFF107SJD.bed.gz")
          shell("gunzip ENCFF607OBG.bed.gz")
          shell("gunzip ENCFF712JXT.bed.gz")
          shell("gunzip ENCFF921FKB.bed.gz")
          shell("gunzip ENCFF898ZLY.bed.gz")
          shell("gunzip ENCFF448ZOJ.bed.gz")
          shell("gunzip ENCFF880CLF.bed.gz")
          shell("gunzip ENCFF618BDP.bed.gz")
          shell("gunzip ENCFF451XZC.bed.gz")
          shell("gunzip ENCFF952WTH.bed.gz")
          shell("gunzip ENCFF396ZQL.bed.gz")
          shell("gunzip ENCFF672HZG.bed.gz")
          shell("gunzip ENCFF698HVQ.bed.gz")
          shell("gunzip ENCFF064ASM.bed.gz")
          shell("gunzip ENCFF226BZE.bed.gz")
          shell("gunzip ENCFF493EPF.bed.gz")
          shell("gunzip ENCFF105LDY.bed.gz")
          shell("gunzip ENCFF082DKP.bed.gz")
          shell("gunzip ENCFF683YQZ.bed.gz")
          shell("gunzip ENCFF779EDF.bed.gz")
          shell("gunzip ENCFF869HZM.bed.gz")
          shell("gunzip ENCFF567ELA.bed.gz")
          shell("gunzip ENCFF132MQB.bed.gz")
          shell("gunzip ENCFF838LKZ.bed.gz")
          shell("gunzip ENCFF536KQX.bed.gz")
          shell("gunzip ENCFF389ZBU.bed.gz")
          shell("gunzip ENCFF151LHR.bed.gz")
          shell("gunzip ENCFF923UMY.bed.gz")
          shell("gunzip ENCFF471ZTS.bed.gz")
          shell("gunzip ENCFF175PDD.bed.gz")
          shell("gunzip ENCFF673HLW.bed.gz")
          shell("gunzip ENCFF558UJR.bed.gz")
          shell("gunzip ENCFF615EAT.bed.gz")
          shell("gunzip ENCFF160KVW.bed.gz")
          shell("gunzip ENCFF227ELA.bed.gz")
          shell("gunzip ENCFF371CYY.bed.gz")
          shell("gunzip ENCFF574UKD.bed.gz")
          shell("gunzip ENCFF569MJS.bed.gz")
          shell("gunzip ENCFF376CIT.bed.gz")
          shell("gunzip ENCFF673VSN.bed.gz")
          shell("gunzip ENCFF657NRR.bed.gz")
          shell("gunzip ENCFF653SLQ.bed.gz")
          shell("gunzip ENCFF229KAY.bed.gz")
          shell("gunzip ENCFF065MQK.bed.gz")
          shell("gunzip ENCFF608MSL.bed.gz")
          shell("gunzip ENCFF540RGI.bed.gz")
          shell("gunzip ENCFF411PMA.bed.gz")
          shell("gunzip ENCFF770ENO.bed.gz")
          shell("gunzip ENCFF872EBX.bed.gz")
          shell("gunzip ENCFF254MJA.bed.gz")
          shell("gunzip ENCFF993GPP.bed.gz")
          shell("gunzip ENCFF627GBZ.bed.gz")
          shell("gunzip ENCFF482BDW.bed.gz")
          shell("gunzip ENCFF012VQA.bed.gz")
          shell("gunzip ENCFF355MNE.bed.gz")
          shell("gunzip ENCFF842JME.bed.gz")
          shell("gunzip ENCFF681MRC.bed.gz")
          shell("gunzip ENCFF555UBC.bed.gz")
          shell("gunzip ENCFF469KDX.bed.gz")
          shell("gunzip ENCFF900JDD.bed.gz")
          shell("gunzip ENCFF967TFR.bed.gz")
          shell("gunzip ENCFF354VWZ.bed.gz")
          shell("Rscript chip_seq.R")
          shell("Rscript summarize_chipseq.R")
          
