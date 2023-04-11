shell.prefix("""
# http://linuxcommand.org/wss0150.php
PROGNAME=$(basename $0)

function error_exit
{{
#   ----------------------------------------------------------------
#   Function for exit due to fatal program error
#       Accepts 1 argument:
#           string containing descriptive error message
#   ----------------------------------------------------------------
    echo "${{PROGNAME}}: ${{1:-"Unknown Error"}}" 1>&2
    exit 1
}}
""")

HTTP_FILES = ["Homo_sapiens.GRCh38.109.gtf",
              "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz"]


GTEX_FILES = ["GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz",
              "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
              "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"]

resources_dir='resources'

rule all:
  input:
    "after_exon_sig_next.RData",
    expand("resources/{g}", g=GTEX_FILES),
    expand("types_{tissue_idx}.RData", tissue_idx = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])


rule downloadReference:
  output:
    gtf = resources_dir+'/' + HTTP_FILES[0]
  params:
    cluster_opts='--mem=36G -t 0:20',
    resources_dir = resources_dir
  shell:
    """
    mkdir -p {params.resources_dir}
    cd {params.resources_dir}
    wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz || error_exit "gtf download failed"
    gunzip Homo_sapiens.GRCh38.109.gtf.gz || error_exit "GTF gunzip failed"
    """
        
rule downloadGtexExpected:
  output:
    "resources/" + GTEX_FILES[0]
  params:
      cluster_opts='--mem=12G -t 0:20',
      resources_dir = resources_dir
  shell:
      """
      mkdir -p {resources_dir}
      cd {resources_dir}
      wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz  || error_exit "GTEX/RSEM download failed!"
      """

rule downloadGtexTpm:
  output:
    "resources/" + GTEX_FILES[1]
  params:
      cluster_opts='--mem=12G -t 0:20',
      resources_dir = resources_dir
  shell:
      """
      mkdir -p {resources_dir}
      cd {resources_dir}
      wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz  || error_exit "GTEX/RSEM TPM download failed!"
      """



rule downloadGtexSampleAttributes:
  output:
    "resources/" + GTEX_FILES[2]
  params:
      cluster_opts='--mem=12G -t 0:20',
      resources_dir = resources_dir
  shell:
      """
      mkdir -p {resources_dir}
      cd {resources_dir}
      wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt  || error_exit "GTEX sample attribs download failed!"
      """

    

    
rule create_signatures:
  input:
    gtf=resources_dir + '/' + HTTP_FILES[0],
    gtex=resources_dir + '/' + GTEX_FILES[0]
  output:
    "after_exon_sig_next.RData"
  params:
    cluster_opts='--mem=64G -t 24:00'
  conda:
    "env/signature_env.yaml"
  message:
    "====create signatures===="
  shell:
    "Rscript scripts/create_signatures.R {input.gtf} {input.gtex} "


rule ge_frac_cor:
  input:
    gtf=resources_dir + '/' + HTTP_FILES[0],
    gtex_transcript_tpm=resources_dir + '/' + GTEX_FILES[1],
    gtex_sample_attr=resources_dir + '/' + GTEX_FILES[2],
    after_exon="after_exon_sig_next.RData",
   # idx="{tissue_idx}"
  output:
    "types_{tissue_idx}.RData"
  params:
    cluster_opts='--mem=64G -t 4:00'
  conda:
    "env/gefrac_env.yaml"
  message:
    "====ge_frac_cor===="
  shell:
    "Rscript scripts/ge_frac_cor.R {input.gtf} {input.gtex_transcript_tpm} {input.gtex_sample_attr} {input.after_exon} {wildcards.tissue_idx} 2"
