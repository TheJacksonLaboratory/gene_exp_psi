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


# Note the base path for all GTEX files is https://storage.googleapis.com/gtex_analysis_v8/
GTEX_FILES = ["GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz",
              "GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz"]

resources_dir='resources'

rule all:
  input:  
    "resources/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz", 
    "resources/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
    "resources/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    "resources/after_exon_sig_next.RData",  
    expand("types_{tissue_idx}.RData", tissue_idx = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])


rule downloadReference:
  output:
    gtf="resources/Homo_sapiens.GRCh38.109.gtf"
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

rule download_gtex:
  output:
      "resources/gtex/{gtex_file}"
  params:
      cluster_opts='--mem=12G -t 0:20'
  shell:
      """
      mkdir -p resources/gtex
      cd resources/gtex
      wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/{wildcards.gtex_file} -O {wildcards.gtex_file}
      """ 
 

    
rule download_gtex_annot:
  output:
      "resources/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  shell:
      """
      mkdir -p resources/
      cd resources/
      wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --output-document GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt || error_exit "gtf annotation download failed"
      """


rule create_signatures:
  input:
    gtf="resources/Homo_sapiens.GRCh38.109.gtf",
    gtex="resources/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz"
  output:
    outfile="resources/after_exon_sig_next.RData"
  params:
    cluster_opts='--mem=64G -t 24:00'
  conda:
    "env/signature_env.yaml"
  message:
    "====create signatures===="
  shell:
    "Rscript scripts/create_signatures.R {input.gtf} {input.gtex} {output.outfile}"


rule ge_frac_cor:
  input:
    gtf="resources/Homo_sapiens.GRCh38.109.gtf",
    gtex_transcript_tpm="resources/gtex/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz",
    gtex_sample_attr="resources/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
    after_exon="resources/after_exon_sig_next.RData",
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
