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
resources_dir='resources'

rule all:
  input:
    "after_exon_sig_next.RData"


rule downloadReference:
  output:
    gtf = resources_dir+'/' + HTTP_FILES[0]
  params:
    resources_dir = resources_dir
  shell:
    """
    mkdir -p {params.resources_dir}
    cd {params.resources_dir}
    wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz || error_exit "gtf download failed"
    gunzip Homo_sapiens.GRCh38.109.gtf.gz || error_exit "GTF gunzip failed"
    """
        
rule downloadGtex:
  output:
      gtex_rsem = resources_dir + '/' +  HTTP_FILES[1]
  params:
      resources_dir = resources_dir
  shell:
      """
      mkdir -p {params.resources_dir}
      cd {params.resources_dir}
      wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz  || error_exit "GTEX/RSEM download failed"
      """
    
rule create_signatures:
  input:
    gtf=resources_dir + '/' + HTTP_FILES[0],
    gtex=resources_dir + '/' + HTTP_FILES[1]
  output:
    "after_exon_sig_next.RData"
  conda:
    "env/signature_env.yaml"
  message:
    "====create signatures===="
  shell:
    "Rscript scripts/create_signatures.R {input.gtf} {input.gtex} "