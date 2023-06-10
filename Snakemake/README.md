# README

This directory contain the files for the Snakemake pipeline. 

we provide a script to submit it to a SLURM system, but it is suitable for any high performance computing that provides a high memory partition.

In order to run it using SLURM, place all the files in the same directory, and submit the file run_types_pipeline.sh using the high_mem partition. 

A singularity container is available on Google drive at this link: 

https://drive.google.com/drive/folders/1ql3OoMwh50wvh3BH7vvCiU8-MTOkvmga?usp=share_link

The following definitions file was used for creating the container:

Bootstrap: docker
From: continuumio/miniconda3

%post

/opt/conda/bin/conda update -n base conda
/opt/conda/bin/conda create --name types-pipeline -c conda-forge -c bioconda snakemake-minimal python=3.11 r-base bedtools anaconda::cmake gcc

%runscript
R --slave -e 'install.packages("BiocManager")'
R --slave -e 'BiocManager::install("BRGenomics")'
R --slave -e 'install.packages("statmod","~/R/x86_64-conda_cos6-linux-gnu-library/4.0")'
R --slave -e 'install.packages("rstatix")'
R --slave -e 'install.packages("data.table")'
R --slave -e 'install.packages("intervals")'
R --slave -e 'BiocManager::install("seqinr")'
R --slave -e 'BiocManager::install("ggpubr")'
R --slave -e 'install.packages("R.utils")'



Any environment that supports the same software would run the pipeline. The output is created in the same directory. Images are outputted as PDF files.