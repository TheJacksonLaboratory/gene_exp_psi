# gene_exp_psi
Gene expression and percent spliced in


## snakemake

For convenience, we provide a snakemake script that recreates the major results of the analysis reported in the main manuscript.
See `Sustainable data analysis with Snakemake <https://f1000research.com/articles/10-33>`_ for an introduction 
and the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_ for details.

## SLURM

we provide a script to submit it to a SLURM system, but it is suitable for any high performance computing that provides a high memory partition.

In order to run it using SLURM, place all the files in the same directory, and submit the file run_types_pipeline.sh using the high_mem partition. 

A singularity container is available on Google drive at this link: 

https://drive.google.com/file/d/1tjEopsK8g6JChiBbI51DXgIKTcgdD0NI/view?usp=sharing

The following definitions file was used for creating the container:
```
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
```


Any environment that supports the same software would run the pipeline. The output is created in the same directory. Images are outputted as PDF files.


# permute

permute is a simple C++ application that performs the permutation analysis described in "Alternative splicing is coupled to gene expression in a subset of variably expressed genes" (See section "Enriched motif testing" in the Methods).

The source code and some additional files can be found in the subdirectory "permutations". To run the code, first download the input file ``motif_locations.zip`` from https://zenodo.org/record/8030743 in the same directory as permute.cc. The motif locations directory contains
the locations and sequences of exons and promoters found to be UHP or DHP (see main manuscript and methods).


To build the program, enter the following command


```
g++ -pthread -o permute permute.cc
```

Alternatively, use the Makefile.


## Running permute

```
./permute <repeat_id={0,1,2}> <num_threads> 1000000
```

Note that 0 stands for core promoter elements (CPE), 1 for transcription factor flexible models (TFFM), and 2 for RNA binding protein (RBP) motifs.
