# Running the snakemake pipeline in a SLURM environment.

There are many ways to do this, but we recommend the following steps. 
We assume that users will run the script on a linux cluster using SLURM.


## Set up Anaconda and install snakemake


```
srun --pty -t 2:00:00 --mem 16g -q batch bash
wget https://repo.anaconda.com/archive/Anaconda3-2023.03-Linux-x86_64.sh
bash Anaconda3-2023.03-Linux-x86_64.sh
```

Place anaconda in the PATH and create a virtual environment.

```
export PATH=/some/path/anaconda3/bin:$PATH
conda create --name esvenv
conda activate esvenv 
```

Install mamba and snakemake. 
See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
for more information.

```
conda install -n base -c conda-forge mamba
```

Now install snakemake.

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```


srun --pty -t 2:00:00 -q batch bash
wget https://repo.anaconda.com/archive/Anaconda3-2023.03-Linux-x86_64.sh
bash Anaconda3-2023.03-Linux-x86_64.sh
export PATH=/some/path/anaconda3/bin:$PATH
conda create –name esvenv
conda activate esvenv 
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake