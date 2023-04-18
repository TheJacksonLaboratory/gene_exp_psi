.. _snakemake:

=========
Snakemake
=========

For convenience, we provide a snakemake script that runs all of the analysis reported in the main manuscript.
See `Sustainable data analysis with Snakemake <https://f1000research.com/articles/10-33>`_ for an introduction 
and the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_ for details.


Some steps of the analysis are probably best performed in an HPC/cluster environment, but we have been able to 
perform the entire analysis on a single desktop with 64Gb RAM and an 8-core Xeon CPU with an Ubuntu linux operating system.
There are several ways of setting up Snakemake. See the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ for 
more information. We present here instructions for setting up Snakemake with miniconda on a linux system for convenience.


.. code-block:: bash
    :caption: download

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh


Follow instructions from the prompt to run the installer to initialize Miniconda3. After this, close the terminal window and open a new shell or run ``source ~/.bashrc``.
```


.. code-block:: bash
    :caption: Initialize conda for bash and update conda if necessary

    conda init bash
    conda update conda




mamba
^^^^^

`mamba <https://github.com/mamba-org/mamba>`_ is a replacement for the `conda` executable that is used by default by Snakemake. 

.. code-block:: bash
    :caption: Install mamba

    conda install -n base -c conda-forge mamba


Install Snakemake
^^^^^^^^^^^^^^^^^

It is sensible to create a conda environment specifically for snakemake:

.. code-block:: bash
    :caption: create conda/mamba virtual environment and install snakemake

    conda create --name snakemake
    conda activate snakemake
    mamba install -c conda-forge -c bioconda snakemake



.. code-block:: bash
    :caption: setup virtual environment

    python3 -m venv gepenv
    source gepenv/bin/activate
    pip install --upgrade pip






As a quick start, the following command will perform a ``dry run`` and show the reasons for performing each step.

.. code-block:: bash
    :caption: snakemake dry run

    snakemake -n -r


