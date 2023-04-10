# gene_exp_psi
Gene expression and percent spliced in


## snakemake

For convenience, we provide a snakemake script that runs all of the analysis reported in the main manuscript.
See `Sustainable data analysis with Snakemake <https://f1000research.com/articles/10-33>`_ for an introduction 
and the `Snakemake documentation <https://snakemake.readthedocs.io/en/stable/>`_ for details.

As a quick start, the following command will perform a ``dry run`` and show the reasons for performing each step.

.. code-block:: bash
    :caption: snakemake dry run

    snakemake -n -r


## Components of the analysis

The snakemake script performs the following steps.

### Downloads

The script downloads the files *GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct.gz* and *Homo_sapiens.GRCh38.109.gtf*, which are
required by analysis scripts. It is possible to download the files manually - in this case, place them in the resources folder or add a softlink to the
location on your system.


### Create signatures

This R scripts creates a data structure to store information about exons and the transcripts they belong to, and stores the information in an R
data file called *after_exon_sig_next.RData*.
