# README

This analysis assesses an association between type 0, type 1, and type 2 exons and the genes they are contained in with RNA binding proteins.


## Data files
All files were downloaded on 25 April 2023 from STRING (https://string-db.org/).

9606.protein.aliases.v11.5.txt.gz		
9606.protein.enrichment.terms.v11.5.txt.gz	
9606.protein.info.v11.5.txt			
9606.protein.physical.links.v11.5.txt

## Analysis files

1) Get proteins annotated to RNA binding (note that the file included inferred annotations).

```
$ zgrep 'GO:0003723' 9606.protein.enrichment.terms.v11.5.txt.gz > RNAbindingAnnotations.tab
```

The resulting file RNAbindingAnnotations.tab is included in this repository.


2) Get significant motifs (for now just group 1) -- run get_sig_motifs.py
this leads to output file group1.txt

3) get list of proteins that interact with an RBP and that do not