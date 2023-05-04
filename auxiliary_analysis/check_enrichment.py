from collections import defaultdict
from argparse import ArgumentParser
from csv import DictReader


def get_gene_symbols(fname):
    """
    Parse the 9606.protein.info.v11.5.txt file, which has the following structure
    #string_protein_id      preferred_name  protein_size    annotation
    9606.ENSP00000000233    ARF5    180     ADP-ribosylation factor 5; GTP-binding protein that functions as (...)
    9606.ENSP00000000412    M6PR    277     Cation-dependent mannose-6-phosphate receptor; Transport of  (...)
    The goal is the get a mapping from the ENSP accession numbers to the gene symbols
    """
    name_d = defaultdict()
    with open(fname) as f:
        for line in f:
            fields = line.split("\t")
            ensp = fields[0]
            symbol = fields[1]
            name_d[ensp] = symbol
    return name_d


def get_rnabp():
    """
    Parse the RNAbindingAnnotations.tab that is located in the same folder and was prepared as described in the
    README.md file. The file contains a list of proteins annotation to "RNA Binding" in Gene Ontology.
    """
    fname = "RNAbindingAnnotations.tab"
    rbp = set()
    with open(fname) as f:
        for line in f:
            fields = line.split("\t")
            ensp = fields[0]
            rbp.add(ensp)
    return rbp


def get_ppis(fname):
    """
    Parse the STRING 9606.protein.physical.links.v11.5.txt file to obtain high quality (score at least 700)
    protein-protein interactions.
    fname -- path to "9606.protein.physical.links.v11.5.txt"
    """
    interactions = []
    with open(fname) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip().split(" ")
            if len(fields) != 3:
                print("Bad line: " + line)
                continue
            score = int(fields[2])
            if score >= 700:
                i = [fields[0], fields[1]]
                interactions.append(i)
    return interactions

def get_significant_and_nonsig_tfbs():
    """
    This script assumes that the command permuted_binder.py has been executed first and the file tf_enrichemnt.tab
    has been created.
    """
    sig_motifs = set()
    non_sig_motifs = set()
    with open("tf_enrichments.tsv") as f:
        reader = DictReader(f, delimiter='\t')
        for row in reader:
            motif = row['motif']
            t1_vs_t0 = row['t1.vs.t0']
            t2_vs_t0 = row['t2.vs.t0']
            t2_vs_t1 = row['t2.vs.t1']
            if ">" in t2_vs_t0 or "<" in t2_vs_t0 or ">" in t1_vs_t0 or "<" in t1_vs_t0:
                sig_motifs.add(motif)
            else:
                non_sig_motifs.add(motif)
    return sig_motifs, non_sig_motifs

def get_rbp_interactions_and_list_of_proteins(interaction_list, ensp2sym_d, rbp_set):
    """
    Go through the interaction list from STRING. Record a set of
    all proteins that take part in interactions, and also record the
    proteins that interact with an RNA-binding protein
    """
    all_proteins = set()
    rbp_binders = set()

    for i in interaction_list:
        p1 = i[0]
        p2 = i[1]
        if p1 in ensp2sym_d and p2 in ensp2sym_d:
            s1 = ensp2sym_d.get(p1)
            s2 = ensp2sym_d.get(p2)
            if p1 in rbp_set:
                rbp_binders.add(s2)
            if p2 in rbp_set:
                rbp_binders.add(s1)
            all_proteins.add(s1)
            all_proteins.add(s2)
    return rbp_binders, all_proteins


parser = ArgumentParser()
parser.add_argument("--protein_info", required=True, help="path to 9606.protein.info.v11.5.txt ")
parser.add_argument("--physical", required=True, help="path to 9606.protein.physical.links.v11.5.txt")
args = parser.parse_args()
print(args.protein_info)


ensp2sym = get_gene_symbols(fname=args.protein_info)
rbp_set = get_rnabp()
interactions = get_ppis(fname=args.physical)

rbp_binders, all_proteins = get_rbp_interactions_and_list_of_proteins(interaction_list=interactions,
                                                                      ensp2sym_d=ensp2sym, rbp_set=rbp_set)



n_prot = len(all_proteins)
n_rbp_binders = len(rbp_binders)

sig_motifs, non_sig_motifs = get_significant_and_nonsig_tfbs()



print(f"Got {len(interactions)} PPIs")
print(f"Got {len(ensp2sym)} ENSEMBL to SYm maps")
print(f"Got {len(rbp_set)} RBPs")
print(f"Got {len(non_sig_motifs)} TFs of which {len(sig_motifs)} bind RBPs")
global_rate = 100 * len(rbp_set) / len(all_proteins)
print(f"Global RBP binding rate {global_rate:.1f}%")
tf_rate = 100 * len(sig_motifs) / len(non_sig_motifs)
print(f"Group 1/2 TF RBP binding rate {tf_rate:.1f}%")
