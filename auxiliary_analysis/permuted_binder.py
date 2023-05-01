import os
from collections import defaultdict


class PermutedBinder:
    """
    Class to store the results of the permuted motifs. FABIAN creates lines like this:
      AR@jaspar2022DetailedTFFMs@TFFM0001.1@TFFM0001.1         : 20.1    >I p=0.0000*     II p=     - 
      ARNT2@jaspar2022DetailedTFFMs@TFFM0853.1@TFFM0853.1      :  9.0    <I p=0.0070     <II p=0.0000*
    0.........0.........0.........0.........0.........0.........0.........0.........0.........0.........0
             10        20        30........40........50........60........70........80........90.......100
    The lines are fixed width.
    Additionally, "BASE' indicates the base of the comparison, e.g., BASE:0 means that we compare group I
    and group II against group 0; BASE:1 means we compare group 2 against group 1. The frequency of a group
    is shown is shown according to BASE (e.g., here the frequency of AR in group 0 is 20.1).
    We need to aggregate the results for all three BASE sections to get all of the results for any motif.
    """
    def __init__(self,  base0, base1, base2) -> None:
        names = set() # use a set to make sure all three names are the same
        name, perc, interp1, pval1, interp2,  pval2 = self._get_components(base0)
        names.add(name)
        self._t0_perc = perc
        self._t1_vs_t0_interp = interp1
        self._t1_vs_t0_pval = pval1
        self._t2_vs_t0_interp = interp2
        self._t2_vs_t0_pval = pval2
        name, perc, interp1, pval1, interp2,  pval2 = self._get_components(base1)
        names.add(name)
        self._t1_perc = perc
        self._t2_vs_t1_interp = interp2
        self._t2_vs_t1_pval = pval2
        name, perc, interp1, pval1, interp2,  pval2 = self._get_components(base2)
        names.add(name)
        self._t2_perc = perc
        if len(names) != 1:
            raise ValueError(f"Multiple names for {name}")
        self._fullname = name.strip()
        fields = name.split("@")
        self._motif_name = fields[0].strip()



    def _get_components(self, line): 
        lfields = line.split(":")
        name = lfields[0].strip() # get everything up to the colon and then remove whitespace
        line2 = lfields[1]
        perc = line2[0:8].strip() # percentage present
        interp1 = line2[8:11].strip()
        pval1 = line2[12:25].strip()
        interp2 = line2[25:28].strip()
        pval2 = line2[28:40].strip()
        return name, perc, interp1, pval1, interp2,  pval2
    

    @staticmethod
    def get_header():
        hfields = ['name', 'motif', 't0.percentage', 't1.percentage', 't2.percentage', 't1.vs.t0', 't2.vs.t0', 't2.vs.t1', 't1.vs.t0.p', 't2.vs.t0.p', 't2.vs.t1.p']
        return "\t".join(hfields)

    def get_row(self):
        fields = [self._fullname, self._motif_name, self._t0_perc, self._t1_perc, self._t2_perc, self._t1_vs_t0_interp, self._t2_vs_t0_interp, self._t2_vs_t1_interp,
                   self._t1_vs_t0_pval, self._t2_vs_t0_pval, self._t2_vs_t1_pval ]
        return "\t".join(fields)

     









def process_file(fname, outname):
    fh = open(outname, 'wt')
    fh.write(PermutedBinder.get_header() + "\n")
    base_d = defaultdict(list)
    inBASE = False
    # Note: files have ordering of BASE0, BASE1, BASE2 - we only check that we get three lines for each motif
    # Some transcription factors have multiple motifs
    with open(fname) as f:
        for line in f:
            if 'BASE:' in line:
                inBASE = True
                continue
            if inBASE and ":" in line:
                fields = line.split(":")
                motif = fields[0].strip()
                base_d[motif].append(line)
    for k, v in base_d.items():
        if len(v) != 3:
            print(f"Error: We got {len(v)} lines for {k}")
            exit(1) # should never happen unless there is a serious error that needs to be fixed.
        binder = PermutedBinder(base0=v[0], base1=v[1], base2=v[2])
        fh.write(binder.get_row() + "\n")
    fh.close()
        



dir_path = os.path.dirname(os.path.realpath(__file__))

outname = "cpe_enrichments.tsv"
process_file( os.path.join(dir_path, "permutation_results/splicing_log_0_head1000000.txt"), outname)

outname = "tf_enrichments.tsv"
process_file(os.path.join(dir_path,"permutation_results/splicing_log_1_head1000000.txt"), outname)

outname = "rbp_enrichments.tsv"
process_file(os.path.join(dir_path,"permutation_results/splicing_log_2_head1000000.txt"), outname)
