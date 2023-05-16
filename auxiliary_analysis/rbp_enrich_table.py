import functools
from collections import defaultdict
from csv import DictReader
import re
import statistics
from my_latex_table import MyLongTable

class RbpMatrix:
    def __init__(self, matrix_id, rbp_symbol, motif) -> None:
        self._matrix_id = matrix_id.strip()
        self._rbp_symbol = rbp_symbol.strip()
        self._motif = motif.strip()

    @property
    def matrix_id(self):
        return self._matrix_id
    
    @property
    def rbp_symbol(self):
        return self._rbp_symbol
    
    @property
    def motif(self):
        return self._motif
    
    


def get_rbp_matrix_dictionary():
    rbp_d = defaultdict(RbpMatrix)
    with open("input/rbp_matrix_list.txt") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) == 3:
                rbpm = RbpMatrix(matrix_id=fields[0], rbp_symbol=fields[1], motif=fields[2])
                rbp_d[rbpm.matrix_id] = rbpm
            elif len(fields) == 4:
                # some of the lines have an extra tab, but we can correct as follows.
                rbpm = RbpMatrix(matrix_id=fields[0], rbp_symbol=fields[2], motif=fields[3])
                rbp_d[rbpm.matrix_id] = rbpm
            else:
                print(f"[ERROR] Malformed line: {line} with {len(fields)} fields")
                print(fields)
                continue
            rbpm = RbpMatrix(matrix_id=fields[0], rbp_symbol=fields[1], motif=fields[2])
            rbp_d[rbpm.matrix_id] = rbpm
    return rbp_d


@functools.total_ordering
class RbpEnrichement:
    def __init__(self, row, rbp_dict) -> None:
        self._name = row['name']
        self._motif = row['motif'].replace(".pfm", "") # e.g., '915_11565747.pfm',
        if self._motif in rbp_dict:
            rmatrix = rbp_dict.get(self._motif)
            self._symbol = rmatrix.rbp_symbol
        else:
            raise ValueError(f"Could not find symbol for {self._motif}")
        self._t0percentage = row['t0.percentage']
        self._t1percentage = row['t1.percentage']
        self._t2percentage = row['t2.percentage']
        self._t1_vs_t0 = row['t1.vs.t0']
        self._t2_vs_t0 = row['t2.vs.t0']
        self._t2_vs_t1 = row['t2.vs.t1']
        self._t1_vs_t0_p = self.get_p(row['t1.vs.t0.p'])
        self._t2_vs_t0_p = self.get_p(row['t2.vs.t0.p'])
        self._t2_vs_t1_p = self.get_p(row['t2.vs.t1.p'])

    def get_p(self, p):
        p = p.replace("p=", "")
        if p == "0.0000*":
            return 0.0
        regex = '^[0-9.]+$'
        if re.search(regex, p):
            return float(p)
        else:
            return 1.0
        
    def get_min_p(self):
        pvals = [self._t1_vs_t0_p, self._t2_vs_t0_p, self._t2_vs_t1_p]
        pvals = sorted(pvals)
        return pvals[0]
    
    def is_significant(self):
        return 0.05 >= self.get_min_p()


    def __lt__(self, other):
        return self.get_min_p() < other.get_min_p()

    def __eq__(self, other):
        return self._motif == other._motif
    
    @property
    def symbol(self):
        return self._symbol
    
    @property
    def name(self):
        return self._name.replace("_", "\\_")
    
    @property
    def t0perc(self):
        return f"{self._t0percentage}\\%"
    
    @property
    def t1perc(self):
        return f"{self._t1percentage}\\%"
    
    @property
    def t2perc(self):
        return f"{self._t2percentage}\\%"
    
    @property
    def t1_vs_t0_p(self):
        if self._t1_vs_t0_p > 0.05:
            return "n.s."
        return str(self._t1_vs_t0_p)
    
    @property
    def t2_vs_t0_p(self):
        if self._t2_vs_t0_p > 0.05:
            return "n.s."
        return str(self._t2_vs_t0_p)
    
    @property
    def t2_vs_t1_p(self):
        if self._t2_vs_t1_p > 0.05:
            return "n.s."
        return str(self._t2_vs_t1_p)
    
    @property
    def t1_t0_perc_diff(self):
        return abs(float(self._t0percentage) - float(self._t1percentage))
    
    @property
    def t2_t0_perc_diff(self):
        return abs(float(self._t0percentage) - float(self._t2percentage))
    
    @property
    def t1_t1_perc_diff(self):
        return abs(float(self._t2percentage) - float(self._t1percentage))
        
        
        
      



def get_rbp_enrichements(rbp_d):
    rbp_enrichment_list = []
    with open("rbp_enrichments.tsv") as f:
        reader = DictReader(f, delimiter="\t")
        for row in reader:
            rbpE = RbpEnrichement(row=row, rbp_dict=rbp_d)
            rbp_enrichment_list.append(rbpE)
    return rbp_enrichment_list


rbp_d = get_rbp_matrix_dictionary()
print(f"Extracted total of {len(rbp_d)} RBP matrices")
rbp_enrichment_list = get_rbp_enrichements(rbp_d=rbp_d)
rbp_enrichment_list = sorted(rbp_enrichment_list)


header_fields = ["motif", "id", "t0", "t1", "t2", "t1 vs t0", "t2 vs t0", "t2 vs t1"]
caption = ""
table = MyLongTable(header=header_fields)
for item in rbp_enrichment_list:
    row = []
    row.append(item.symbol)
    row.append(item.name)
    row.append(item.t0perc)
    row.append(item.t1perc)
    row.append(item.t2perc)
    row.append(item.t1_vs_t0_p)
    row.append(item.t2_vs_t0_p)
    row.append(item.t2_vs_t1_p)
    table.add_line(row=row)

# print(table.get_table())


rbp_set = set()
matrix_set = set()
rbp_sig_set = set()
t1_t0_diff = []
t1_t0_sig_set = set()
t2_t0_sig_set = set()
t2_t0_diff = []
t2_t1_diff = []
for item in rbp_enrichment_list:
    rbp_set.add(item.symbol)
    matrix_set.add(item.name)
    if item._t1_vs_t0_p <= 0.05:
        t1_t0_diff.append(item.t1_t0_perc_diff)
        rbp_sig_set.add(item.symbol)
        t1_t0_sig_set.add(item.symbol)
    if item._t2_vs_t0_p <= 0.05:
        t2_t0_diff.append(item.t2_t0_perc_diff)
        rbp_sig_set.add(item.symbol)
        t2_t0_sig_set.add(item.symbol)
    if item._t2_vs_t1_p <= 0.05:
        t2_t1_diff.append(item.t2_t1_perc_diff)
        rbp_sig_set.add(item.symbol)



print(f"Total of {len(rbp_set)} RBPs of which {len(rbp_sig_set)} showed a significant difference")
print(f"Total of {len(matrix_set)} matrices")
mt1t0 = statistics.mean(t1_t0_diff)

print(f"T1/T0 n = {len(t1_t0_diff)}, mean {mt1t0} with total {len(t1_t0_sig_set)} significant RBPs")
mt2t0 = statistics.mean(t2_t0_diff)
print(f"T2/T0 n = {len(t1_t0_diff)}, mean {mt2t0} with total {len(t2_t0_sig_set)} significant RBPs")
if len(t2_t1_diff) == 0:
    print("No significant differences between T2 and T1")
else:
    mt2t1 = statistics.mean(t2_t1_diff)
    print(f"T2/T1 n = {len(t2_t1_diff)}, mean {mt2t1}")





