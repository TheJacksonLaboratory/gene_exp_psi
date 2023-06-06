import functools
from collections import defaultdict
from csv import DictReader
import re
import statistics
from my_latex_table import MyLatexTable, MyLongTable
from enriched_item import EnrichedItem






def get_enrichements(fname):
    cpe_enrichment_list = []
    with open(fname) as f:
        reader = DictReader(f, delimiter="\t")
        for row in reader:
            item = EnrichedItem(row=row)
            cpe_enrichment_list.append(item)
    return cpe_enrichment_list



cpe_enrichment_list = get_enrichements("cpe_enrichments.tsv")
cpe_enrichment_list = sorted(cpe_enrichment_list)


header_fields = ["motif", "id", "t0", "t1", "t2", "t1 vs t0", "t2 vs t0", "t2 vs t1"]
caption = ""
table = MyLatexTable(header=header_fields)
for item in cpe_enrichment_list:
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
for item in cpe_enrichment_list:
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




print(table.get_table())


print(f"Total of {len(rbp_set)} CPEs of which {len(rbp_sig_set)} showed a significant difference")
print(f"Total of {len(matrix_set)} matrices")
if len(t1_t0_diff) == 0:
    print("No significant differences between T1 and T10")
else:
    mt1t0 = statistics.mean(t1_t0_diff)
    print(f"T1/T0 n = {len(t1_t0_diff)}, mean {mt1t0} with total {len(t1_t0_sig_set)} significant RBPs")
if len(t2_t0_diff) == 0:
    print("No significant differences between T2  and T10")
else:
    mt2t0 = statistics.mean(t2_t0_diff)
    print(f"T2/T0 n = {len(t1_t0_diff)}, mean {mt2t0} with total {len(t2_t0_sig_set)} significant RBPs")
if len(t2_t1_diff) == 0:
    print("No significant differences between T2 and T1")
else:
    mt2t1 = statistics.mean(t2_t1_diff)
    print(f"T2/T1 n = {len(t2_t1_diff)}, mean {mt2t1}")


print("\n\n\n\n\n&&&&&&&&&&&&&&&&&&&&&&&&&\n\n\n\n")

## TFBS

tf_enrichment_list = get_enrichements("tf_enrichments.tsv")
tf_enrichment_list = sorted(tf_enrichment_list)


header_fields = ["motif", "id", "t0", "t1", "t2", "t1 vs t0", "t2 vs t0", "t2 vs t1"]
caption = ""
table = MyLongTable(header=header_fields)
for item in tf_enrichment_list:
    row = []
    # item.symbol is something like this
    # ZNF460@jaspar2022DetailedTFFMs@TFFM0642.1@TFFM0642.1
    fields = item.symbol.split("@")
    symbol = fields[0]
    name = fields[-1]
    row.append(symbol)
    row.append(name)
    row.append(item.t0perc)
    row.append(item.t1perc)
    row.append(item.t2perc)
    row.append(item.t1_vs_t0_p)
    row.append(item.t2_vs_t0_p)
    row.append(item.t2_vs_t1_p)
    table.add_line(row=row)

# print(table.get_table())


tf_set = set()
matrix_set = set()
tf_sig_set = set()
t1_t0_diff = []
t1_t0_sig_set = set()
t2_t0_sig_set = set()
t2_t0_diff = []
t2_t1_diff = []
for item in tf_enrichment_list:
    tf_set.add(item.symbol)
    matrix_set.add(item.name)
    if item._t1_vs_t0_p <= 0.05:
        t1_t0_diff.append(item.t1_t0_perc_diff)
        tf_sig_set.add(item.symbol)
        t1_t0_sig_set.add(item.symbol)
    if item._t2_vs_t0_p <= 0.05:
        t2_t0_diff.append(item.t2_t0_perc_diff)
        tf_sig_set.add(item.symbol)
        t2_t0_sig_set.add(item.symbol)
    if item._t2_vs_t1_p <= 0.05:
        t2_t1_diff.append(item.t2_t1_perc_diff)
        tf_sig_set.add(item.symbol)




print(table.get_table())


print(f"Total of {len(tf_set)} TFBSs of which {len(tf_sig_set)} showed a significant difference")
print(f"Total of {len(matrix_set)} matrices")
if len(t1_t0_diff) == 0:
    print("No significant differences between T1 and T10")
else:
    mt1t0 = statistics.mean(t1_t0_diff)
    print(f"T1/T0 n = {len(t1_t0_diff)}, mean {mt1t0} with total {len(t1_t0_sig_set)} significant TFs")
if len(t2_t0_diff) == 0:
    print("No significant differences between T2  and T10")
else:
    mt2t0 = statistics.mean(t2_t0_diff)
    print(f"T2/T0 n = {len(t1_t0_diff)}, mean {mt2t0} with total {len(t2_t0_sig_set)} significant TFs")
if len(t2_t1_diff) == 0:
    print("No significant differences between T2 and T1")
else:
    mt2t1 = statistics.mean(t2_t1_diff)
    print(f"T2/T1 n = {len(t2_t1_diff)}, mean {mt2t1}")


## Show top ten TFBS
## according to maximum difference


tf_enrichment_list.sort(key=lambda r: r.max_diff, reverse=True)


header_fields = ["motif", "id", "t0", "t1", "t2", "t1 vs t0", "t2 vs t0", "t2 vs t1"]
caption = ""
table = MyLatexTable(header=header_fields)
for item in tf_enrichment_list[:10]:
    row = []
     # item.symbol is something like this
    # ZNF460@jaspar2022DetailedTFFMs@TFFM0642.1@TFFM0642.1
    fields = item.symbol.split("@")
    symbol = fields[0]
    name = fields[-1]
    row.append(symbol)
    row.append(name)
    row.append(item.t0perc)
    row.append(item.t1perc)
    row.append(item.t2perc)
    row.append(item.t1_vs_t0_p)
    row.append(item.t2_vs_t0_p)
    row.append(item.t2_vs_t1_p)
    table.add_line(row=row)

print(table.get_table())