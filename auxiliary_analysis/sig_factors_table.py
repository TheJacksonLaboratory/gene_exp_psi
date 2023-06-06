from csv import DictReader

class Entry:
    def __init__(self, row):
        self.name = row['name']
        self.motif = row['motif']
        self.t0percentage = float(row['t0.percentage'])
        self.t1percentage = float(row['t1.percentage'])
        self.t2percentage = float(row['t2.percentage'])
        self.t1vst0 = row['t1.vs.t0']
        self.t2vst0 = row['t2.vs.t0']
        self.t2vst1 = row['t2.vs.t1']
        self.t1vst0p = Entry.get_p(row['t1.vs.t0.p'])
        self.t2vst0p = Entry.get_p(row['t2.vs.t0.p'])
        self.t2vst1p = Entry.get_p(row['t2.vs.t1.p'])

    @staticmethod
    def get_p(p_str):
        p_str = p_str.replace("p=", "")
        p_str = p_str.replace("*", "")
        try:
            return float(p_str)
        except:
            return 1.0

    def min_p(self):
        numbers = [self.t1vst0p, self.t2vst0p, self.t2vst1p]
        numbers.sort()
        return numbers[0]

    def max_difference(self):
        d1 = abs(self.t0percentage - self.t1percentage)
        d2 = abs(self.t0percentage - self.t2percentage)
        return max(d1, d2)


    def __lt__(self, other):
        if self.min_p() == other.min_p():
            return self.max_difference() > other.max_difference()
        return self.min_p() < other.min_p()

    def __str__(self):
        return f"{self.motif}: T0: {self.t0percentage}; T1: {self.t1percentage} T2: {self.t2percentage}, maxdiff: {self.max_difference()}, minp {self.min_p()}"

with open("tf_enrichments.tsv") as f:
    reader = DictReader(f, delimiter='\t')
    entries = []
    for row in reader:
        entries.append(Entry(row=row))


entries_sorted = sorted(entries)

for e in entries_sorted:
    print(e)