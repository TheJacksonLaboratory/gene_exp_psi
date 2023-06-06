import functools
import re

@functools.total_ordering
class EnrichedItem:
    def __init__(self, row) -> None:
        self._name = row['name']
        self._motif = row['motif'].replace(".pfm", "") # e.g., '915_11565747.pfm',
        self._symbol = self._name
        self._t0percentage = row['t0.percentage']
        self._t1percentage = row['t1.percentage']
        self._t2percentage = row['t2.percentage']
        # calculate max difference
        d1 = abs(float(self._t1percentage) - float(self._t0percentage))
        d2 = abs(float(self._t2percentage) - float(self._t0percentage))
        d3 = abs(float(self._t2percentage) - float(self._t1percentage))
        self._max_diff = max(d1, max(d2,d3))
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
    
    @property
    def max_diff(self):
        return self._max_diff
        