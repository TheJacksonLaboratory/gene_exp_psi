




class MyLongTable:
    def __init__(self, header, headerformat=None, caption="to do", centered=True):
        self._caption = caption
        self._centered = centered
        if not isinstance(header, list):
            raise ValueError(f"header argument must be a list but was {type(header)}: {header}")
        self._header = header
        # default cell format is l(eft)
        if headerformat is None:
            self._headerformat = "l"*len(header)
        self._hlines = []
        self._lines = []
        self._flines = []
        self._label = "label"


    def _headerlines(self):
        if self._centered:
            self._hlines.append("\\begin{center}\n")
        self._hlines.append("\\begin{longtable}{" + self._headerformat + "}\n")
        self._hlines.append("\\caption{" + self._caption +"}\n")
        self._hlines.append("\\label{tab:" + self._label + "}\n")
        self._hlines.append("\\hline\n")
        hfields = []
        for h in self._header:
            hfields.append("\\multicolumn{1}{c}{\\textbf{" + h + "}}" )
        self._hlines.append(" & ".join(hfields) + "\\\\")
        self._hlines.append("\\hline\n")
        self._hlines.append("\\endfirsthead\n")
        self._hlines.append("\\multicolumn{" + str(len(self._header)) + "}{c}%\n")
        self._hlines.append("{{\\bfseries \\tablename\ \\thetable{} -- continued from previous page}} \\\\ \n")
        self._hlines.append("\\hline\n")
        self._hlines.append(" & ".join(hfields) + "\\\\")
        self._hlines.append("\\hline\n")
        self._hlines.append("\\endhead \n")
        self._hlines.append("\\hline\n")
        self._hlines.append("\\multicolumn{" + str(len(self._header)) + "}{r}{{Continued on next page}} \\\\ \n")
        self._hlines.append("\\hline\n") 
        self._hlines.append("\\endfoot \n")
        self._hlines.append("\\hline\n") 
        self._hlines.append("\\hline\n")
        self._hlines.append("\\endfoot \n")

    def _footerlines(self):
        self._flines.append("\\end{longtable}\n")
        if self._centered:
            self._flines.append("\\end{center}\n")

    def add_line(self, row):
        if len(row) != len(self._header):
            raise ValueError(f"Expected a line with {len(self._header)} fields but got {len(row)}: {row}")
        self._lines.append(" & ".join(row) + "\\\\ \n")

    def get_table(self):
        rows = []
        self._headerlines()
        self._footerlines()
        rows.extend(self._hlines)
        rows.extend(self._lines)
        rows.extend(self._flines)
        return "".join(rows)
    
    def write_table(self, fh):
        lines = self.get_table()
        for line in lines:
            fh.write(line)




class MyLatexTable:
    def __init__(self, header, headerformat=None, caption="to do", centered=True):
        self._caption = caption
        self._centered = centered
        if not isinstance(header, list):
            raise ValueError(f"header argument must be a list but was {type(header)}: {header}")
        self._header = header
        # default cell format is l(eft)
        if headerformat is None:
            self._headerformat = "l"*len(header)
        self._hlines = []
        self._lines = []
        self._flines = []
        self._label = "label"

    def _headerlines(self):
        self._hlines.append("\\begin{table}\n")
        if self._centered:
            self._hlines.append("\\centering\n")
        self._hlines.append("\\begin{tabular}{" + self._headerformat + "}\n")
        self._hlines.append("\\hline\n")
        hfields = []
        for h in self._header:
            hfields.append("\\textbf{" + h + "}" )
        self._hlines.append(" & ".join(hfields) + "\\\\")
        self._hlines.append("\\hline\n")
       
    def add_line(self, row):
        if len(row) != len(self._header):
            raise ValueError(f"Expected a line with {len(self._header)} fields but got {len(row)}: {row}")
        self._lines.append(" & ".join(row) + "\\\\ \n")

    
    def _footerlines(self):
        self._flines.append("\\hline \n")
        self._flines.append("\\end{tabular}\n")
        self._hlines.append("\\caption{" + self._caption +"}\n")
        self._hlines.append("\\label{tab:" + self._label + "}\n")
        self._flines.append("\\end{table}\n")

    def get_table(self):
        rows = []
        self._headerlines()
        self._footerlines()
        rows.extend(self._hlines)
        rows.extend(self._lines)
        rows.extend(self._flines)
        return "".join(rows)
    
    def write_table(self, fh):
        lines = self.get_table()
        for line in lines:
            fh.write(line)
        



