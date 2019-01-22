#!/usr/bin/env python3
# Output results.
#
# Version 1.0.0
# 2019.01.13
#
# Author: Liu Jia
#


import os
import numpy as np


class Output(object):
    def __init__(self, output_file):
        self._file = output_file
        with open(self._file, "w") as f:
            f.write("Deform info.\n\n")
            f.write("Version 1.5.0\n")
            f.write("Calculation start.\n\n")
    
    def title(self, info):
        with open(self._file, "a") as f:
            f.write("===== ===== ===== ===== {} ===== ===== ===== =====\n".format(info[0].upper() + info[1:]))
    
    def section(self, info):
        with open(self._file, "a") as f:
            f.write("=== === {} === ===\n".format(info[0].upper() + info[1:]))
    
    def part(self, info):
        with open(self._file, "a") as f:
            f.write("Pr) {}\n".format(info[0].upper() + info[1:]))

    def text(self, info, end="\n"):
        with open(self._file, "a") as f:
            f.write("{}{}".format(info, end))

    def matrix(self, info):
        with open(self._file, "a") as f:
            m = np.array(info)
            s = m.shape
            for i in range(s[0]):
                f.write("[")
                for j in range(s[1]):
                    fmt = "{:.3E}".format(m[i][j])
                    if fmt[0] == "-":
                        f.write("{}   ".format(fmt))
                    else:
                        f.write(" {}   ".format(fmt))
                f.write("]\n")

    def draft(self, info):
        with open(self._file, "a") as f:
            f.write("--- -- - {} - -- ---\n".format(info[0].upper() + info[1:]))