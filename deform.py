#!/usr/bin/env python3
# Running deformation calc.
#
# Version 1.2.0
# 2019.01.15
#
# Author: Liu Jia
#

import sys
import os
from slib.calc import Calc
from slib.material import Material
from slib.output import Output
from slib.solve_S import solve_S_matrix
from slib.solve_trS import solve_trS_matrix
from slib.solve_ABD import solve_ABD_matrix


class Deformation(object):
    def __init__(self, calc_file, out_file="output.dat"):
        if not os.path.isfile(calc_file):
            raise ValueError("calculation file doesn't exist: {}.".format(calc_file))
        
        self._calc = Calc(calc_file)
        self._material = Material(self._calc.material)
        self._op = Output(out_file)

        self._op.title("system")
        if self._calc.system == "b":
            self._op.text("bulk.")
        elif self._calc.system == "l":
            self._op.text("laminate.")
        elif self._calc.system == "s":
            self._op.text("stacking laminates.")
        else:
            raise ValueError("Unexpected system: {}.".format(self._calc.system))
        self._op.text("")
        
        self._op.text("length    = {} m;".format(self._calc.length))
        self._op.text("width     = {} m;".format(self._calc.width))
        if self._calc.system == "b":
            self._op.text("thickness = {} m.".format(self._calc.thickness))
        if self._calc.system == "l":
            self._op.text("theta = {} deg.".format(self._calc.rotate_deg))
        if self._calc.system == "s":
            self._op.text("thickness = {} m;".format(self._calc.thickness))
            self._op.text("theta = {} deg.".format(self._calc.rotate_deg))
    
        self._op.text("")

        if self._calc.system == "b":
            solve_S_matrix(self._calc, self._material, self._op)
        elif self._calc.system == "l":
            solve_trS_matrix(self._calc, self._material, self._op)
        else:
            solve_ABD_matrix(self._calc, self._material, self._op)
        
        self._op.text("")
        self._op.title("end")
        self._op.text("Have a nice day.")
        

if __name__ == "__main__":
    Df = Deformation("calc.dat")
