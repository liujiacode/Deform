#!/usr/bin/env python3
# Loading and paraphrasing calculation inputs.
#
# Version 1.4.0
# 2019.01.13
#
# Author: Liu Jia
#


import sys
import os
import re
from slib.pstrlst import del_comment, del_stri_ht_space
from sympy import Symbol
import numpy as np


class Calc(object):
    def __init__(self, input_file):
        if not os.path.isfile(input_file):
            raise ValueError("Calculation input file doesn't exist: {}.".format(input_file))
        
        self.material = "material.dat"
        self.system = "b"
        self.lay_up_sequence = None
        
        self.x_normal_load = Symbol("x_nl")
        self.y_normal_load = Symbol("y_nl")
        self.z_normal_load = Symbol("z_nl")

        self.x_normal_stress = Symbol("x_nr")
        self.y_normal_stress = Symbol("y_nr")
        self.z_normal_stress = Symbol("z_nr")
        
        self.yz_shear_stress = Symbol("yz_sr")
        self.xz_shear_stress = Symbol("xz_sr")
        self.xy_shear_stress = Symbol("xy_sr")
        
        self.x_normal_strain = Symbol("x_ns")
        self.y_normal_strain = Symbol("y_ns")
        self.z_normal_strain = Symbol("z_ns")

        self.yz_shear_load = Symbol("yz_sl")
        self.xz_shear_load = Symbol("xz_sl")
        self.xy_shear_load = Symbol("xy_sl")
        
        self.yz_shear_strain = Symbol("yz_ss")
        self.xz_shear_strain = Symbol("xz_ss")
        self.xy_shear_strain = Symbol("xy_ss")

        self.x_bending_moment = Symbol("x_bm")
        self.y_bending_moment = Symbol("y_bm")
        self.xy_twisting_moment = Symbol("xy_tm")
        
        self.x_surface_curvature = Symbol("x_sc")
        self.y_surface_curvature = Symbol("y_sc")
        self.xy_twisting_curvature = Symbol("xy_tc")
        

        self.length = None
        self.width = None
        self.thickness = None

        self.constrain_x = False
        self.constrain_y = False
        self.constrain_z = False

        self.rotate_deg = 0

        self.delta_thm = 0
        self.delta_moi = 0
        

        # Load property file.
        props = open(input_file)
        for line in props:
            pline = del_comment(line)
            
            # Skip comment line.
            if not pline:
                continue
            
            # Paraphras properties.
            words = pline.split("=")
            if len(words) != 2:
                raise ValueError("Expect one equation in each line: {}.".format(pline))
            variable = del_stri_ht_space(words[0])
            value = del_stri_ht_space(words[1])

            # Loading.
            if variable == "material":
                if os.path.isfile(value):
                    self.material = value
                else:
                    raise ValueError("Material doesn't exist: {}.".format(value))

            elif variable == "system":
                if value[0].lower() in ("b", "l", "s"):
                    self.system = value[0].lower()
                else:
                    raise ValueError("Unknown system: {}.".format(value))
            
            elif variable == "x_normal_load":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.x_normal_load = self._N_or_n(value)
            elif variable == "y_normal_load":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.y_normal_load = self._N_or_n(value)
            elif variable == "z_normal_load":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.z_normal_load = self._N_or_n(value)
            
            elif variable == "x_normal_strain":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.x_normal_strain = self._N_or_n(value)
            elif variable == "y_normal_strain":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.y_normal_strain = self._N_or_n(value)
            elif variable == "z_normal_strain":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.z_normal_strain = self._N_or_n(value)
            
            elif variable == "x_normal_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.x_normal_stress = self._N_or_n(value)
            elif variable == "y_normal_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.y_normal_stress = self._N_or_n(value)
            elif variable == "z_normal_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.z_normal_stress = self._N_or_n(value)
            
            elif variable == "yz_shear_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.yz_shear_stress = self._N_or_n(value)
            elif variable == "xz_shear_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xz_shear_stress = self._N_or_n(value)
            elif variable == "xy_shear_stress":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xy_shear_stress = self._N_or_n(value)
            
            elif variable in ("yz_shear_load", "zy_shear_load"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.yz_shear_load = self._N_or_n(value)
            elif variable in ("xz_shear_load", "zx_shear_load"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xz_shear_load = self._N_or_n(value)
            elif variable in ("xy_shear_load", "yx_shear_load"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xy_shear_load = self._N_or_n(value)
            
            elif variable in ("yz_shear_strain", "zy_shear_strain"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.yz_shear_strain = self._N_or_n(value)
            elif variable in ("xz_shear_strain", "zx_shear_strain"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xz_shear_strain = self._N_or_n(value)
            elif variable in ("xy_shear_strain", "yx_shear_strain"):
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xy_shear_strain = self._N_or_n(value)
            
            elif variable == "x_bending_moment":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.x_bending_moment = self._N_or_n(value)
            elif variable == "y_bending_moment":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.y_bending_moment = self._N_or_n(value)
            elif variable == "xy_twisting_moment":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xy_twisting_moment = self._N_or_n(value)
            
            elif variable == "x_surface_curvature":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.x_surface_curvature = self._N_or_n(value)
            elif variable == "y_surface_curvature":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.y_surface_curvature = self._N_or_n(value)
            elif variable == "xy_twisting_curvature":
                if isinstance(self._N_or_n(value), (int, float)):
                    self.xy_twisting_curvature = self._N_or_n(value)
            
            elif variable == "length":
                self.length = float(value)
            elif variable == "width":
                self.width = float(value)
            elif variable == "thickness":
                if "/" in value:
                    self.thickness = [float(i) for i in value.split("/")]
                else:
                    self.thickness = float(value)
            
            elif variable == "constrain_x":
                self.constrain_x = self._T_or_F(value)
            elif variable == "constrain_y":
                self.constrain_y = self._T_or_F(value)
            elif variable == "constrain_z":
                self.constrain_z = self._T_or_F(value)

            elif variable == "rotate_deg":
                if "/" in value:
                    self.rotate_deg = [float(i) for i in value.split("/")]
                else:
                    self.rotate_deg = float(value)
            
            elif variable == "delta_thm":
                self.delta_thm = float(value)
            elif variable == "delta_moi":
                self.delta_moi = float(value)

            else:
                print("Warnning: unknown variable: {}.".format(variable))
                continue
            
        props.close()

        # Confirm variables.
        if not self.length or self.length <= 0:
            raise ValueError("Unexpect length: {}.".format(self.length))
        if not self.width or self.width <= 0:
            raise ValueError("Unexpect width: {}.".format(self.width))
        
        # Constrain in directions.
        if self.constrain_x:
            print("Note: x direction is constrained.")
            print("x strain is zero.")
            self.x_normal_strain = 0
            self.x_normal_load = Symbol("x_nl")
            self.x_normal_stress = Symbol("x_nr")
        if self.constrain_y:
            print("Note: y direction is constrained.")
            print("y strain is zero.")
            self.y_normal_strain = 0
            self.y_normal_load = Symbol("y_nl")
            self.y_normal_stress = Symbol("y_nr")
        if self.constrain_z:
            print("Note: z direction is constrained.")
            print("z strain is zero.")
            self.z_normal_strain = 0
            self.z_normal_load = Symbol("z_nl")
            self.z_normal_stress = Symbol("z_nr")

        # if self.constrain_x or self.constrain_y or self.constrain_z:
        #     self.yz_shear_load = self.xz_shear_load = self.xy_shear_load = 0
        #     self.yz_shear_stress = self.xz_shear_stress = self.xy_shear_stress = 0
        #     self.yz_shear_strain = self.xz_shear_strain = self.xy_shear_strain = 0
        
        if self.system == "l":
            print("Note: plane stress and strain conditions are used.")
            print("z stresses are zeros (not considered).")
            self.z_normal_strain = Symbol("z_ns")
            self.yz_shear_strain = 0
            self.xz_shear_strain = 0
            self.z_normal_load = 0
            self.yz_shear_load = 0
            self.xz_shear_load = 0
            self.z_normal_stress = 0
            self.yz_shear_stress = 0
            self.xz_shear_stress = 0
        
        if self.system == "s":
            print("Note: classical lamination theory and Kickhoff hypothesis are used.")
            print("z stresses and strains are zeros (not considered).")
            self.z_normal_strain = 0
            self.yz_shear_strain = 0
            self.xz_shear_strain = 0
            self.z_normal_load = 0
            self.yz_shear_load = 0
            self.xz_shear_load = 0
            self.z_normal_stress = 0
            self.yz_shear_stress = 0
            self.xz_shear_stress = 0
            print("delta thm and moi are zeros (not considered).")
            self.delta_moi = 0
            self.delta_thm = 0

            if isinstance(self.rotate_deg, (int, float)):
                self.rotate_deg =np.array([self.rotate_deg,])
            else:
                self.rotate_deg = np.array(self.rotate_deg)
            if isinstance(self.thickness, (int, float)):
                scale_num = len(self.rotate_deg)
                scale_len = scale_num * self.thickness
                self.thickness = np.linspace(-1 * scale_len / 2, scale_len / 2, num=scale_num + 1)
            else:
                self.thickness = np.array(self.thickness)
            
            if len(self.thickness) != len(self.rotate_deg) + 1:
                raise ValueError("length of thicknesses and of rotate degrees are not coincident.")
            
    # Return None or number.
    def _lay_up(self, stri):
        lsti = stri.split("/")
        lay = [float(i) for i in lsti]
        return lay
    
    # Return None or number.
    def _N_or_n(self, stri):
        if stri[0].lower() in ("n", "u"):
            return None
        else:
            return float(stri)
    
    def _T_or_F(self, stri):
        if stri[0].lower() == "t":
            return True
        elif stri[0].lower() == "f":
            return False
        else:
            raise ValueError("Unexpect bool variable: {}.".format(stri))
    