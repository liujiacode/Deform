#!/usr/bin/env python3
# Running deformation calc.
#
# Version 1.0.0
# 2019.01.13
#
# Author: Liu Jia
#

import sys
import os
from slib.calc import Calc
from slib.material import Material
from slib.output import Output
from sympy import solve, Matrix
import math
import numpy as np


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
        
        self._op.text("length    = {};".format(self._calc.length))
        self._op.text("width     = {};".format(self._calc.width))
        if self._calc.system == "b":
            self._op.text("thickness = {}.".format(self._calc.thickness))
        if self._calc.system == "l":
            self._op.text("theta = {}.".format(self._calc.rotate_deg))
        if self._calc.system == "s":
            self._op.text("thickness = {}.".format(self._calc.thickness))
            self._op.text("theta = {}.".format(self._calc.rotate_deg))
    
        self._op.text("")

        if self._calc.system == "b":
            self._solve_S_matrix(self._calc, self._material)
        elif self._calc.system == "l":
            self._solve_trS_matrix(self._calc, self._material)
        else:
            self._solve_ABD_matrix(self._calc, self._material)
        
        self._op.text("")
        self._op.title("end")
        self._op.text("Have a nice day.")
        
    def _solve_S_matrix(self, c, m):
        self._op.title("solving S matrix")
        self._op.section("known")
        self._op.part("normal strains")
        self._op.text("epsilon_x = {};".format(c.x_normal_strain), end=" ")
        self._op.text("epsilon_y = {};".format(c.y_normal_strain), end=" ")
        self._op.text("epsilon_z = {};".format(c.z_normal_strain))
        self._op.text("")
        self._op.part("normal loads")
        self._op.text("sigma_x = {} N;".format(c.x_normal_load), end=" ")
        self._op.text("sigma_y = {} N;".format(c.y_normal_load), end=" ")
        self._op.text("sigma_z = {} N;".format(c.z_normal_load))
        self._op.text("")
        self._op.part("shear strains")
        self._op.text("gamma_yz = {};".format(c.yz_shear_load), end=" ")
        self._op.text("gamma_xz = {};".format(c.xz_shear_load), end=" ")
        self._op.text("gamma_xy = {};".format(c.xy_shear_load))
        self._op.text("")
        self._op.part("shear loads")
        self._op.text("tau_yz = {} N;".format(c.yz_shear_strain), end=" ")
        self._op.text("tau_xz = {} N;".format(c.xz_shear_strain), end=" ")
        self._op.text("tau_xy = {} N;".format(c.xy_shear_strain))
        self._op.text("")
        if c.delta_thm:
            self._op.part("thermal conditions")
            self._op.text("alpha_1 = {} 1/C;".format(m.alpha_1), end=" ")
            self._op.text("alpha_2 = {} 1/C;".format(m.alpha_2), end=" ")
            self._op.text("alpha_3 = {} 1/C;".format(m.alpha_3), end=" ")
            self._op.text("delta_T = {} C;".format(c.delta_thm), end=" ")
            self._op.text("")
        if c.delta_moi:
            self._op.part("moisture conditions")
            self._op.text("beta_1 = {} 1/C;".format(m.beta_1), end=" ")
            self._op.text("beta_2 = {} 1/C;".format(m.beta_2), end=" ")
            self._op.text("beta_3 = {} 1/C;".format(m.beta_3), end=" ")
            self._op.text("delta_M = {} C;".format(c.delta_moi), end=" ")
            self._op.text("")
        self._op.part("S_matrix (Pa)")
        S_matrix = m.get_S_matrix()
        self._op.matrix(S_matrix)
        self._op.text("")

        self._op.section("solution")

        stress = Matrix([[c.x_normal_load / (c.width     * c.thickness)],
                         [c.y_normal_load / (c.thickness * c.length   )],
                         [c.z_normal_load / (c.length    * c.width    )],
                         [c.yz_shear_load / (c.width     * c.thickness)],
                         [c.xz_shear_load / (c.thickness * c.length   )],
                         [c.xy_shear_load / (c.length    * c.width    )],])
        
        strain = Matrix([[c.x_normal_strain - m.alpha_1 * c.delta_thm - m.beta_1 * c.delta_moi],
                         [c.y_normal_strain - m.alpha_2 * c.delta_thm - m.beta_2 * c.delta_moi],
                         [c.z_normal_strain - m.alpha_3 * c.delta_thm - m.beta_3 * c.delta_moi],
                         [c.yz_shear_strain                                                   ],
                         [c.xz_shear_strain                                                   ],
                         [c.xy_shear_strain                                                   ],])
        solution = solve(Matrix(S_matrix) * stress - strain)
        
        self._op.part("solved unknowns")
        for i in solution:
            self._op.text("{} = {};".format(i, solution[i]))
        self._op.text("")

        final_c = {}
        x_ns = c.x_normal_strain
        y_ns = c.y_normal_strain
        z_ns = c.z_normal_strain
        
        if not isinstance(x_ns, (int, float)):
            x_ns = solution[c.x_normal_strain]
        if not isinstance(y_ns, (int, float)):
            y_ns = solution[c.y_normal_strain]
        if not isinstance(z_ns, (int, float)):
            z_ns = solution[c.z_normal_strain]
        
        final_c["x_ns"] = x_ns
        final_c["y_ns"] = y_ns 
        final_c["z_ns"] = z_ns

        yz_ss = c.yz_shear_strain
        xz_ss = c.xz_shear_strain
        xy_ss = c.xy_shear_strain

        if not isinstance(yz_ss, (int, float)):
            yz_ss = solution[c.yz_shear_strain]
        if not isinstance(xz_ss, (int, float)):
            xz_ss = solution[c.xz_shear_strain]
        if not isinstance(xy_ss, (int, float)):
            xy_ss = solution[c.xy_shear_strain]
        
        final_c["yz_ss"] = yz_ss
        final_c["xz_ss"] = xz_ss 
        final_c["xy_ss"] = xy_ss

        f_l = c.length * (1 + x_ns)
        f_w = c.width * (1 + y_ns)
        f_t = c.thickness * (1 + z_ns)
        f_yz = 90 + yz_ss * 180 / math.pi
        f_xz = 90 + xz_ss * 180 / math.pi
        f_xy = 90 + xy_ss * 180 / math.pi
        
        self._op.title("final size")
        self._op.section("final length")
        self._op.text("length    = {} m;".format(f_l))
        self._op.text("widht     = {} m;".format(f_w))
        self._op.text("thickness = {} m;".format(f_t))
        self._op.text("")

        self._op.section("final angle")
        self._op.text("yz angle = {} deg;".format(f_yz))
        self._op.text("xz angle = {} deg;".format(f_xz))
        self._op.text("xy angle = {} deg;".format(f_xy))

        return final_c
    
    def _solve_trS_matrix(self, c, m):
        self._op.title("solving transform reduced S matrix")
        self._op.section("known")
        self._op.part("normal strains")
        self._op.text("epsilon_x = {};".format(c.x_normal_strain), end=" ")
        self._op.text("epsilon_y = {};".format(c.y_normal_strain))
        self._op.text("")
        self._op.part("normal loads")
        self._op.text("sigma_x = {} N;".format(c.x_normal_load), end=" ")
        self._op.text("sigma_y = {} N;".format(c.y_normal_load))
        self._op.text("")
        self._op.part("shear strains")
        self._op.text("gamma_xy = {};".format(c.xy_shear_load))
        self._op.text("")
        self._op.part("shear loads")
        self._op.text("tau_xy = {} N;".format(c.xy_shear_strain))
        self._op.text("")
        if c.delta_thm:
            self._op.part("thermal conditions")
            self._op.text("alpha_1 = {} 1/C;".format(m.alpha_1), end=" ")
            self._op.text("alpha_2 = {} 1/C;".format(m.alpha_2), end=" ")
            self._op.text("delta_T = {} C;".format(c.delta_thm))
            self._op.text("")
        if c.delta_moi:
            self._op.part("moisture conditions")
            self._op.text("beta_1 = {} 1/C;".format(m.beta_1), end=" ")
            self._op.text("beta_2 = {} 1/C;".format(m.beta_2), end=" ")
            self._op.text("delta_M = {} C;".format(c.delta_moi))
            self._op.text("")
        
        self._op.part("rotate degree")
        self._op.text("theta = {} deg;".format(c.rotate_deg))
        self._op.text("")
        self._op.part("transformed reduced S_matrix (Pa)")
        trS_matrix = m.get_trans_reduced_S_matrix(rotate_deg = c.rotate_deg)
        self._op.matrix(trS_matrix)
        self._op.text("")

        self._op.section("solution")

        stress = Matrix([[c.x_normal_load / (c.width     * c.thickness)],
                         [c.y_normal_load / (c.thickness * c.length   )],
                         [c.xy_shear_load / (c.length    * c.width    )],])
        
        strain = Matrix([[c.x_normal_strain - m.alpha_1 * c.delta_thm - m.beta_1 * c.delta_moi],
                         [c.y_normal_strain - m.alpha_2 * c.delta_thm - m.beta_2 * c.delta_moi],
                         [c.xy_shear_strain                                                   ],])
        solution = solve(Matrix(trS_matrix) * stress - strain)
        
        self._op.part("solved unknowns")
        for i in solution:
            self._op.text("{} = {};".format(i, solution[i]))
        self._op.text("")

        final_c = {}
        x_ns = c.x_normal_strain
        y_ns = c.y_normal_strain
        
        if not isinstance(x_ns, (int, float)):
            x_ns = solution[c.x_normal_strain]
        if not isinstance(y_ns, (int, float)):
            y_ns = solution[c.y_normal_strain]
        
        final_c["x_ns"] = x_ns
        final_c["y_ns"] = y_ns 

        xy_ss = c.xy_shear_strain

        if not isinstance(xy_ss, (int, float)):
            xy_ss = solution[c.xy_shear_strain]
        
        final_c["xy_ss"] = xy_ss

        f_l = c.length * (1 + x_ns)
        f_w = c.width * (1 + y_ns)
        f_xy = 90 + xy_ss * 180 / math.pi
        
        self._op.title("final size")
        self._op.section("final length")
        self._op.text("length = {} m;".format(f_l))
        self._op.text("widht  = {} m;".format(f_w))
        self._op.text("")
        
        self._op.section("final angle")
        self._op.text("xy angle = {} deg;".format(f_xy))

        return final_c
    
    def _solve_ABD_matrix(self, c, m):
        self._op.title("solving ABD matrix")
        self._op.section("known")
        self._op.part("normal strains")
        self._op.text("epsilon_x = {};".format(c.x_normal_strain), end=" ")
        self._op.text("epsilon_y = {};".format(c.y_normal_strain))
        self._op.text("")
        self._op.part("normal stress")
        self._op.text("N_x = {} N/m2;".format(c.x_normal_load), end=" ")
        self._op.text("N_y = {} N/m2;".format(c.y_normal_load))
        self._op.text("")
        self._op.part("shear strains")
        self._op.text("gamma_xy = {};".format(c.xy_shear_load))
        self._op.text("")
        self._op.part("shear stress")
        self._op.text("N_xy = {} N/m2;".format(c.xy_shear_strain))
        self._op.text("")
        self._op.part("bending moment")
        self._op.text("M_x = {} N/m;".format(c.x_bending_moment), end=" ")
        self._op.text("M_y = {} N/m;".format(c.y_bending_moment))
        self._op.text("")
        self._op.part("surface curvature")
        self._op.text("kappa_x = {} 1/m;".format(c.x_surface_curvature), end=" ")
        self._op.text("kappa_y = {} 1/m;".format(c.y_surface_curvature))
        self._op.text("")
        self._op.part("twisting moment")
        self._op.text("M_xy = {} N/m;".format(c.xy_twisting_moment))
        self._op.text("")
        self._op.part("twisting curvature")
        self._op.text("kappa_xy = {} 1/m;".format(c.xy_twisting_curvature))
        self._op.text("")
        
        self._op.part("rotate degree")
        self._op.text("theta = {} deg;".format(c.rotate_deg))
        self._op.text("")
        self._op.part("ABD_matrix (**)")
        ABD_matrix = m.get_ABD_matrix(c.thickness, c.rotate_deg)
        self._op.matrix(ABD_matrix)
        self._op.text("")

        self._op.section("solution")

        stress = Matrix([[c.x_normal_load],
                         [c.y_normal_load],
                         [c.xy_shear_load],
                         [c.x_bending_moment],
                         [c.y_bending_moment],
                         [c.xy_twisting_moment],])
        
        strain = Matrix([[c.x_normal_strain],
                         [c.y_normal_strain],
                         [c.xy_shear_strain],
                         [c.x_surface_curvature],
                         [c.y_surface_curvature],
                         [c.xy_twisting_curvature],])
        solution = solve(Matrix(ABD_matrix) * strain - stress)

        self._op.part("solved unknowns")
        for i in solution:
            self._op.text("{} = {};".format(i, solution[i]))
        self._op.text("")

        final_c = {}
        x_ns = c.x_normal_strain
        y_ns = c.y_normal_strain
        
        if not isinstance(x_ns, (int, float)):
            x_ns = solution[c.x_normal_strain]
        if not isinstance(y_ns, (int, float)):
            y_ns = solution[c.y_normal_strain]
        
        final_c["x_ns"] = x_ns
        final_c["y_ns"] = y_ns 

        xy_ss = c.xy_shear_strain

        if not isinstance(xy_ss, (int, float)):
            xy_ss = solution[c.xy_shear_strain]
        
        final_c["xy_ss"] = xy_ss

        x_sc = c.x_surface_curvature
        y_sc = c.y_surface_curvature
        
        if not isinstance(x_sc, (int, float)):
            x_sc = solution[c.x_surface_curvature]
        if not isinstance(y_sc, (int, float)):
            y_sc = solution[c.y_surface_curvature]
        
        final_c["x_sc"] = x_sc
        final_c["y_sc"] = y_sc

        xy_tc = c.xy_twisting_curvature

        if not isinstance(xy_tc, (int, float)):
            xy_tc = solution[c.xy_twisting_curvature]
        
        final_c["xy_tc"] = xy_tc

        epsilon_xs = []
        epsilon_ys = []
        gamma_xys = []

        upper_sigma_xs = []
        upper_sigma_ys = []
        upper_tau_xys = []

        lower_sigma_xs = []
        lower_sigma_ys = []
        lower_tau_xys = []

        self._op.title("off axis strains solution")
        for i in range(len(c.thickness)):
            self._op.section("scale {}".format(i + 1))
            
            z_scale = c.thickness[i]
            self._op.text("z scale = {}.".format(z_scale))
            
            epsilon_x = x_ns + z_scale * x_sc
            epsilon_y = y_ns + z_scale * y_sc
            gamma_xy = xy_ss + z_scale * xy_tc

            epsilon_xs.append(epsilon_x)
            epsilon_ys.append(epsilon_y)
            gamma_xys.append(gamma_xy)

            self._op.part("off strain")
            self._op.text("off epsilon_x = {}.".format(epsilon_x))
            self._op.text("off epsilon_y = {}.".format(epsilon_y))
            self._op.text("off gamma_xy = {}.".format(gamma_xy))
            self._op.text("")

            if i > 0:
                theta = c.rotate_deg[i - 1]
                upper_trQ = m.get_trans_reduced_Q_matrix(rotate_deg=theta)
                stress = upper_trQ.dot(np.matrix([[epsilon_x], [epsilon_y], [gamma_xy],]))

                upper_sigma_xs.append(stress[0, 0])
                upper_sigma_ys.append(stress[1, 0])
                upper_tau_xys.append(stress[2, 0])

                self._op.part("off upper layer stress")
                self._op.text("theta = {}.".format(theta))
                self._op.text("off sigma_x = {}.".format(stress[0, 0]))
                self._op.text("off sigma_y = {}.".format(stress[1, 0]))
                self._op.text("off tau_xy = {}.".format(stress[2, 0]))
                self._op.text("")
            if i < len(c.thickness) - 1:
                theta = c.rotate_deg[i]
                upper_trQ = m.get_trans_reduced_Q_matrix(rotate_deg=theta)
                stress = upper_trQ.dot(np.matrix([[epsilon_x], [epsilon_y], [gamma_xy],]))

                lower_sigma_xs.append(stress[0, 0])
                lower_sigma_ys.append(stress[1, 0])
                lower_tau_xys.append(stress[2, 0])

                self._op.part("off lower layer stress")
                self._op.text("theta = {}.".format(theta))
                self._op.text("off sigma_x = {}.".format(stress[0, 0]))
                self._op.text("off sigma_y = {}.".format(stress[1, 0]))
                self._op.text("off tau_xy = {}.".format(stress[2, 0]))
                self._op.text("")
        
        self._op.title("on axis strains solution")
        for i in range(len(c.thickness)):
            self._op.section("scale {}".format(i + 1))
            
            z_scale = c.thickness[i]
            self._op.text("z scale = {}.".format(z_scale))

            if i > 0:
                theta = c.rotate_deg[i - 1]
                m = math.cos(theta * math.pi / 180)
                n = math.sin(theta * math.pi / 180)

                T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                      [ n ** 2,    m ** 2, -2 * m * n      ],
                                      [-1 * m * n, m * n,   m ** 2 - n ** 2],])
                
                self._op.part("on upper layer strain")

                self._op.text("theta = {}.".format(theta))
                strain = T_matrix.dot(np.matrix([[epsilon_xs[i]], [epsilon_ys[i]], [1 / 2 * gamma_xys[i]]]))
                self._op.text("on epsilon_x = {}.".format(strain[0, 0]))
                self._op.text("on epsilon_y = {}.".format(strain[1, 0]))
                self._op.text("on gamma_xy = {}.".format(strain[2, 0] * 2))
                self._op.text("")
            
            if i < len(c.thickness) - 1:
                theta = c.rotate_deg[i]
                m = math.cos(theta * math.pi / 180)
                n = math.sin(theta * math.pi / 180)

                T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                      [ n ** 2,    m ** 2, -2 * m * n      ],
                                      [-1 * m * n, m * n,   m ** 2 - n ** 2],])
                
                self._op.part("on lower layer strain")

                self._op.text("theta = {}.".format(theta))
                strain = T_matrix.dot(np.matrix([[epsilon_xs[i]], [epsilon_ys[i]], [1 / 2 * gamma_xys[i]]]))
                self._op.text("on epsilon_x = {}.".format(strain[0, 0]))
                self._op.text("on epsilon_y = {}.".format(strain[1, 0]))
                self._op.text("on gamma_xy = {}.".format(strain[2, 0] * 2))
                self._op.text("")
            
            if i > 0:
                theta = c.rotate_deg[i - 1]
                m = math.cos(theta * math.pi / 180)
                n = math.sin(theta * math.pi / 180)

                T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                      [ n ** 2,    m ** 2, -2 * m * n      ],
                                      [-1 * m * n, m * n,   m ** 2 - n ** 2],])
                
                self._op.part("on upper layer stress")

                self._op.text("theta = {}.".format(theta))
                stress = T_matrix.dot(np.matrix([[upper_sigma_xs[i - 1]], [upper_sigma_ys[i - 1]], [upper_tau_xys[i - 1]]]))
                self._op.text("on sigma_x = {}.".format(stress[0, 0]))
                self._op.text("on sigma_y = {}.".format(stress[1, 0]))
                self._op.text("on tau_xy = {}.".format(stress[2, 0]))
                self._op.text("")
            
            if i < len(c.thickness) - 1:
                theta = c.rotate_deg[i]
                m = math.cos(theta * math.pi / 180)
                n = math.sin(theta * math.pi / 180)

                T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                      [ n ** 2,    m ** 2, -2 * m * n      ],
                                      [-1 * m * n, m * n,   m ** 2 - n ** 2],])
                
                self._op.part("on lower layer stress")

                self._op.text("theta = {}.".format(theta))
                stress = T_matrix.dot(np.matrix([[upper_sigma_xs[i]], [upper_sigma_ys[i]], [upper_tau_xys[i]]]))
                self._op.text("on sigma_x = {}.".format(stress[0, 0]))
                self._op.text("on sigma_y = {}.".format(stress[1, 0]))
                self._op.text("on tau_xy = {}.".format(stress[2, 0]))
                self._op.text("")

        return final_c

        

if __name__ == "__main__":
    Df = Deformation("calc.dat")
