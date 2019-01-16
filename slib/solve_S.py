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
from sympy import solve, Matrix
import math
import numpy as np


def solve_S_matrix(c, m, op):
    op.title("solving S matrix")

    op.section("knowns and unknowns")
    op.part("normal strains")
    op.text("epsilon_x = {};".format(c.x_normal_strain), end=" ")
    op.text("epsilon_y = {};".format(c.y_normal_strain), end=" ")
    op.text("epsilon_z = {};".format(c.z_normal_strain))
    op.text("")
    op.part("normal loads")
    op.text("F_x = {} N;".format(c.x_normal_load), end=" ")
    op.text("F_y = {} N;".format(c.y_normal_load), end=" ")
    op.text("F_z = {} N;".format(c.z_normal_load))
    op.text("")
    op.part("normal stresses")
    op.text("sigma_x = {} Pa;".format(c.x_normal_stress), end=" ")
    op.text("sigma_y = {} Pa;".format(c.y_normal_stress), end=" ")
    op.text("sigma_z = {} Pa;".format(c.z_normal_stress))
    op.text("")
    op.part("shear strains")
    op.text("gamma_yz = {};".format(c.yz_shear_load), end=" ")
    op.text("gamma_xz = {};".format(c.xz_shear_load), end=" ")
    op.text("gamma_xy = {};".format(c.xy_shear_load))
    op.text("")
    op.part("shear loads")
    op.text("F_yz = {} N;".format(c.yz_shear_strain), end=" ")
    op.text("F_xz = {} N;".format(c.xz_shear_strain), end=" ")
    op.text("F_xy = {} N;".format(c.xy_shear_strain))
    op.text("")
    op.part("shear stresses")
    op.text("tau_yz = {} Pa;".format(c.yz_shear_stress), end=" ")
    op.text("tau_xz = {} Pa;".format(c.xz_shear_stress), end=" ")
    op.text("tau_xy = {} Pa;".format(c.xy_shear_stress))
    op.text("")

    if c.delta_thm:
        op.part("thermal conditions")
        op.text("alpha_1 = {} 1/C;".format(m.alpha_1), end=" ")
        op.text("alpha_2 = {} 1/C;".format(m.alpha_2), end=" ")
        op.text("alpha_3 = {} 1/C;".format(m.alpha_3), end=" ")
        op.text("delta_T = {} C;".format(c.delta_thm))
        op.text("")
    if c.delta_moi:
        op.part("moisture conditions")
        op.text("beta_1 = {} 1/C;".format(m.beta_1), end=" ")
        op.text("beta_2 = {} 1/C;".format(m.beta_2), end=" ")
        op.text("beta_3 = {} 1/C;".format(m.beta_3), end=" ")
        op.text("delta_M = {} C;".format(c.delta_moi))
        op.text("")

    op.part("S_matrix (Pa)")
    S_matrix = m.get_S_matrix()
    op.matrix(S_matrix)
    op.text("")
    op.part("C_matrix (1/Pa)")
    C_matrix = m.get_C_matrix()
    op.matrix(C_matrix)
    op.text("")

    op.section("solution")
    # solution 1.
    op.part("solve stress vs load")
    stress_1 = Matrix([[c.x_normal_load / (c.width     * c.thickness)],
                       [c.y_normal_load / (c.thickness * c.length   )],
                       [c.z_normal_load / (c.length    * c.width    )],
                       [c.yz_shear_load / (c.width     * c.thickness)],
                       [c.xz_shear_load / (c.thickness * c.length   )],
                       [c.xy_shear_load / (c.length    * c.width    )],])
    
    stress_2 = Matrix([[c.x_normal_stress],
                       [c.y_normal_stress],
                       [c.z_normal_stress],
                       [c.yz_shear_stress],
                       [c.xz_shear_stress],
                       [c.xy_shear_stress],])
    
    solution = solve(stress_1 - stress_2)
    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    # x_nr = c.x_normal_stress if isinstance(c.x_normal_stress, (int, float)) else solution[c.x_normal_stress]
    # y_nr = c.y_normal_stress if isinstance(c.y_normal_stress, (int, float)) else solution[c.y_normal_stress]
    # z_nr = c.z_normal_stress if isinstance(c.z_normal_stress, (int, float)) else solution[c.z_normal_stress]

    x_nr = solution[c.x_normal_stress] if c.x_normal_stress in solution else c.x_normal_stress
    y_nr = solution[c.y_normal_stress] if c.y_normal_stress in solution else c.y_normal_stress
    z_nr = solution[c.z_normal_stress] if c.z_normal_stress in solution else c.z_normal_stress
    
    op.part("normal stresses")
    op.text("sigma_x = {} Pa;".format(x_nr), end=" ")
    op.text("sigma_y = {} Pa;".format(y_nr), end=" ")
    op.text("sigma_z = {} Pa;".format(z_nr))
    op.text("")

    # yz_sr = c.yz_shear_stress if isinstance(c.yz_shear_stress, (int, float)) else solution[c.yz_shear_stress]
    # xz_sr = c.xz_shear_stress if isinstance(c.xz_shear_stress, (int, float)) else solution[c.xz_shear_stress]
    # xy_sr = c.xy_shear_stress if isinstance(c.xy_shear_stress, (int, float)) else solution[c.xy_shear_stress]

    yz_sr = solution[c.yz_shear_stress] if c.yz_shear_stress in solution else c.yz_shear_stress
    xz_sr = solution[c.xz_shear_stress] if c.xz_shear_stress in solution else c.xz_shear_stress
    xy_sr = solution[c.xy_shear_stress] if c.xy_shear_stress in solution else c.xy_shear_stress

    op.part("shear stresses")
    op.text("tau_yz = {} Pa;".format(yz_sr), end=" ")
    op.text("tau_xz = {} Pa;".format(xz_sr), end=" ")
    op.text("tau_xy = {} Pa;".format(xy_sr))
    op.text("")

    # solution 2.
    op.part("solved strain vs stress")
    strain = Matrix([[c.x_normal_strain - m.alpha_1 * c.delta_thm - m.beta_1 * c.delta_moi],
                     [c.y_normal_strain - m.alpha_2 * c.delta_thm - m.beta_2 * c.delta_moi],
                     [c.z_normal_strain - m.alpha_3 * c.delta_thm - m.beta_3 * c.delta_moi],
                     [c.yz_shear_strain                                                   ],
                     [c.xz_shear_strain                                                   ],
                     [c.xy_shear_strain                                                   ],])
    
    stress = Matrix([[x_nr],
                     [y_nr],
                     [z_nr],
                     [yz_sr],
                     [xz_sr],
                     [xy_sr],])
    
    solution = solve(Matrix(S_matrix) * stress - strain)
    
    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    x_ns = c.x_normal_strain if isinstance(c.x_normal_strain, (int, float)) else solution[c.x_normal_strain]
    y_ns = c.y_normal_strain if isinstance(c.y_normal_strain, (int, float)) else solution[c.y_normal_strain]
    z_ns = c.z_normal_strain if isinstance(c.z_normal_strain, (int, float)) else solution[c.z_normal_strain]

    op.part("normal strains")
    op.text("epsilon_x = {};".format(x_ns), end=" ")
    op.text("epsilon_y = {};".format(y_ns), end=" ")
    op.text("epsilon_z = {};".format(z_ns))
    op.text("")

    yz_ss = c.yz_shear_strain if isinstance(c.yz_shear_strain, (int, float)) else solution[c.yz_shear_strain]
    xz_ss = c.xz_shear_strain if isinstance(c.xz_shear_strain, (int, float)) else solution[c.xz_shear_strain]
    xy_ss = c.xy_shear_strain if isinstance(c.xy_shear_strain, (int, float)) else solution[c.xy_shear_strain]

    op.part("shear strains")
    op.text("gamma_yz = {};".format(yz_ss), end=" ")
    op.text("gamma_xz = {};".format(xz_ss), end=" ")
    op.text("gamma_xy = {};".format(xy_ss))
    op.text("")
    
    # solution 3.
    f_l = c.length * (1 + x_ns)
    f_w = c.width * (1 + y_ns)
    f_t = c.thickness * (1 + z_ns)
    f_yz = 90 + yz_ss * 180 / math.pi
    f_xz = 90 + xz_ss * 180 / math.pi
    f_xy = 90 + xy_ss * 180 / math.pi
    
    op.title("final size")
    op.section("final length")
    op.text("length    = {} m;".format(f_l))
    op.text("widht     = {} m;".format(f_w))
    op.text("thickness = {} m;".format(f_t))
    op.text("")
    op.section("final angle")
    op.text("yz angle = {} deg;".format(f_yz))
    op.text("xz angle = {} deg;".format(f_xz))
    op.text("xy angle = {} deg;".format(f_xy))
