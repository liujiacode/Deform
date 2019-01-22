#!/usr/bin/env python3
# Running deformation calc.
#
# Version 1.4.0
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


def solve_trS_matrix(c, m, op):
    op.title("solving transform reduced S matrix")

    op.section("knowns and unknowns")
    op.part("normal strains")
    op.text("epsilon_x = {};".format(c.x_normal_strain), end=" ")
    op.text("epsilon_y = {};".format(c.y_normal_strain))
    op.text("")
    op.part("normal loads")
    op.text("F_x = {} N;".format(c.x_normal_load), end=" ")
    op.text("F_y = {} N;".format(c.y_normal_load))
    op.text("")
    op.part("normal stresses")
    op.text("sigma_x = {} Pa;".format(c.x_normal_stress), end=" ")
    op.text("sigma_y = {} Pa;".format(c.y_normal_stress))
    op.text("")
    op.part("shear strains")
    op.text("gamma_xy = {};".format(c.xy_shear_strain))
    op.text("")
    op.part("shear loads")
    op.text("F_xy = {} N;".format(c.xy_shear_load))
    op.text("")
    op.part("shear stresses")
    op.text("tau_xy = {} Pa;".format(c.xy_shear_stress))
    op.text("")
    op.part("rotate degree")
    op.text("theta = {} deg;".format(c.rotate_deg))
    op.text("")
    op.part("m = cos theta")
    _cm = math.cos(c.rotate_deg * math.pi / 180)
    op.text("m = {};".format(_cm))
    op.text("")
    op.part("n = sin theta")
    _cn = math.sin(c.rotate_deg * math.pi / 180)
    op.text("n = {};".format(_cn))
    op.text("")

    if c.delta_thm:
        op.part("thermal conditions")
        op.text("alpha_1 = {} 1/C;".format(m.alpha_1))
        op.text("alpha_2 = {} 1/C;".format(m.alpha_2))
        op.text("")
        alpha_x = _cm ** 2 * m.alpha_1 + _cn ** 2 * m.alpha_2
        op.text("alpha_x = m ^ 2 * alpha_1 + n ^ 2 * alpha_2 = {} 1/C;".format(alpha_x))
        alpha_y = _cn ** 2 * m.alpha_1 + _cm ** 2 * m.alpha_2
        op.text("alpha_y = n ^ 2 * alpha_1 + m ^ 2 * alpha_2 = {} 1/C;".format(alpha_y))
        alpha_xy = 2 * (m.alpha_1 - m.alpha_2) * _cm * _cn
        op.text("alpha_xy = 2(alpha_1 - alpha_2)mn = {} 1/C;".format(alpha_xy))
        op.text("")
        op.text("delta_T = {} C;".format(c.delta_thm))
        op.text("")
    else:
        alpha_x = alpha_y = alpha_xy = 0

    if c.delta_moi:
        op.part("moisture conditions")
        op.text("beta_1 = {} 1/%M;".format(m.beta_1))
        op.text("beta_2 = {} 1/%M;".format(m.beta_2))
        op.text("")
        beta_x = _cm ** 2 * m.beta_1 + _cn ** 2 * m.beta_2
        op.text("beta_x = m ^ 2 * beta_1 + n ^ 2 * beta_2 = {} 1/%M;".format(beta_x))
        beta_y = _cn ** 2 * m.beta_1 + _cm ** 2 * m.beta_2
        op.text("beta_y = n ^ 2 * beta_1 + m ^ 2 * beta_2 = {} 1/%M;".format(beta_y))
        beta_xy = 2 * (m.beta_1 - m.beta_2) * _cm * _cn
        op.text("beta_xy = 2(beta_1 - beta_2)mn = {} 1/%M;".format(beta_xy))
        op.text("")
        op.text("delta_M = {} %M;".format(c.delta_moi))
        op.text("")
    else:
        beta_x = beta_y = beta_xy = 0

    op.part("reduced S_matrix (1/Pa)")
    rS_matrix = m.get_reduced_S_matrix()
    op.matrix(rS_matrix)
    op.text("")
    op.part("reduced Q_matrix (Pa)")
    rQ_matrix = m.get_reduced_Q_matrix()
    op.matrix(rQ_matrix)
    op.text("")

    op.part("transformation matrix")
    T_matrix = np.matrix([[ _cm ** 2,      _cn ** 2,     2 * _cm * _cn     ],
                          [ _cn ** 2,      _cm ** 2,    -2 * _cm * _cn     ],
                          [-1 * _cm * _cn, _cm * _cn,   _cm ** 2 - _cn ** 2],])
    op.matrix(T_matrix)
    op.text("")
    op.part("inversed transformation matrix")
    Tr_matrix = np.linalg.inv(T_matrix)
    op.matrix(Tr_matrix)
    op.text("")

    op.part("transformed reduced S_matrix (1/Pa)")
    trS_matrix = m.get_trans_reduced_S_matrix(rotate_deg = c.rotate_deg)
    op.matrix(trS_matrix)
    op.text("")
    op.part("transformed reduced Q_matrix (Pa)")
    trQ_matrix = m.get_trans_reduced_Q_matrix(rotate_deg = c.rotate_deg)
    op.matrix(trQ_matrix)
    op.text("")

    op.section("solution")
    # solution 1.
    op.part("first solve stress vs load")
    stress_1 = Matrix([[c.x_normal_load / (c.width     * c.thickness)],
                       [c.y_normal_load / (c.thickness * c.length   )],
                       [c.xy_shear_load / (c.length    * c.width    )],])
    
    stress_2 = Matrix([[c.x_normal_stress],
                       [c.y_normal_stress],
                       [c.xy_shear_stress],])
    
    solution = solve(stress_1 - stress_2)
    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    # x_nr = c.x_normal_stress if isinstance(c.x_normal_stress, (int, float)) else solution[c.x_normal_stress]
    # y_nr = c.y_normal_stress if isinstance(c.y_normal_stress, (int, float)) else solution[c.y_normal_stress]

    x_nr = solution[c.x_normal_stress] if c.x_normal_stress in solution else c.x_normal_stress
    y_nr = solution[c.y_normal_stress] if c.y_normal_stress in solution else c.y_normal_stress
    
    op.part("normal stresses")
    op.text("sigma_x = {} Pa;".format(x_nr), end=" ")
    op.text("sigma_y = {} Pa;".format(y_nr))
    op.text("")

    # xy_sr = c.xy_shear_stress if isinstance(c.xy_shear_stress, (int, float)) else solution[c.xy_shear_stress]

    xy_sr = solution[c.xy_shear_stress] if c.xy_shear_stress in solution else c.xy_shear_stress

    op.part("shear stresses")
    op.text("tau_xy = {} Pa;".format(xy_sr))
    op.text("")

    # solution 2.
    op.part("then solved strain vs stress")
    strain = Matrix([[c.x_normal_strain - alpha_x * c.delta_thm - beta_x * c.delta_moi  ],
                     [c.y_normal_strain - alpha_y * c.delta_thm - beta_y * c.delta_moi  ],
                     [c.xy_shear_strain - alpha_xy * c.delta_thm - beta_xy * c.delta_moi],])

    stress = Matrix([[x_nr],
                     [y_nr],
                     [xy_sr],])

    solution = solve(Matrix(trS_matrix) * stress - strain)

    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    x_ns = solution[c.x_normal_strain] if c.x_normal_strain in solution  else c.x_normal_strain 
    y_ns = solution[c.y_normal_strain] if c.y_normal_strain in solution  else c.y_normal_strain 
    
    op.part("normal strains")
    op.text("epsilon_x = {};".format(x_ns), end=" ")
    op.text("epsilon_y = {};".format(y_ns))
    op.text("")

    xy_ss = solution[c.xy_shear_strain] if c.xy_shear_strain in solution else c.xy_shear_strain 

    op.part("shear strains")
    op.text("gamma_xy = {};".format(xy_ss))
    op.text("")

    # solution 3 (optional)
    op.part("thirdly solve special z strains (ref. chapter conti. 5, P6)")
    op.text("epsilon_z = alpha_3 * delta_T + beta_3 * delta_M + (S_13 * m ^ 2 + S_23 * n ^ 2)  * sigma_x + (S_13 * n ^ 2 + S_23 * m ^ 2) * sigma_y + 2(S_13 - S_23) * mn * tau_xy")
    op.text("")
    if c.delta_thm:
        op.part("thermal conditions")
        op.text("alpha_3 = {} 1/C;".format(m.alpha_3), end=" ")
        op.text("delta_T = {} C;".format(c.delta_thm))
        op.text("")
    if c.delta_moi:
        op.part("moisture conditions")
        op.text("beta_3 = {} 1/%M;".format(m.beta_3), end=" ")
        op.text("delta_M = {} %M;".format(c.delta_moi))
        op.text("")

    S_matrix = m.get_S_matrix()
    S_13 = S_matrix[0, 2]
    S_23 = S_matrix[1, 2]
    op.part("partical S_matrix")
    op.text("S_13 = {} 1/Pa;".format(S_13))
    op.text("S_23 = {} 1/Pa;".format(S_23))
    op.text("")

    op.part("z strains")
    x_nr = solution[x_nr] if x_nr in solution else x_nr
    y_nr = solution[y_nr] if y_nr in solution else y_nr
    xy_sr = solution[xy_sr] if xy_sr in solution else xy_sr

    epsilon_z = m.alpha_3 * c.delta_thm + m.beta_3 * c.delta_moi + (S_13 * _cm ** 2 + S_23 * _cn ** 2)  * x_nr + (S_13 * _cn ** 2 + S_23 * _cm ** 2) * y_nr + 2 * (S_13 - S_23) * _cm * _cn * xy_sr
    op.text("epsilon_z = {};".format(epsilon_z))
    op.text("")

    # solution 4.
    f_l = c.length * (1 + x_ns)
    f_w = c.width * (1 + y_ns)
    f_t = c.thickness * (1 + epsilon_z)
    f_xy = 90 + xy_ss * 180 / math.pi
    
    op.title("final size")
    op.section("final length")
    op.text("length    = {} m;".format(f_l))
    op.text("widht     = {} m;".format(f_w))
    op.text("")
    op.section("final angle")
    op.text("xy angle = {} deg;".format(f_xy))
    op.text("")
    op.section("final z")
    op.text("thickness = {} m;".format(f_t))
