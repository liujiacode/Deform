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


def solve_ABD_matrix(c, m, op):
    op.title("solving ABD matrix")

    op.section("knowns and unknowns")
    op.part("stresses")
    op.text("sigma_x = {} Pa;".format(c.x_normal_stress), end=" ")
    op.text("sigma_y = {} Pa;".format(c.y_normal_stress), end=" ")
    op.text("tau_xy = {} Pa;".format(c.xy_shear_stress))
    op.text("")
    op.part("normal load resultants")
    op.text("N_x = {} N/m;".format(c.x_normal_load), end=" ")
    op.text("N_y = {} N/m;".format(c.y_normal_load), end=" ")
    op.text("N_xy = {} N/m;".format(c.xy_shear_load))
    op.text("")
    op.part("bending moment resultants")
    op.text("M_x = {} Nm/m;".format(c.x_bending_moment), end=" ")
    op.text("M_y = {} Nm/m;".format(c.y_bending_moment), end=" ")
    op.text("M_xy = {} Nm/m;".format(c.xy_twisting_moment))
    op.text("")
    op.part("strains")
    op.text("epsilon_x = {} Pa;".format(c.x_normal_strain), end=" ")
    op.text("epsilon_y = {} Pa;".format(c.y_normal_strain), end=" ")
    op.text("gamma_xy = {} Pa;".format(c.xy_shear_strain))
    op.text("")
    op.part("curvatures")
    op.text("kappa_x = {} 1/m;".format(c.x_surface_curvature), end=" ")
    op.text("kappa_y = {} 1/m;".format(c.y_surface_curvature), end=" ")
    op.text("kappa_xy = {} 1/m;".format(c.xy_twisting_curvature))
    op.text("")

    op.part("rotate degree")
    op.text("theta = {} deg;".format(c.rotate_deg))
    op.text("")

    op.part("ABD_matrix (**)")
    ABD_matrix = m.get_ABD_matrix(c.thickness, c.rotate_deg)
    op.matrix(ABD_matrix)
    op.text("")
    op.part("abd_matrix (**)")
    abd_matrix = m.get_abd_matrix(c.thickness, c.rotate_deg)
    op.matrix(abd_matrix)
    op.text("")

    op.section("solution")
    # solution 1.
    op.part("first solve load vs stress")
    integ_z1 = c.thickness[-1] - c.thickness[0]
    integ_z2 = 1 / 2 * (c.thickness[-1] ** 2 - c.thickness[0] ** 2 )
    op.text("integrate1(z) = {}.".format(integ_z1))
    op.text("integrate2(z) = {}.".format(integ_z2))
    op.text("")

    load_1 = Matrix([[c.x_normal_stress * (integ_z1)],
                     [c.y_normal_stress * (integ_z1)],
                     [c.xy_shear_stress * (integ_z1)],
                     [c.x_normal_stress * (integ_z2)],
                     [c.y_normal_stress * (integ_z2)],
                     [c.xy_shear_stress * (integ_z2)],])
    
    load_2 = Matrix([[c.x_normal_load     ],
                     [c.y_normal_load     ],
                     [c.xy_shear_load     ],
                     [c.x_bending_moment  ],
                     [c.y_bending_moment  ],
                     [c.xy_twisting_moment],])
    
    solution = solve(load_1 - load_2)
    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    # x_nl = c.x_normal_load if isinstance(c.x_normal_load, (int, float)) else solution[c.x_normal_load]
    # y_nl = c.y_normal_load if isinstance(c.y_normal_load, (int, float)) else solution[c.y_normal_load]
    # xy_sl = c.xy_shear_load if isinstance(c.xy_shear_load, (int, float)) else solution[c.xy_shear_load]

    x_nl = solution[c.x_normal_load] if c.x_normal_load in solution else c.x_normal_load
    y_nl = solution[c.y_normal_load] if c.y_normal_load in solution else c.y_normal_load
    xy_sl = solution[c.xy_shear_load] if c.xy_shear_load in solution else c.xy_shear_load
    
    op.part("normal load resultants")
    op.text("N_x = {} N/m;".format(x_nl), end=" ")
    op.text("N_y = {} N/m;".format(y_nl), end=" ")
    op.text("N_xy = {} N/m;".format(xy_sl))
    op.text("")

    # x_bm = c.x_bending_moment if isinstance(c.x_bending_moment, (int, float)) else solution[c.x_bending_moment]
    # y_bm = c.y_bending_moment if isinstance(c.y_bending_moment, (int, float)) else solution[c.y_bending_moment]
    # xy_tm = c.xy_twisting_moment if isinstance(c.xy_twisting_moment, (int, float)) else solution[c.xy_twisting_moment]

    x_bm = solution[c.x_bending_moment] if c.x_bending_moment in solution else c.x_bending_moment
    y_bm = solution[c.y_bending_moment] if c.y_bending_moment in solution else c.y_bending_moment
    xy_tm = solution[c.xy_twisting_moment] if c.xy_twisting_moment in solution else c.xy_twisting_moment
    
    op.part("bending moment resultants")
    op.text("M_x = {} Nm/m;".format(x_bm), end=" ")
    op.text("M_y = {} Nm/m;".format(y_bm), end=" ")
    op.text("M_xy = {} Nm/m;".format(xy_tm))
    op.text("")

    # solution 2.
    op.part("then solved strain vs load")
    strain = Matrix([[c.x_normal_strain      ],
                     [c.y_normal_strain      ],
                     [c.xy_shear_strain      ],
                     [c.x_surface_curvature  ],
                     [c.y_surface_curvature  ],
                     [c.xy_twisting_curvature],])
    
    load = Matrix([[x_nl ],
                   [y_nl ],
                   [xy_sl],
                   [x_bm ],
                   [y_bm ],
                   [xy_tm],])

    solution = solve(Matrix(ABD_matrix) * strain - load)

    for i in solution:
        op.text("{} = {};".format(i, solution[i]))
    op.text("")

    x_ns = solution[c.x_normal_strain] if c.x_normal_strain in solution else c.x_normal_strain
    y_ns = solution[c.y_normal_strain] if c.y_normal_strain in solution else c.y_normal_strain
    xy_ss = solution[c.xy_shear_strain] if c.xy_shear_strain in solution else c.xy_shear_strain

    op.part("strains")
    op.text("epsilon_x = {} Pa;".format(x_ns), end=" ")
    op.text("epsilon_y = {} Pa;".format(y_ns), end=" ")
    op.text("gamma_xy = {} Pa;".format(xy_ss))
    op.text("")

    x_sc = solution[c.x_surface_curvature] if c.x_surface_curvature in solution else c.x_surface_curvature
    y_sc = solution[c.y_surface_curvature] if c.y_surface_curvature in solution else c.y_surface_curvature
    xy_tc = solution[c.xy_twisting_curvature] if c.xy_twisting_curvature in solution else c.xy_twisting_curvature

    op.part("curvatures")
    op.text("kappa_x = {} 1/m;".format(x_sc), end=" ")
    op.text("kappa_y = {} 1/m;".format(y_sc), end=" ")
    op.text("kappa_xy = {} 1/m;".format(xy_tc))
    op.text("")

    # solution 3.
    epsilon_xs = []
    epsilon_ys = []
    gamma_xys = []
    upper_sigma_xs = []
    upper_sigma_ys = []
    upper_tau_xys = []
    lower_sigma_xs = []
    lower_sigma_ys = []
    lower_tau_xys = []

    op.title("off axis strains solution")
    for i in range(len(c.thickness)):
        op.section("scale {}".format(i + 1))
        
        z_scale = c.thickness[i]
        op.text("z scale = {}.".format(z_scale))
        op.text("")
        
        epsilon_x = x_ns + z_scale * x_sc
        epsilon_y = y_ns + z_scale * y_sc
        gamma_xy = xy_ss + z_scale * xy_tc

        epsilon_xs.append(epsilon_x)
        epsilon_ys.append(epsilon_y)
        gamma_xys.append(gamma_xy)

        op.part("off strain")
        op.text("off epsilon_x = {}.".format(epsilon_x))
        op.text("off epsilon_y = {}.".format(epsilon_y))
        op.text("off gamma_xy = {}.".format(gamma_xy))
        op.text("")

        if i > 0:
            theta = c.rotate_deg[i - 1]
            upper_trQ = m.get_trans_reduced_Q_matrix(rotate_deg=theta)
            stress = upper_trQ.dot(np.matrix([[epsilon_x], [epsilon_y], [gamma_xy],]))
            
            op.draft("upper layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(upper_trQ)
            op.text("[stress] = [trQ][strain]")
            op.draft("end")
            op.text("")

            upper_sigma_xs.append(stress[0, 0])
            upper_sigma_ys.append(stress[1, 0])
            upper_tau_xys.append(stress[2, 0])
            
            op.part("off upper layer stress")
            op.text("theta = {}.".format(theta))
            op.text("off sigma_x = {}.".format(stress[0, 0]))
            op.text("off sigma_y = {}.".format(stress[1, 0]))
            op.text("off tau_xy = {}.".format(stress[2, 0]))
            op.text("")
        if i < len(c.thickness) - 1:
            theta = c.rotate_deg[i]
            lower_trQ = m.get_trans_reduced_Q_matrix(rotate_deg=theta)
            stress = lower_trQ.dot(np.matrix([[epsilon_x], [epsilon_y], [gamma_xy],]))
            
            op.draft("lower layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(lower_trQ)
            op.text("[stress] = [trQ][strain]")
            op.draft("end")
            op.text("")

            lower_sigma_xs.append(stress[0, 0])
            lower_sigma_ys.append(stress[1, 0])
            lower_tau_xys.append(stress[2, 0])
            
            op.part("off lower layer stress")
            op.text("theta = {}.".format(theta))
            op.text("off sigma_x = {}.".format(stress[0, 0]))
            op.text("off sigma_y = {}.".format(stress[1, 0]))
            op.text("off tau_xy = {}.".format(stress[2, 0]))
            op.text("")
    
    op.title("on axis strains solution")
    for i in range(len(c.thickness)):
        op.section("scale {}".format(i + 1))
        
        z_scale = c.thickness[i]
        op.text("z scale = {}.".format(z_scale))
        op.text("")
        
        if i > 0:
            theta = c.rotate_deg[i - 1]
            m = math.cos(theta * math.pi / 180)
            n = math.sin(theta * math.pi / 180)
            T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                  [ n ** 2,    m ** 2, -2 * m * n      ],
                                  [-1 * m * n, m * n,   m ** 2 - n ** 2],])
            
            op.part("on upper layer strain")
            op.text("theta = {}.".format(theta))
            strain = T_matrix.dot(np.matrix([[epsilon_xs[i]], [epsilon_ys[i]], [1 / 2 * gamma_xys[i]]]))
            
            op.draft("upper layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(T_matrix)
            op.text("[1/2 strain] = [T][1/2 strain]")
            op.draft("end")
            op.text("")
            
            op.text("on epsilon_x = {}.".format(strain[0, 0]))
            op.text("on epsilon_y = {}.".format(strain[1, 0]))
            op.text("on gamma_xy = {}.".format(strain[2, 0] * 2))
            op.text("")
        
        if i < len(c.thickness) - 1:
            theta = c.rotate_deg[i]
            m = math.cos(theta * math.pi / 180)
            n = math.sin(theta * math.pi / 180)
            T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                  [ n ** 2,    m ** 2, -2 * m * n      ],
                                  [-1 * m * n, m * n,   m ** 2 - n ** 2],])
            
            op.part("on lower layer strain")
            op.text("theta = {}.".format(theta))
            strain = T_matrix.dot(np.matrix([[epsilon_xs[i]], [epsilon_ys[i]], [1 / 2 * gamma_xys[i]]]))

            op.draft("lower layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(T_matrix)
            op.text("[1/2 strain] = [T][1/2 strain]")
            op.draft("end")
            op.text("")

            op.text("on epsilon_x = {}.".format(strain[0, 0]))
            op.text("on epsilon_y = {}.".format(strain[1, 0]))
            op.text("on gamma_xy = {}.".format(strain[2, 0] * 2))
            op.text("")
        
        if i > 0:
            theta = c.rotate_deg[i - 1]
            m = math.cos(theta * math.pi / 180)
            n = math.sin(theta * math.pi / 180)
            T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                  [ n ** 2,    m ** 2, -2 * m * n      ],
                                  [-1 * m * n, m * n,   m ** 2 - n ** 2],])
            
            op.part("on upper layer stress")
            op.text("theta = {}.".format(theta))
            stress = T_matrix.dot(np.matrix([[upper_sigma_xs[i - 1]], [upper_sigma_ys[i - 1]], [upper_tau_xys[i - 1]]]))
            
            op.draft("upper layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(T_matrix)
            op.text("[stress] = [T][stress]")
            op.draft("end")
            op.text("")

            op.text("on sigma_x = {}.".format(stress[0, 0]))
            op.text("on sigma_y = {}.".format(stress[1, 0]))
            op.text("on tau_xy = {}.".format(stress[2, 0]))
            op.text("")
        
        if i < len(c.thickness) - 1:
            theta = c.rotate_deg[i]
            m = math.cos(theta * math.pi / 180)
            n = math.sin(theta * math.pi / 180)
            T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                                  [ n ** 2,    m ** 2, -2 * m * n      ],
                                  [-1 * m * n, m * n,   m ** 2 - n ** 2],])
            
            op.part("on lower layer stress")
            op.text("theta = {}.".format(theta))
            stress = T_matrix.dot(np.matrix([[upper_sigma_xs[i]], [upper_sigma_ys[i]], [upper_tau_xys[i]]]))

            op.draft("lower layer info")
            op.text("theta = {}.".format(theta))
            op.matrix(T_matrix)
            op.text("[stress] = [T][stress]")
            op.draft("end")
            op.text("")

            op.text("on sigma_x = {}.".format(stress[0, 0]))
            op.text("on sigma_y = {}.".format(stress[1, 0]))
            op.text("on tau_xy = {}.".format(stress[2, 0]))
            op.text("")
            