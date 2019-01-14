#!/usr/bin/env python3
# Loading and paraphrasing fiber/matrix engineering properties.
#
# Version 1.0.0
# 2019.01.13
#
# Author: Liu Jia
#


import sys
import os
import re
from slib.pstrlst import del_comment, del_stri_ht_space
import numpy as np
import math


class Material(object):
    def __init__(self, material_file):
        if not os.path.isfile(material_file):
            raise ValueError("Material property file doesn't exist: {}.".format(material_file))
        
        self.E_1 = None
        self.E_2 = None
        self.E_3 = None

        self.nu_23 = None
        self.nu_13 = None
        self.nu_12 = None

        self.G_23 = None
        self.G_13 = None
        self.G_12 = None

        self.alpha_1 = None
        self.alpha_2 = None
        self.alpha_3 = None

        self.beta_1 = None
        self.beta_2 = None
        self.beta_3 = None

        # Loading property file.
        props = open(material_file)
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
            if variable == "E_1":
                self.E_1 = self._N_or_n(value)
            elif variable == "E_2":
                self.E_2 = self._N_or_n(value)
            elif variable == "E_3":
                self.E_3 = self._N_or_n(value)

            elif variable == "nu_23":
                self.nu_23 = self._N_or_n(value)
            elif variable == "nu_13":
                self.nu_13 = self._N_or_n(value)
            elif variable == "nu_12":
                self.nu_12 = self._N_or_n(value)

            elif variable == "G_23":
                self.G_23 = self._N_or_n(value)
            elif variable == "G_13":
                self.G_13 = self._N_or_n(value)
            elif variable == "G_12":
                self.G_12 = self._N_or_n(value)
            
            elif variable == "alpha_1":
                self.alpha_1 = self._N_or_n(value)
            elif variable == "alpha_2":
                self.alpha_2 = self._N_or_n(value)
            elif variable == "alpha_3":
                self.alpha_3 = self._N_or_n(value)
            
            elif variable == "beta_1":
                self.beta_1 = self._N_or_n(value)
            elif variable == "beta_2":
                self.beta_2 = self._N_or_n(value)
            elif variable == "beta_3":
                self.beta_3 = self._N_or_n(value)

            else:
                print("Warnning: unknown variable: {}.".format(variable))
                continue
            
        props.close()

    # Return None or number.
    def _N_or_n(self, stri):
        if stri[0].lower() in ("n", "u"):
            return None
        else:
            return float(stri)
    
    def get_S_matrix(self):
        S_11 = 1 / self.E_1
        S_12 = -1 * self.nu_12 / self.E_1
        S_13 = -1 * self.nu_13 / self.E_1
        
        S_22 = 1 / self.E_2
        S_23 = -1 * self.nu_23 / self.E_2
        
        S_33 = 1 / self.E_3

        S_44 = 1 / self.G_23
        S_55 = 1 / self.G_13
        S_66 = 1 / self.G_12

        S_matrix = np.matrix([[S_11, S_12, S_13, 0,    0,    0   ],
                              [S_12, S_22, S_23, 0,    0,    0   ],
                              [S_13, S_23, S_33, 0,    0,    0   ],
                              [0,    0,    0,    S_44, 0,    0   ],
                              [0,    0,    0,    0,    S_55, 0   ],
                              [0,    0,    0,    0,    0,    S_66],])
        return S_matrix
    
    def get_reduced_S_matrix(self):
        S_11 = 1 / self.E_1
        S_12 = -1 * self.nu_12 / self.E_1
        
        S_22 = 1 / self.E_2
        
        S_66 = 1 / self.G_12

        rS_matrix = np.matrix([[S_11, S_12, 0   ],
                               [S_12, S_22, 0   ],
                               [0,    0,    S_66],])
        return rS_matrix
    
    # def get_trans_reduced_S_matrix(self, rotate_deg=0):
    #     S_11 = 1 / self.E_1
    #     S_12 = -1 * self.nu_12 / self.E_1
    #     
    #     S_22 = 1 / self.E_2
    #     
    #     S_66 = 1 / self.G_12
    # 
    #     m = math.cos(rotate_deg * math.pi / 180)
    #     n = math.sin(rotate_deg * math.pi / 180)
    # 
    #     bS_11 = S_11 * m ** 4 + (2 * S_12 + S_66) * m ** 2 * n ** 2 + S_22 * n ** 4
    #     bs_12 = (S_11 + S_22 - S_66) * m ** 2 * n ** 2 + S_12 * (m ** 4 + n ** 4)
    #     bs_16 = (2 * S_11 - 2 * S_12 - S_66) * m ** 3 * n - (2 * S_22 - 2 * S_12 - S_66) * m * n ** 3
    #     bs_22 = S_11 * n ** 4 + (2 * S_12 + S_66) * m ** 2 * n ** 2 + S_22 * m ** 4
    #     bs_26 = (2 * S_11 - 2 * S_12 - S_66) * m * n ** 3 - (2 * S_22 - 2 * S_12 - S_66) * m ** 3 * n
    #     bs_66 = 2 * (2 * S_11 + 2 * S_22 - 4 * S_12 - S_66) * m ** 2 * n ** 2 + S_66 * (m ** 4 + n ** 4)
    # 
    #     trS_matrix = np.matrix([[bS_11, bs_12, bs_16],
    #                             [bs_12, bs_22, bs_26],
    #                             [bs_16, bs_26, bs_66],])
    #     return trS_matrix

    def get_trans_reduced_S_matrix(self, rotate_deg=0):
        S_11 = 1 / self.E_1
        S_12 = -1 * self.nu_12 / self.E_1
        
        S_22 = 1 / self.E_2
        
        S_66 = 1 / self.G_12

        m = math.cos(rotate_deg * math.pi / 180)
        n = math.sin(rotate_deg * math.pi / 180)

        T_matrix = np.matrix([[ m ** 2,    n ** 2,  2 * m * n      ],
                              [ n ** 2,    m ** 2, -2 * m * n      ],
                              [-1 * m * n, m * n,   m ** 2 - n ** 2],])
        
        Tr_matrix = np.linalg.inv(T_matrix)

        rS_matrix = np.matrix([[S_11, S_12, 0           ],
                              [S_12, S_22, 0           ],
                              [0,    0,    1 / 2 * S_66],])
        
        twice_matrix = np.matrix([[1, 0, 0],
                                  [0, 1, 0],
                                  [0, 0, 2],])
        
        trS_matrix = twice_matrix.dot(Tr_matrix).dot(rS_matrix).dot(T_matrix)

        return trS_matrix
    
    def get_ABD_matrix(self, thickness, rotate_deg):
        def _Aij(i, j):
            total = 0
            for k in range(len(rotate_deg)):
                t = thickness[k + 1] - thickness[k]
                r = rotate_deg[k]
                trQ_matrix = self.get_trans_reduced_Q_matrix(rotate_deg=r)
                total = total + trQ_matrix[i - 1, j - 1] * t
            return total
        
        def _Bij(i, j):
            total = 0
            for k in range(len(rotate_deg)):
                t = thickness[k + 1] ** 2 - thickness[k] ** 2
                r = rotate_deg[k]
                trQ_matrix = self.get_trans_reduced_Q_matrix(rotate_deg=r)
                total = total + 1 / 2 * trQ_matrix[i - 1, j - 1] * t
            return total
        
        def _Dij(i, j):
            total = 0
            for k in range(len(rotate_deg)):
                t = thickness[k + 1] ** 3 - thickness[k] ** 3
                r = rotate_deg[k]
                trQ_matrix = self.get_trans_reduced_Q_matrix(rotate_deg=r)
                total = total + 1 / 3 * trQ_matrix[i - 1, j - 1] * t
            return total

        ABD_matrix = np.matrix([[_Aij(1, 1), _Aij(1, 2), _Aij(1, 3), _Bij(1, 1), _Bij(1, 2), _Bij(1, 3)],
                                [_Aij(2, 1), _Aij(2, 2), _Aij(2, 3), _Bij(2, 1), _Bij(2, 2), _Bij(2, 3)],
                                [_Aij(3, 1), _Aij(3, 2), _Aij(3, 3), _Bij(3, 1), _Bij(3, 2), _Bij(3, 3)],
                                [_Bij(1, 1), _Bij(1, 2), _Bij(1, 3), _Dij(1, 1), _Dij(1, 2), _Dij(1, 3)],
                                [_Bij(2, 1), _Bij(2, 2), _Bij(2, 3), _Dij(2, 1), _Dij(2, 2), _Dij(2, 3)],
                                [_Bij(3, 1), _Bij(3, 2), _Bij(3, 3), _Dij(3, 1), _Dij(3, 2), _Dij(3, 3)],])
        
        return ABD_matrix

    def get_C_matrix(self):
         S_matrix = self.get_S_matrix()
         C_matrix = np.linalg.inv(S_matrix)
         return C_matrix

    def get_reduced_Q_matrix(self):
        rS_matrix = self.get_reduced_S_matrix()
        rQ_matrix = np.linalg.inv(rS_matrix)
        return rQ_matrix
    
    def get_trans_reduced_Q_matrix(self, rotate_deg=0):
        trS_matrix = self.get_trans_reduced_S_matrix(rotate_deg)
        trQ_matrix = np.linalg.inv(trS_matrix)
        return trQ_matrix
    
    def get_abd_matrix(self, thickness, rotate_deg):
        ABD_matrix = self.get_ABD_matrix(thickness, rotate_deg)
        abd_matrix = np.linalg.inv(ABD_matrix)
        return abd_matrix
    