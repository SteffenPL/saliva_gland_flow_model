#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:43:17 2018

@author: plunder
"""
import sympy as sp
import numpy as np

class Units:
    m = 1  # meter
    μm = 1e-6
    nm = 1e-12
    mM = 1  # mili molar (mM = mol/m³ = 1e-3 mol/L)
    J = 1  # joule
    nM = 1e-6  # nano molar
    K = 1   # kelvin
    s = 1   # seconds
    μM = 1e-3  # micro molar
    g = 1  # gram
    mol = 1
    C = 1  # °C... hard to translate... use K-C = constant
    Simens = 1


class shelley_model:

    def __init__(self):

        s = self

        # Init Parameters
        s.na_B = 140.7 * Units.mM
        s.cl_B = 125. * Units.mM
        s.k_B = 5 * Units.mM
        s.pH = 7.35
        s.hco3_B = 42.9 * Units.mM
        s.co2_B = 1.28 * Units.mM
        s.h_B = 10**(-s.pH)

        # TODO: adjust diameter according to speed

        s.Area = 20.9 * Units.nm**2
        s.V_W = 18.05e-6 * Units.m**3/Units.mol
        s.chi_C = 60 * Units.mM
        s.Psi = 0.2 * Units.mM
        s.R = 8.3145 * Units.J / (Units.g * Units.K)
        s.T = 310 * Units.K
        s.F = 96.485 * Units.C / Units.mol

        s.ca_C = 60 * Units.μM
        s.w_A = 2.090 * Units.μm**3  # todo: find if this really is s.w_A

        s.r_NKA = 1.305e-3 / (Units.mM**3 * Units.s)
        s.beta_NKA = 6.47e-4 / Units.mM

        s.gamma_1 = 157.55 / Units.s**2
        s.gamma_2 = 2e-5 / (Units.mM**4 * Units.s**2)
        s.gamma_3 = 1.03 / Units.s
        s.gamma_4 = 1.38e-6 / (Units.mM**4*Units.s)

        s.k1 = 1.4159e3/Units.s
        s.k2 = 2.5296e9/Units.s
        s.k3 = 5.8599/Units.s
        s.k4 = 9.3124e7/Units.s
        s.k5 = -6.039e-1/(Units.mM*Units.s)
        s.k6 = 1.0004/Units.s

        s.k_1 = 1.4294e11/Units.s
        s.k_2 = 1.7857e2/Units.s
        s.k_3 = 1.0652e8/Units.s
        s.k_4 = 5.1418/Units.s
        s.k_5 = 1e8/(Units.mM*Units.s)
        s.k_6 = -1.9352e-1/(Units.mM*Units.s)

        uCond = Units.Simens/Units.m**2

        s.G_CFTR = 90 * uCond
        s.alpha_CFTR = 1 * Units.m**(-2)

        s.G_ENaC = 54.6 * uCond
        s.alpha_ENaC = 1 * Units.m**(-2)

        s.G_BK = 10 * uCond
        # alpha_BK = 1

        s.alpha_NKA = 5e-8 * Units.m**(-2)
        s.alpha_NHE_A = 0.1 * Units.m**(-2)
        s.alpha_NHE_B = 1 * Units.m**(-2)
        s.alpha_AE2_A = 1 * Units.m**(-2)
        s.alpha_AE2_B = 5 * Units.m**(-2)
        s.alpha_NBC_A = 0 * Units.m**(-2)
        s.alpha_NBC_B = 1 * Units.m**(-2)

        s.alpha_NKCC1 = 1e-9 * Units.m**(-2)  # taken from transport.m line 264

        s.G_na_P = 0.4 * uCond
        s.G_cl_P = 1 * uCond
        s.G_k_P = 3 * uCond  # this one is maybe 5?
        s.p_co2 = 5e-6 * 1./Units.s


        # TODO find exact values!!!

        s.L_A = 1e-6 * Units.m/Units.s
        s.L_B = 1e-6 * Units.m/Units.s
        s.eta = 1.49
        s.K_m_BK = 0.26e-3 * Units.mM
        s.G_k_B = 0. * uCond  # float('inf')
        s.alpha_k_B = 0.  # float('inf')


    def ode(self,x):
        """returns np.ndarray (with shape=(13,))

            parameters arrangement:
                x = [w_C,
                na_A, cl_A, k_A, h_A, hco3_A, co2_A,
                na_C, cl_C, k_C, h_C, hco3_C, co2_C]

                where w_C denotes the volume of the duct cell,
                *_A denotes the concentration on the apical side
                and *_C denotes the concentration in the duct cells



        .. note::
        .. seealso:: :class:`MainClass2`
        .. warning:: all concentrations must be positive (since log(*_{A,C]}) is may computed in the model.)
        .. todo:: transform parameters and internal state variables as class parameters
        """

        # to shorten the notations, we set
        s = self

        # unkowns
        w_C, na_A, cl_A, k_A, h_A, hco3_A, co2_A, na_C, cl_C, k_C, h_C, hco3_C, co2_C = x

        con_A = [na_A, cl_A, k_A, h_A, hco3_A, co2_A]
        con_B = [s.na_B, s.cl_B, s.k_B, s.h_B, s.hco3_B]
        con_C = [na_C, cl_C, k_C, h_C, hco3_C, co2_C]


        z_cl = -1
        z_hco3 = -1
        z_na = +1
        z_k = +1

        s.J_w_A = s.L_A*s.V_W*(sum(con_A)-sum(con_C) + (-s.chi_C/w_C + s.Psi))
        s.J_w_B = s.L_B*s.V_W*(sum(con_C)-sum(con_B) + (s.chi_C/w_C))
        s.rhs_w_c = s.Area*(s.J_w_B - s.J_w_A)

        # to resolve the algebraic constraint, we need to
        # compute all potentials first, then find the
        # current volumes and after that compute the currents

        # nerst equation
        s.V_cl_A = s.R*s.T/(z_cl*s.F) * np.log(cl_A/cl_C)


        # CFTR channel transports bicarbonate
        s.V_hco3_A = s.R*s.T/(z_hco3*s.F)*np.log(hco3_A / hco3_C)

        # sodium channels
        s.V_na_A = s.R*s.T/(z_na*s.F)*np.log(na_A/na_C)

        # potassium channels

        # open prob via Hill Function, calcium activated
        ca_C_eta = s.ca_C**s.eta
        s.P_BK = ca_C_eta /(ca_C_eta  + s.K_m_BK**s.eta)

        s.V_k_A = s.R*s.T/(z_k*s.F) * np.log(k_A/k_C)

        # on the basolateral side, Is.k1 channel
        s.V_k_B = s.R*s.T/(z_k*s.F)*np.log(s.k_B/k_C)

        # paracelluar potentials

        def paracelluar_potentials(G_n_P,n_A,n_C,n_z):
            return s.R*s.T/(n_z*s.F)*np.log(n_A/n_C)

        def paracelluar_currents(G_n_P, n_A, n_C, V_n_P):
            I_n_P = G_n_P * s.Area * (V_t - V_n_P)
            return I_n_P

        s.V_na_P = paracelluar_potentials(s.G_na_P, na_A, na_C, z_na)
        s.V_cl_P = paracelluar_potentials(s.G_cl_P, cl_A, cl_C, z_cl)
        s.V_k_P = paracelluar_potentials(s.G_k_P, k_A, k_C, z_k)


        # Na-K_ATPase (NKA)
        s.u_NKA_B = s.r_NKA*(s.k_B**2*na_C**3)/(s.k_B**2 + s.beta_NKA*na_C**3)
        s.J_NKA_B = s.alpha_NKA*s.u_NKA_B

        # membrane potential
        # algebraic constraint!
        # rhs_kirchhoff_A = s.I_ENaC+s.I_na_P+s.I_CFTR+s.I_cl_P+s.I_BK+s.I_k_P+s.I_CFTR_B
        # rhs_kirchhoff_B = -s.F*w_C*s.J_NKA_B-s.I_na_P-s.I_cl_P+s.I_k_B-s.I_k_P

        # compute the matrix s.t. B*[s.V_A,s.V_B] == b equals the constraint

        B = np.matrix([[0.,0.],[0.,0.]])


        B[0,0] = s.Area*(s.G_CFTR*s.alpha_CFTR
                 +1/4*s.G_CFTR*s.alpha_CFTR
                 +    s.G_ENaC*s.alpha_ENaC
                 +            s.P_BK*s.G_BK
                 +(s.G_na_P+s.G_cl_P+s.G_k_P))

        B[0,1] = -s.Area*(s.G_na_P+s.G_cl_P+s.G_k_P)
        B[1,0] = -s.Area*(s.G_na_P+s.G_cl_P+s.G_k_P)
        B[1,1] = s.Area*(s.G_k_B*s.alpha_k_B +(s.G_na_P+s.G_cl_P+s.G_k_P))

        b = np.array([0.,0.])

        b[0] = s.Area*(s.G_na_P*s.V_na_P+s.G_cl_P*s.V_cl_P+s.G_k_P*s.V_k_P)
        b[1] = (s.F*w_C*s.J_NKA_B + s.G_k_B*s.alpha_k_B*s.Area
                -s.Area*(s.G_na_P*s.V_na_P+s.G_cl_P*s.V_cl_P+s.G_k_P*s.V_k_P))

        x = np.linalg.solve(B,b)
        s.V_A, s.V_B = x

        # ionic currents
        s.I_CFTR = s.G_CFTR*s.alpha_CFTR*s.Area*(s.V_A-s.V_cl_A)

        # CFTR channel transports bicarbonate
        s.I_CFTR_B = 1/4*s.G_CFTR*s.alpha_CFTR*s.Area*(s.V_A-s.V_hco3_A)

        # sodium channels
        s.I_ENaC = s.G_ENaC*s.alpha_ENaC*s.Area*(s.V_A-s.V_na_A)

        # patassium channels

        # BK channel
        s.I_BK = s.P_BK * s.G_BK * s.Area * (s.V_A-s.V_k_A)

        # on the basolateral side, Is.k1 channel
        s.V_k_B = s.R*s.T/(z_k*s.F)*np.log(s.k_B/k_C)
        s.I_k_B = s.G_k_B * s.alpha_k_B*s.Area*(s.V_B - s.V_k_B)


        # Na-K-2Cl cotransporter NKCC1

        s.S_C = na_C*cl_C**2*k_C
        s.u_NKCC1 = (s.gamma_1 - s.gamma_2 * s.S_C)/(s.gamma_3+s.gamma_4*s.S_C)
        s.J_NKCC1 = s.alpha_NKCC1*s.u_NKCC1

        # acid transporters

        # Na-H exchanger NHE
        def chem_reaction_ex(p1, p2, p3, p4, k1, k2, k_1, k_2):
            return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1+k2*p2+k_1*p3+k_2*p4)

        s.J_NHE_A = s.alpha_NHE_A * chem_reaction_ex(na_A, h_C, na_C, h_A,
                                                 s.k1, s.k2, s.k_1, s.k_2)

        s.J_NHE_B = s.alpha_NHE_B * chem_reaction_ex(s.na_B, h_C, na_C, s.h_B,
                                                 s.k1, s.k2, s.k_1, s.k_2)

        # anion exchangers (AE2)

        s.J_AE2_A = s.alpha_AE2_A * chem_reaction_ex(cl_A, hco3_C, cl_C, hco3_A,
                                                 s.k3, s.k4, s.k_3, s.k_4)
        s.J_AE2_B = s.alpha_AE2_B * chem_reaction_ex(s.cl_B, hco3_C, cl_C, s.hco3_B,
                                                 s.k3, s.k4, s.k_3, s.k_4)

        # NA-HCO3 cotransporter (NBC)

        def chem_reaction_ex_2(p1, p2, p3, p4, k1, k2, k_1, k_2):
            return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1*p2+k2*k_1+k_2*p3*p4)

        s.J_NBC_A = s.alpha_NBC_A * chem_reaction_ex_2(na_C, hco3_C, na_A, hco3_A,
                                                   s.k5, s.k6, s.k_5, s.k_6)
        s.J_NBC_B = s.alpha_NBC_B * chem_reaction_ex_2(na_C, hco3_C, s.na_B, s.hco3_B,
                                                   s.k5, s.k6, s.k_5, s.k_6)

        # paracelluar currents
        V_t = s.V_A - s.V_B


        s.I_na_P = paracelluar_currents(s.G_na_P, na_A, na_C, na_C)
        s.I_cl_P = paracelluar_currents(s.G_cl_P, cl_A, cl_C, cl_C)
        s.I_k_P = paracelluar_currents(s.G_k_P, k_A, k_C, k_C)

        # carbon dioxide transport
        s.J_CDF_A = s.p_co2*(co2_C-co2_A)
        s.J_CDF_B = s.p_co2*(co2_C-s.co2_B)

        # bicarbonate buffer
        # h2o + co2 <-> hco3 + h

        s.J_buf_C = w_C*(s.k5*co2_C - s.k_5*hco3_C*h_C)
        s.J_buf_A = w_C*(s.k5*co2_A - s.k_5*hco3_A*h_A)

        rhs_na_c = -s.I_ENaC + s.F*w_C*(s.J_NKCC1 - 3*s.J_NKA_B
                                    + s.J_NHE_A + s.J_NHE_B
                                    + s.J_NBC_A + s.J_NBC_B)
        rhs_cl_c = s.I_CFTR + s.F*w_C*(s.J_AE2_A+s.J_AE2_B+2*s.J_NKCC1)
        rhs_k_c = -s.I_BK - s.I_k_B + s.F*w_C*(s.J_NKCC1 + 2*s.J_NKA_B)
        rhs_h_c = s.J_buf_C - s.J_NHE_A - s.J_NHE_B
        rhs_hco3_c = s.I_CFTR_B + s.F*w_C*(s.J_buf_C
                                       - s.J_AE2_A - s.J_AE2_B
                                       + s.J_NBC_A + s.J_NBC_B)
        rhs_co2_c = -s.J_buf_C - s.J_CDF_A - s.J_CDF_B

        rhs_na_a = s.I_ENaC + s.I_na_P - s.F*s.w_A*(s.J_NHE_A + s.J_NBC_A)
        rhs_cl_A = -s.I_CFTR - s.I_cl_P - s.F*s.w_A*s.J_AE2_A
        rhs_k_a = s.I_BK + s.I_k_B + s.I_k_P
        rhs_h_a = s.J_buf_A + s.J_NHE_A
        rhs_hco3_a = s.F*s.w_A*(s.J_buf_A + s.J_AE2_A - s.J_NBC_A) - s.I_CFTR_B
        rhs_co2_a = -s.J_buf_A + s.J_CDF_A + s.J_CDF_B


        return np.array( [s.rhs_w_c,
                    rhs_na_c, rhs_cl_c, rhs_k_c,
                    rhs_h_c, rhs_hco3_c, rhs_co2_c,
                    rhs_na_a, rhs_cl_A, rhs_k_a,
                    rhs_h_a, rhs_hco3_a, rhs_co2_a] )
