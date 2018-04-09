#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:43:17 2018

@author: plunder
"""
# import sympy as sp
import numpy as np
import numdifftools.nd_algopy as nda
import matplotlib.pyplot as plt

class Units:
    m = 1  # meter
    μm = 1e-6
    nm = 1e-9
    mM = 1  # millimolar (mM = mol/m³ = 1e-3 mol/L)
    μM = 1e-3  # micro molar
    nM = 1e-6  # nano molar
    J = 1  # joule
    K = 1   # kelvin
    s = 1   # seconds
    g = 1  # gram
    mol = 1
    C = 1  # °C... hard to translate... use K-C = constant
    Simens = 1

    w_dimless = μm**2  # currently not used
    w_dim = 1 # = μm**2 if w_dimless == 1

# helper functions
def chem_reaction_ex(p1, p2, p3, p4, k1, k2, k_1, k_2):
    return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1+k2*p2+k_1*p4+k_2*p3)


def chem_reaction_ex_2(p1, p2, p3, p4, k1, k2, k_1, k_2):
    return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1*p2+k2*+k_1+k_2*p3*p4)


class shelley_model:

    def __init__(self):

        s = self

        # Settings
        s.active_channels = {
                'CFTR' : False,
                'CFTR_B' : False,
                'ENaC' : False,
                'BK' : False,
                'IK1' : False,
                'NKA' : False,
                'NKCC1' : False,
                'NHE' : False,
                'AE2' : False,
                'NBC' : False,
                'P' : False,  # paracelluar
                'CDF' : False,
                }

        for k in s.active_channels.keys():
            s.active_channels[k] = True

        #s.active_channels['CFTR_B'] = True
        #s.active_channels['P'] = True
        #s.active_channels['ENaC'] = True
        #s.active_channels['BK'] = True
        #s.active_channels['IK1'] = True
        #s.active_channels['CDF'] = True
        s.active_channels['NHE'] = True

        s.active_buffers = {
                'CO2' : True
                }

        # Intersitial compartment
        s.na_B = 140.7 * Units.mM
        s.cl_B = 125. * Units.mM
        s.k_B = 5 * Units.mM
        s.pH = 7.2
        s.hco3_B = 42.9 * Units.mM
        s.co2_B = 1.28 * Units.mM
        s.h_B = 10**(-s.pH) * 1e3 * Units.mM

        # TODO: adjust diameter according to speed

        s.Area = 20.9 * Units.nm**2
        s.V_W = 18.05e-6 * Units.m**3/Units.mol
        s.chi_C = 0.00 * Units.mol # TODO: check value # fix value
        s.chi_A = 0.00 * Units.mol # TODO: check value
        # fix value, needed for electro neutrality
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
        s.kbuf = 0.03/Units.s

        s.k_1 = 1.4294e11/Units.s
        s.k_2 = 1.7857e2/Units.s
        s.k_3 = 1.0652e8/Units.s
        s.k_4 = 5.1418/Units.s
        s.k_5 = 1e8/(Units.mM*Units.s)
        s.k_6 = -1.9352e-1/(Units.mM*Units.s)
        s.k_buf = 20./(Units.mM*Units.s)

        # 'unit of conductance' aka uCond
        uCond = Units.Simens/Units.m**2

        s.G_CFTR = 90 * uCond
        s.alpha_CFTR = 1 * Units.m**(-2)

        s.G_CFTR_B = 90 * uCond
        s.alpha_CFTR_B = 1 * Units.m**(-2)

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
        s.alpha_BK = 1 * Units.m**(-2)

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
        s.G_IK1 = 1 * uCond
        s.alpha_IK1 = 1. * Units.m**(-2)

        s.names = ['w_C', '[Na^{+}]_A', '[Cl^{-}]_A', '[K^{+}]_A', '[H^{+}]_A',
                   '{[HCO^{-}_3]}_A', '[{CO_2}]_A',
                   '[Na^{+}]_C', '[Cl^{-}]_C', '[K^{+}]_C', '[H^{+}]_C',
                   '{[HCO^{-}_3]}_C', '{[CO_2]}_C']


    def fix_initial_conditions(self, x_0):
        # electroneutrality
        (w_C,
         na_A, cl_A, k_A, h_A, hco3_A, co2_A,
         na_C, cl_C, k_C, h_C, hco3_C, co2_C) = x_0

        test = w_C*Units.w_dim*(na_C - cl_C + h_C + k_C - hco3_C) - self.chi_C
        if test < 0:
            na_C = self.chi_C/(w_C*Units.w_dim) - (- cl_C + h_C + k_C - hco3_C)
        else:
            cl_C = -self.chi_C/(w_C*Units.w_dim) + (na_C + h_C + k_C - hco3_C)


        test = self.w_A*Units.w_dim*(na_A - cl_A + h_A + k_A - hco3_A) - self.chi_A
        if test < 0:
            na_A = self.chi_A/(self.w_A*Units.w_dim) - (- cl_A + h_A + k_A - hco3_A)
        else:
            cl_A = -self.chi_A/(self.w_A*Units.w_dim) + (na_A + h_A + k_A - hco3_A)


        # co2 equilibrium
        co2_A = self.k_buf/self.kbuf*( hco3_A * h_A )
        co2_C = self.k_buf/self.kbuf*( hco3_C * h_C )

        return np.array([w_C,
                         na_A, cl_A, k_A, h_A, hco3_A, co2_A,
                         na_C, cl_C, k_C, h_C, hco3_C, co2_C])


    def activate_channel(self,name,state=True):
        self.active_channels[name] = state

    def activate_buffer(self,name,state=True):
        self.active_buffers[name] = state

    def ode(self, x):
        """returns np.ndarray (with shape=(13,))

            parameters arrangement:
                x = [w_C,
                na_A, cl_A, k_A, h_A, hco3_A, co2_A,
                na_C, cl_C, k_C, h_C, hco3_C, co2_C]

                where w_C denotes the volume of the duct cell,
                *_A denotes the concentration on the apical side
                and *_C denotes the concentration in the duct cells



        .. note::
        .. warning:: all concentrations must be positive
          (since log(*_{A,C]}) is may computed in the model.)
        """
        # are we in an automatic differentiation run?
        auto_diff = not type(x[0]) is np.float64


        # to shorten the notations, we set
        s = self

        # unkowns
        (w_C,
         na_A, cl_A, k_A, h_A, hco3_A, co2_A,
         na_C, cl_C, k_C, h_C, hco3_C, co2_C) = x


        (s.w_C,
         s.na_A, s.cl_A, s.k_A, s.h_A, s.hco3_A, s.co2_A,
         s.na_C, s.cl_C, s.k_C, s.h_C, s.hco3_C, s.co2_C) = x

        con_A = [na_A, cl_A, k_A, h_A, hco3_A, co2_A]
        con_B = [s.na_B, s.cl_B, s.k_B, s.h_B, s.hco3_B]
        con_C = [na_C, cl_C, k_C, h_C, hco3_C, co2_C]

        z_cl = -1
        z_hco3 = -1
        z_na = +1
        z_k = +1

        # return values,
        s.rhs_w_C = 0

        s.rhs_na_C = 0
        s.rhs_cl_C = 0
        s.rhs_k_C = 0
        s.rhs_h_C = 0
        s.rhs_hco3_C = 0
        s.rhs_co2_C = 0

        s.rhs_na_A = 0
        s.rhs_cl_A = 0
        s.rhs_k_A = 0
        s.rhs_h_A = 0
        s.rhs_hco3_A = 0
        s.rhs_co2_A = 0

        # compute the matrix s.t. B*[s.V_A,s.V_B] == b equals the constraint

        B = nda.algopy.zeros(shape=(2,2),dtype=x[0])
        b = nda.algopy.zeros(shape=(2,1),dtype=x[0])

        s.J_w_A = Units.w_dim**(-1)*s.L_A*s.V_W*(sum(con_A)-sum(con_C) + (-s.chi_C/(w_C*Units.w_dim) + s.Psi))
        s.J_w_B = Units.w_dim**(-1)*s.L_B*s.V_W*(sum(con_C)-sum(con_B) + (s.chi_C/(w_C*Units.w_dim)))

        # to resolve the algebraic constraint, we need to
        # compute all potentials first, then find the
        # current volumes and after that compute the currents

        if s.active_channels['CFTR']:
            # nerst equation
            s.V_cl_A = s.R*s.T/(z_cl*s.F) * np.log(cl_A/cl_C)

            factor = s.G_CFTR*s.alpha_CFTR*s.Area
            b[0] = b[0] + factor*s.V_cl_A
            B[0, 0] += factor

            # init the current with the factor, s.t. later only the multiplication
            # with V_A - V_k_A is needed
            s.I_CFTR = factor

        if s.active_channels['CFTR_B']:
            # CFTR channel transports bicarbonate
            s.V_hco3_A = s.R*s.T/(z_hco3*s.F)*np.log(hco3_A / hco3_C)

            factor = 0.25*s.G_CFTR_B*s.alpha_CFTR_B*s.Area
            b[0] += factor*s.V_hco3_A
            B[0, 0] += factor

            s.I_CFTR_B = factor

        if s.active_channels['ENaC']:
            # sodium channels
            s.V_na_A = s.R*s.T/(z_na*s.F)*np.log(na_A/na_C)

            factor = s.G_ENaC*s.alpha_ENaC*s.Area
            b[0] += factor*s.V_na_A
            B[0, 0] += factor

            s.I_ENaC = factor

        if s.active_channels['BK']:
            # potassium channels
            # open prob via Hill Function, calcium activated
            ca_C_eta = s.ca_C**s.eta
            s.P_BK = ca_C_eta / (ca_C_eta + s.K_m_BK**s.eta)

            s.V_k_A = s.R*s.T/(z_k*s.F) * np.log(k_A/k_C)

            factor = s.P_BK*s.G_BK*s.alpha_BK*s.Area
            b[0] += factor*s.V_k_A
            B[0, 0] += factor

            s.I_BK = factor

        if s.active_channels['IK1']:
            # on the basolateral side, Is.k1 channel
            s.V_k_B = s.R*s.T/(z_k*s.F)*np.log(s.k_B/k_C)

            factor = s.G_IK1*s.alpha_IK1*s.Area
            b[1] += factor*s.V_k_B
            B[1, 1] += factor

            s.I_IK1 = factor

        if s.active_channels['P']:
            # paracelluar potentials
            s.V_na_P = s.R*s.T/(z_na*s.F)*np.log(na_A/na_C)
            s.V_cl_P = s.R*s.T/(z_cl*s.F)*np.log(cl_A/cl_C)
            s.V_k_P = s.R*s.T/(z_k*s.F)*np.log(k_A/k_C)

            factor = s.G_na_P*s.Area
            b[0] += factor*s.V_na_P
            b[1] -= b[0]
            B[0, 0] += factor
            B[0, 1] -= factor
            B[1, 0] += factor
            B[1, 1] -= factor
            s.I_na_P = factor

            factor = s.G_cl_P*s.Area
            b[0] += factor*s.V_cl_P
            b[1] -= b[0]
            B[0, 0] += factor
            B[0, 1] -= factor
            B[1, 0] += factor
            B[1, 1] -= factor
            s.I_cl_P = factor

            factor = s.G_k_P*s.Area
            b[0] += factor*s.V_k_P
            b[1] -= b[0]
            B[0, 0] += factor
            B[0, 1] -= factor
            B[1, 0] += factor
            B[1, 1] -= factor
            s.I_k_P = factor

        if s.active_channels['NKA']:
            # Na-K_ATPase (NKA)
            s.u_NKA_B = s.r_NKA*(s.k_B**2*na_C**3)/(s.k_B**2 + s.beta_NKA*na_C**3)
            s.J_NKA_B = s.alpha_NKA*s.u_NKA_B

            b[0] += s.F*(w_C*Units.w_dim)*s.J_NKA_B

        # membrane potential
        # algebraic constraint!
        # solve for V_A, V_B to get an ODE instead of a DAE
        #if nda.algopy.det(B) ==
        #V_AB = nda.algopy.solve(B,b)

        #s.V_A, s.V_B = V_AB
        s.V_A = 0
        s.V_B = 0
        # now we compute the currents and build up the RHS
        # of the equations, mainly by summing up the currents

        pre_factor = 1./(s.F*w_C*Units.w_dim)

        if s.active_channels['CFTR']:
            # ionic currents
            s.I_CFTR *= s.V_A-s.V_cl_A

            J = pre_factor * s.I_CFTR
            s.rhs_cl_C += J
            s.rhs_cl_A -= J

        if s.active_channels['CFTR_B']:
            # CFTR channel transports bicarbonate
            s.I_CFTR_B *= s.V_A-s.V_hco3_A

            J = pre_factor * s.I_CFTR_B
            s.rhs_hco3_C += J
            s.rhs_hco3_A -= J

        if s.active_channels['ENaC']:
            # sodium channels
            s.I_ENaC *= s.V_A-s.V_na_A

            J = pre_factor * s.I_ENaC
            s.rhs_na_C -= J
            s.rhs_na_A += J


        if s.active_channels['BK']:
            # patassium channels
            # BK channel
            s.I_BK *= s.V_A-s.V_k_A

            J = pre_factor * s.I_BK
            s.rhs_k_C -= J
            s.rhs_k_A += J

        if s.active_channels['IK1']:
            # on the basolateral side, Is.k1 channel
            s.I_IK1 *= s.V_B - s.V_k_A

            J = pre_factor * s.I_IK1
            s.rhs_k_C -= J
            s.rhs_k_A += J

        if s.active_channels['NKA']:
            s.rhs_na_C -= 3*s.J_NKA_B
            s.rhs_k_C += 2*s.J_NKA_B

        if s.active_channels['NKCC1']:
            # Na-K-2Cl cotransporter NKCC1
            s.S_C = na_C*cl_C**2*k_C
            s.u_NKCC1 = (s.gamma_1 - s.gamma_2 * s.S_C)/(s.gamma_3+s.gamma_4*s.S_C)
            s.J_NKCC1 = s.alpha_NKCC1*s.u_NKCC1

            s.rhs_na_C += s.J_NKCC1
            s.rhs_na_A -= s.J_NKCC1

            s.rhs_cl_C += 2*s.J_NKCC1
            s.rhs_cl_A -= 2*s.J_NKCC1

            s.rhs_k_C += s.J_NKCC1
            s.rhs_k_A -= s.J_NKCC1

        # acid transporters

        if s.active_channels['NHE']:
            # Na-H exchanger NHE
            s.J_NHE_A = s.alpha_NHE_A * chem_reaction_ex(na_A, h_C, na_C, h_A,
                                                     s.k1, s.k2, s.k_1, s.k_2)

            s.J_NHE_B = s.alpha_NHE_B * chem_reaction_ex(s.na_B, h_C, na_C, s.h_B,
                                                     s.k1, s.k2, s.k_1, s.k_2)

            s.rhs_na_C += s.J_NHE_A + s.J_NHE_B
            s.rhs_na_A -= s.J_NHE_A

            s.rhs_h_C -= s.J_NHE_A + s.J_NHE_B
            s.rhs_h_A += s.J_NHE_A


        if s.active_channels['AE2']:
            # anion exchangers (AE2)

            s.J_AE2_A = s.alpha_AE2_A * chem_reaction_ex(cl_A, hco3_C, cl_C, hco3_A,
                                                     s.k3, s.k4, s.k_3, s.k_4)
            s.J_AE2_B = s.alpha_AE2_B * chem_reaction_ex(s.cl_B, hco3_C, cl_C, s.hco3_B,
                                                     s.k3, s.k4, s.k_3, s.k_4)

            s.rhs_k_C += s.J_NKCC1
            s.rhs_k_A -= s.J_NKCC1

        if s.active_channels['NBC']:
            # NA-HCO3 cotransporter (NBC)

            s.J_NBC_A = s.alpha_NBC_A * chem_reaction_ex_2(na_C, hco3_C, na_A, hco3_A,
                                                       s.k5, s.k6, s.k_5, s.k_6)
            s.J_NBC_B = s.alpha_NBC_B * chem_reaction_ex_2(na_C, hco3_C, s.na_B, s.hco3_B,
                                                       s.k5, s.k6, s.k_5, s.k_6)

            s.rhs_na_C += self.J_NBC_A + self.J_NBC_B
            s.rhs_na_A -= self.J_NBC_A

            s.rhs_hco3_C += self.J_NBC_A + self.J_NBC_B
            s.rhs_hco3_A -= self.J_NBC_A

        # paracelluar currents
        V_t = s.V_A - s.V_B

        if s.active_channels['P']:
            s.I_na_P *= V_t - s.V_na_P
            s.I_cl_P *= V_t - s.V_cl_P
            s.I_k_P *= V_t - s.V_k_P

            s.rhs_na_A += pre_factor*s.I_na_P
            s.rhs_cl_A -= pre_factor*s.I_cl_P
            s.rhs_k_A += pre_factor*s.I_k_P

        if s.active_channels['CDF']:
            # carbon dioxide transport
            s.J_CDF_A = s.p_co2*(co2_C-co2_A)
            s.J_CDF_B = s.p_co2*(co2_C-s.co2_B)

            s.rhs_co2_C -= s.J_CDF_A + s.J_CDF_B
            s.rhs_co2_A += s.J_CDF_A

        if s.active_buffers['CO2']:
            # bicarbonate buffer
            # h2o + co2 <-> hco3 + h

            s.J_buf_C = w_C*Units.w_dim*(s.kbuf*co2_C - s.k_buf*hco3_C*h_C)
            s.J_buf_A = w_C*Units.w_dim*(s.kbuf*co2_A - s.k_buf*hco3_A*h_A)

            s.rhs_hco3_C += s.J_buf_C
            s.rhs_hco3_A += s.J_buf_A

            s.rhs_co2_C -= s.J_buf_C
            s.rhs_co2_A -= s.J_buf_A

            s.rhs_h_C += s.J_buf_C
            s.rhs_h_A += s.J_buf_A

        # test if the wanted equalities hold!
        # s.test_A = s.I_ENaC + s.I_na_P + s.I_CFTR + s.I_cl_P + s.I_BK + s.I_k_P + s.I_CFTR_B
        # s.test_B = -s.F*w_C*s.J_NKA_B - (s.I_na_P + s.I_cl_P + s.I_k_P) + s.I_IK1


        s.rhs_w_c = s.Area*(s.J_w_B - s.J_w_A)

        # product rule (to convert the change in mass into change in concentrations)
        s.rhs_na_C += s.rhs_na_C + s.rhs_w_C / w_C * na_C
        s.rhs_cl_C += s.rhs_cl_C + s.rhs_w_C / w_C * cl_C
        s.rhs_k_C += s.rhs_k_C + s.rhs_w_C / w_C * k_C
        s.rhs_hco3_C += s.rhs_hco3_C + s.rhs_w_C / w_C * hco3_C

        s.rhs = nda.algopy.zeros(shape=(13,),dtype=x[0])
        s.rhs[0] = s.rhs_w_C
        s.rhs[1] = s.rhs_na_C
        s.rhs[2] = s.rhs_cl_C
        s.rhs[3] = s.rhs_k_C
        s.rhs[4] = s.rhs_h_C
        s.rhs[5] = s.rhs_hco3_C
        s.rhs[6] = s.rhs_co2_C
        s.rhs[7] = s.rhs_na_A
        s.rhs[8] = s.rhs_cl_A
        s.rhs[9] = s.rhs_k_A
        s.rhs[10] = s.rhs_h_A
        s.rhs[11] = s.rhs_hco3_A
        s.rhs[12] = s.rhs_co2_A



        if not auto_diff:
            if np.any(np.isnan(s.rhs)):
                return np.zeros(shape=(13,))

        return s.rhs


    def plot_state(s):
        x = [ chem + region for chem in range(0,6) for region in [0,0.2,0.4] ]

        names = ['na','cl','k','h','hco3','co2']
        y = [ getattr(s,name + '_' + region ) for name in names for region in ['A','C','B']]

        plt.bar(x, y, width=0.15)
        plt.xticks( np.arange(0,6) , names )

        plt.show()



    def plot_change(s):
        x = [ chem + region for chem in range(0,6) for region in [0,0.2] ]

        names = ['na','cl','k','h','hco3','co2']
        y = [ getattr(s,'rhs_' + name + '_' + region ) for name in names for region in ['A','C']]

        plt.bar(x, y, width=0.15)
        plt.xticks( np.arange(0,6) , names )

        plt.show()

