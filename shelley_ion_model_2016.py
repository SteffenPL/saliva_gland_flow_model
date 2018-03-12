#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 16:38:58 2018

@author: plunder
"""

from sympy import *
#from sympy.physics import units as u

init_printing(use_latex=True)
F = symbols('F')
w_A = symbols('w_A') # water flux into the lumen
w_C = symbols('w_C') # water flux into the duct cell
na_A, cl_A, k_A, h_A, hco3_A, co2_A = symbols('[Na^{+}]_A,\
                            [Cl^{-}]_A,\
                            [K^{+}]_A,\
                            [H^{+}]_A,\
                            [HCO_3^{-}]_A\
                            [CO_2]_A')

na_B, cl_B, k_B, h_B, hco3_B, co2_B = symbols('[Na^{+}]_B,\
                            [Cl^{-}]_B,\
                            [K^{+}]_B,\
                            [H^{+}]_B,\
                            [HCO_3^{-}]_B\
                            [CO_2]_B')

na_C, cl_C, k_C, h_C, hco3_C, co2_C = symbols('[Na^{+}]_C,\
                            [Cl^{-}]_C,\
                            [K^{+}]_C,\
                            [H^{+}]_C,\
                            [HCO_3^{-}]_C\
                            [CO_2]_C')



con_A = [na_A,cl_A,k_A,h_A,hco3_A]
con_B = [na_B,cl_B,k_B,h_B,hco3_B]
con_C = [na_C,cl_C,k_C,h_C,hco3_C]


# Na+ term
I_ENaC = symbols('I_{ENaC}')
J_NKCC1, J_NKA_B = symbols('J_{NKCC1}, J_{{NKA}_B}')
J_NHE_A, J_NHE_B, J_NBC_A, J_NBC_B = symbols('J_{{NHE}_A}, J_{{NHE}_B}, J_{{NBC}_A}, J_{{NBC}_B}')

#Cl- term
I_CFTR = symbols('I_{CFTR}')
J_AE2_A, J_AE2_B = symbols('J_{{AE2}_A}, J_{{AE2}_B}')

# K+ term
I_BK, I_K_B = symbols('I_{BK}, I_{K_B}')

#H+ term
J_buf_A, J_buf_C = symbols('J_{\\text{buf}_A}, J_{\\text{buf}_C}')

#HCO3 term
I_CFTR_B = symbols('I_{CFTR,B}')

#CO2 term
J_CDF_A, J_CDF_B, J_CDF_P = symbols('J_{CDF_A},J_{CDF_B}, J_{CDF_P}')

I_na_P, I_cl_P, I_k_P = symbols('I^{Na}_P, I^{Cl}_P, I^K_P')



#odes

# water flux
A = symbols('A')
J_w_A, J_w_B = symbols('J^w_A, J^w_B')
L_A, L_B, V_W = symbols('L_A, L_B, V_W')
chi_C, Psi = symbols('\\chi_C, \\Psi^{-}')

J_w_A = L_A*V_W*( sum(con_A)-sum(con_C) + ( -chi_C/w_C + Psi ))
J_w_B = L_B*V_W*( sum(con_C)-sum(con_B) + ( chi_C/w_C ))
rhs_w_c = A*(J_w_B - J_w_A)

# ionic currents
G_CFTR, alpha_CFTR = symbols('G_{CFTR}, \\alpha_{CFTR}')
V_A, V_cl_A, V_hco3_A, V_na_A, V_k_A = symbols('V_A, V^{Cl}_A, V^{HCO_3}_A, V^{Na}_A, V^K_A')
V_B, V_cl_B, V_hco3_B, V_na_B, V_k_B = symbols('V_B, V^{Cl}_B, V^{HCO_3}_B, V^{Na}_B, V^K_B')
I_CFTR = G_CFTR*alpha_CFTR*A*(V_A-V_cl_A)

R, T, F, z = symbols('R, T, F, z')

# nerst equation
V_cl_A = R*T/(z*F) * log(cl_A/cl_C)

# CFTR channel transports bicarbonate
I_CFTR_B = 1/4*G_CFTR*alpha_CFTR*A*(V_A-V_hco3_A)
V_hco3_A = R*T/(z*F)*log( hco3_A / hco3_C )

# sodium channels
G_ENaC, alpha_ENaC = symbols('G_{ENaC}, \\alpha_{ENaC}')
I_ENaC = G_ENaC*alpha_ENaC*A*(V_A-V_na_A)
V_na_A = R*T/(z*F)*log( na_A/na_C )

# patassium channels

# open prob via Hill Function, calcium activated
P_BK, G_BK, K_m_BK, eta = symbols('P_BK, G_BK, K^m_{BK}, \\eta')
ca_C = symbols('[Ca^{2+}]_C')
P_BK = 1/( 1 + (K_m_BK/ca_C)**eta )

I_BK = P_BK * G_BK * A* (V_A-V_k_A)
V_k_A = R*T/(z*F) * log( k_A/k_C )

# on the basolateral side, IK1 channel
G_k_B, alpha_k_B = symbols('G_{K_B}, \\alpha_{K_B}')
I_k_B = G_k_B * alpha_k_B*A*(V_B - V_k_B)
V_k_B = R*T/(z*F)*log(k_B/k_C)

# Na-K_ATPase (NKA)
ny_NKA, r_NKA, alpha_NKA, beta_NKA = symbols('ny_{NKA}, r_{NKA}, \\alpha_{NKA}, \\beta_{NKA}')
u_NKA = symbols('u_{NKA}')

ny_NKA = r_NKA*(k_B**2*na_C**3)/(k_B**2 + beta_NKA*na_C**3)
J_NKA = alpha_NKA*u_NKA

# Na-K-2Cl cotransporter NKCC1

gamma_1, gamma_2, gamma_3, gamma_4 = symbols('gamma_1,gamma_2,gamma_3,gamma_4')
alpha_NKCC1, S_C = symbols('\\alpha_{NKCC1}, S_C')

S_C = na_C*cl_C**2*k_C
u_NKCC1 = (gamma_1 - gamma_2 * S_C)/(gamma_3+gamma_4*S_C)
J_NKCC1 = alpha_NKCC1*u_NKCC1

# acid transporters
#
# Find out later!
#

# Na-H exchanger NHE
k1,k2,k_1,k_2 = symbols('k^+_1, k^+_1, k^+_1, k^+_1')
h_e, na_e = symbols('H_e, Na_e')
def chem_reaction_ex(p1,p2,p3,p4,k1,k2,k_1,k_2):
    return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1+k2*p2+k_1*p3+k_2*p4)

alpha_NHE = symbols('\\alpha_{NHE}')

J_NHE = alpha_NHE * chem_reaction_ex(na_e,h_C,na_C,h_e,k1,k2,k_1,k_2)

# anion exchangers (AE2)

k3,k4,k_3,k_4 = symbols('k^+_3, k^+_4, k^+_3, k^+_4')
alpha_AE2 = symbols('\\alpha_{AE2}')
cl_e, hco3_e = symbols('[Cl]_e, [HCO^{-}_{3}]_e')
J_AE2 = alpha_AE2 * chem_reaction_ex(cl_e,hco3_C,cl_C,hco3_e,k3,k4,k_3,k_4)

# NA-HCO3 cotransporter (NBC)

def chem_reaction_ex_2(p1,p2,p3,p4,k1,k2,k_1,k_2):
    return (k1*k2*p1*p2 - k_1*k_2*p3*p4)/(k1*p1*p2+k2*k_1+k_2*p3*p4)

k5,k6,k_5,k_6 = symbols('k^{+}_5, k^{+}_6, k^{+}_5, k^+{}}_6')
alpha_NBC = symbols('\\alpha_{NBC}')
J_NBC = alpha_NBC * chem_reaction_ex_2(na_e,hco3_e,na_C,hco3_C,k5,k6,k_5,k_6)

# paracelluar currents
V_t = V_A - V_B
def paracelluar_currents(G_n_P,V_n_P,n_A,n_C):
    I_n_P = G_n_P * A * (V_t - V_n_P)
    V_n_P = R*T/(z*F)*log( n_A/n_C)
    return (I_n_P,V_n_P)

G_na_P, G_cl_P, G_k_P = symbols('G^{Na}_P, G^{Cl}_P, G^{K}_P')
V_na_P, V_cl_P, V_k_P = symbols('V^{Na}_P, V^{Cl}_P, V^{K}_P')
I_na_P, V_na_P = paracelluar_currents(G_na_P,V_na_P,na_A,na_C)
I_cl_P, V_cl_P = paracelluar_currents(G_cl_P,V_cl_P,cl_A,cl_C)
I_k_P, V_k_P = paracelluar_currents(G_k_P,V_k_P,k_A,k_C)

# carbon dioxide transport
p_co2 = symbols('p_{C0_2}')
J_CDF_A = p_co2*(co2_C-co2_A)
J_CDF_B = p_co2*(co2_C-co2_B)

# bicarbonate buffer
# h2o + co2 <-> hco3 + h

kbuf,k_buf = symbols('k^{+},k^{-}')
J_buf_C = w_C*(k5*co2_C - k_5*hco3_C*h_C)
J_buf_A = w_C*(k5*co2_A - k_5*hco3_A*h_A)

# membrane potential
# algebraic constraint!
rhs_kirchhoff_A = I_ENaC+I_na_P+I_CFTR+I_cl_P+I_BK+I_k_P+I_CFTR_B
rhs_kirchhoff_B = -F*w_C*J_NKA_B-I_na_P-I_cl_P+I_k_B-I_k_P

rhs_na_c = -I_ENaC + F*w_C*(J_NKCC1 - 3*J_NKA_B + J_NHE_A + J_NHE_B + J_NBC_A + J_NBC_B)
rhs_cl_c = I_CFTR + F*w_C*(J_AE2_A+J_AE2_B+2*J_NKCC1)
rhs_k_c = -I_BK - I_K_B + F*w_C*( J_NKCC1 + 2*J_NKA_B)
rhs_h_c = J_buf_C - J_NHE_A - J_NHE_B
rhs_hco3_c = I_CFTR_B + F*w_C*(J_buf_C - J_AE2_A - J_AE2_B + J_NBC_A + J_NBC_B)
rhs_co2_c = -J_buf_C - J_CDF_A - J_CDF_B

rhs_na_a = I_ENaC + I_na_P - F*w_A*(J_NHE_A + J_NBC_A)
rhs_cl_a = -I_CFTR - I_cl_P - F*w_A*J_AE2_A
rhs_k_a = I_BK + I_K_B + I_k_P
rhs_h_a = J_buf_A + J_NHE_A
rhs_hco3_a = F*w_A*(J_buf_A + J_AE2_A - J_NBC_A) - I_CFTR_B
rhs_co2_a = -J_buf_A + J_CDF_A + J_CDF_P

rhs_expr = [rhs_w_c,rhs_kirchhoff_A,rhs_kirchhoff_B,rhs_na_c,rhs_cl_c,
            rhs_k_c,rhs_h_c,rhs_hco3_c,rhs_co2_c,rhs_na_a,rhs_cl_a,
            rhs_k_a,rhs_h_a,rhs_hco3_a,rhs_co2_a]


for rhs in rhs_expr:
    display(rhs)