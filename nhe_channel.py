#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 18:43:52 2018

@author: plunder
"""

import numpy as np
import numdifftools.nd_algopy as nda
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, Radau, BDF, LSODA
from shelley_ion_model_system import shelley_model, Units, chem_reaction_ex

plt.interactive(True)


s = shelley_model()

na_B = 100.
h_B = 10**(-4.35)

for dim in [2,3]:

    if dim == 2:
        def f(t,x):
            (na_A,h_A,na_C,h_C) = x
            #J = s.alpha_NHE_A*(s.k1*s.k2*na_A * h_C - s.k_1*s.k_2*na_C * h_A) / (s.k1*na_A + s.k2*h_C + s.k_1*h_A + s.k_2*na_C)

            J = s.alpha_NHE_A * chem_reaction_ex(na_A, h_C, na_C, h_A,
                                                     s.k1,s.k2,s.k_1,s.k_2)

            J_2 = s.alpha_NHE_B * chem_reaction_ex(na_B, h_C, na_C, h_B,
                                                     s.k1,s.k2,s.k_1,s.k_2)
            return [-J,J,J+J_2,-J-J_2]

    else:
        def f(t,x):
            (na_A,h_A,na_C,h_C,na_B,h_B) = x

            J = s.alpha_NHE_A * chem_reaction_ex(na_A, h_C, na_C, h_A,
                                                     s.k1,s.k2,s.k_1,s.k_2)

            J_2 = s.alpha_NHE_B * chem_reaction_ex(na_B, h_C, na_C, h_B,
                                                     s.k1,s.k2,s.k_1,s.k_2)

            return [-J,J,J+J_2,-J-J_2,-J_2,J_2]



    if dim == 2:
        x_0 = [1,0.1+10**(-4.2),1,10**(-4.2)]
    else:
        x_0 = [1,0.1+10**(-4.2),1,10**(-4.2),1,10**(-4.2)]
    t0 = 0
    t_bound = 1.e-2
    ode = BDF( f, y0 = x_0, t0 = t0, t_bound = t_bound )
    ode.max_step = t_bound/10000
    T = np.array([t0])
    Y = np.array(x_0)
    Y = Y[:,np.newaxis]

    F = np.array(f(0,x_0))
    F = F[:,np.newaxis]

    if dim == 2:
        ode.max_step = 1.e-7


    p_last = 0
    p_plot = 0.1

    while ode.t < t_bound:

        p = ode.t/t_bound
        if p - p_last > p_plot:
            print("Progress %.2f\n" % p)
            p_last = p

        try:
            ode.step()
        except ValueError:
            break
        T = np.append(T,ode.t)
        Y = np.append(Y,ode.y[:,np.newaxis],axis=1)
        # x_cur = ode.y[:,np.newaxis]
        F = np.append(F, np.array(f(0,ode.y))[:,np.newaxis],axis=1)

    plt.plot(T,Y.T)
    params = {'text.usetex':False}
    plt.rcParams.update(params)
    plt.legend(['na_A','h_A','na_C','h_C'])
    plt.show()


    plt.plot(T[10:],F[:,10:].T)
    params = {'text.usetex':False}
    plt.rcParams.update(params)
    plt.legend(['na_A','h_A','na_C','h_C'])
    plt.show()