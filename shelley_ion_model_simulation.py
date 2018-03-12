#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:39:27 2018

@author: plunder
"""

import numpy as np

from shelley_ion_model_system import shelley_model, Units

s = shelley_model()



from scipy.integrate import solve_ivp

t_span = (0,1)

w_C = 0.001 * Units.nm**2
na_A, na_C = (0., 0.)
cl_A, cl_C = (0., 0.)
k_A, k_C = (0., 0.)
h_A, h_C = (0., 0.)
hco3_A, hco3_C = (0., 0.)
co2_A, co2_C = (0., 0.)

# [w_C,
#  na_A, cl_A, k_A, h_A, hco3_A, co2_A,
#  na_C, cl_C, k_C, h_C, hco3_C, co2_C]

x_0 = [w_C,
        na_A, cl_A, k_A, h_A, hco3_A, co2_A,
        na_C, cl_C, k_C, h_C, hco3_C, co2_C]

x_0 = np.abs(0.1+np.random.randn(13))

s.ode(x_0)

solve_ivp(lambda t,x: s.ode(x), t_span, x_0)
