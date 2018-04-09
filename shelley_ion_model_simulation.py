#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 17:39:27 2018

@author: plunder
"""

import numpy as np
import numdifftools.nd_algopy as nda
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, Radau, BDF, LSODA
from shelley_ion_model_system import shelley_model, Units

# init matplotlib
params = {'legend.fontsize': 'x-large','figure.figsize': (8, 6),
         'axes.labelsize': 'x-large','axes.titlesize':'x-large',
         'xtick.labelsize':'x-large','ytick.labelsize':'x-large',
         'text.usetex':True}
plt.rcParams.update(params)


# here the modeling starts:

s = shelley_model()

mM = Units.mM

# initial conditions
w_C = 2.090 * Units.w_dimless
na_A = 141.6 * mM
na_C = 141.6 * mM
cl_A = 112.6 * mM
cl_C = 112.6 * mM
k_A = 35.6 * mM
k_C = 35.6 * mM
h_A = 10**(-7.2) * 1e3 * mM
h_C = 10**(-7.2) * 1e3 * mM
hco3_A = 42.9 * mM
hco3_C = 42.9 * mM
co2_A = 1.288 * mM
co2_C = 1.228 * mM


x_0 = [w_C,
        na_A, cl_A, k_A, h_A, hco3_A, co2_A,
        na_C, cl_C, k_C, h_C, hco3_C, co2_C]

x_0 = s.fix_initial_conditions(x_0)
x_cur = x_0

# note forward mode ist faster than reverse mode in our case!
ode_jac = nda.Jacobian(s.ode, method='forward')


t0 = 0
t_bound = 0.1

ode = Radau( lambda t,x: s.ode(x), y0 = x_0, t0 = t0, t_bound = t_bound, jac = lambda t,x: ode_jac(x) )

T = np.array([t0])
Y = np.array(x_0)
Y = Y[:,np.newaxis]


#ode.max_step = 1e-7


while ode.t < t_bound:

    try:
        ode.step()
    except ValueError:
        print("Value Error in ode.stop()")
        break
    except RuntimeError as ex:
        print("Runtime Error in ode.step()")
        break
    T = np.append(T,ode.t)
    Y = np.append(Y,ode.y[:,np.newaxis],axis=1)
    x_cur = ode.y


#y = solve_ivp(lambda t, x: s.ode(x), t_span, x_0,method='BDF',jac=lambda t,x: ode_jac(x))
#print("Result: %s" % y.message )


# plot
plt.plot(T,Y[0,:],label='$'+s.names[0]+'$')
plt.legend()
plt.show()

cols = [[1,0,0],[1,0.5,0],[0,1,0],[0,0,1],[1,0,1],[0,0,0]]

for i in range(1,7):
    plt.plot(T,Y[i][:],color=cols[i-1],label='$'+s.names[i]+'$')


plt.legend(loc=2, ncol=3, borderaxespad=0.)
plt.show()

for i in range(7,13):
    plt.plot(T,Y[i][:],color=cols[i-7],label='$'+s.names[i]+'$')


plt.legend(loc=2, ncol=3, borderaxespad=0.)
plt.show()
