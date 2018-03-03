#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: plunder

@brief: we model a simple potential flow, where a influx appears
at all points of the interior and and outflux is possible at the plane x=1.
For the concentrations we
"""

from fenics import *

# Parameters
waterInFlux = 0.4

# reaction coefficients


# Create classes for defining parts of the boundaries and the interior
# of the domain

class Hole(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],1)


class ActiveArea(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0],(0.,0.1)) and between(x[1],(0,0.2))

# Initialize sub-domain instances
hole = Hole()
activeArea = ActiveArea()

# Define mesh
#mesh = UnitSquareMesh(64, 64)
mesh = Mesh('nontrivial_mesh_large.xml')


# Initialize mesh function for interior domains
domains = MeshFunction("size_t",mesh,2)
domains.set_all(0)
activeArea.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = MeshFunction("size_t", mesh,1)
boundaries.set_all(0)
hole.mark(boundaries, 1)

# Define new measures associated with the interior domains and
# exterior boundaries

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space and basis functions
V = FunctionSpace(mesh, "CG", 2)
phi = TrialFunction(V)  # potential
v  = TestFunction(V)

# Define Dirichlet boundary conditions at top and bottom boundaries
bcs = [DirichletBC(V, 0.0, boundaries, 1)]

# Define input data
a0 = Constant(1.0)
a1 = Constant(0.01)
f_water = Constant(waterInFlux)

# Define the reaction equation

# Define variational form
F_pot = inner(a0*grad(phi), grad(v))*dx - f_water*v*dx

# Separate left and right hand sides of equation
a, L = lhs(F_pot), rhs(F_pot)

# Solve problem
phi = Function(V)
solve(a == L, phi, bcs)

# Evaluate integral of normal gradient over top boundary
n = FacetNormal(mesh)
result_flux = assemble(dot(grad(phi), n)*ds(1))

# test conservation of water
expected_flux = -waterInFlux*assemble(1*dx)

print( "relative error of conservation of water %f" % ((result_flux-expected_flux)/expected_flux) )

# Save result
output = File("/tmp/res_potential_flow.pvd")
output<<phi





# since we know the velocity field now, we can continue to solve the reaction advection equation.

# Steady State  --- ONLY DIFFUSION_CONVECTION EQUATION JET!!!

#TODO Unstable for no diffusion!

clInFlux = 0.0
f_cl = Constant(clInFlux)

cl = TrialFunction(V) # Ca- concentration
na = TrialFunction(V) # Na+ concentration
k  = TrialFunction(V)  # K+ concentration

v_cl = TestFunction(V)
v_na = TestFunction(V)
v_k  = TestFunction(V)

n = FacetNormal(mesh)
F_steady = inner(0.01*grad(cl),grad(v_cl))*dx + inner(cl*grad(phi),grad(v_cl))*dx - inner(grad(phi),n)*v_cl*cl*ds - f_cl*v_cl*dx - 2.*v_cl*dx(1)

# Boundary conditions
bcs = [DirichletBC(V, 0.0, boundaries, 1)]

# Separate left and right hand sides of equation
a_1, L_1 = lhs(F_steady), rhs(F_steady)

# Solve problem
cl = Function(V)
solve(a_1 == L_1, cl, bcs)

# Save result
output = File("/tmp/res_cl.pvd")
output<<cl


