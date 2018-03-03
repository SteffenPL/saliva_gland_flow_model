#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: plunder
"""

from fenics import *

# Parameters
waterInFlux = 0.1


# Create classes for defining parts of the boundaries and the interior
# of the domain

class Hole(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0],1)

# Initialize sub-domain instances
hole = Hole()

# Define mesh
#mesh = UnitSquareMesh(64, 64)
mesh = Mesh('nontrivial_mesh_large.xml')


# Initialize mesh function for interior domains
domains = MeshFunction("size_t",mesh,2)
domains.set_all(0)
#obstacle.mark(domains, 1)

# Initialize mesh function for boundary domains
boundaries = MeshFunction("size_t", mesh,1)
boundaries.set_all(0)
hole.mark(boundaries, 1)

# Define function space and basis functions
V = FunctionSpace(mesh, "CG", 2)
u = TrialFunction(V)
v = TestFunction(V)

# Define Dirichlet boundary conditions at top and bottom boundaries
bcs = [DirichletBC(V, 5.0, boundaries, 1)]

# Define new measures associated with the interior domains and
# exterior boundaries

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define input data
a0 = Constant(1.0)
a1 = Constant(0.01)
f = Constant(waterInFlux)

# Define variational form
F = (inner(a0*grad(u), grad(v))*dx(0)
     - f*v*dx(0))

# Separate left and right hand sides of equation
a, L = lhs(F), rhs(F)

# Solve problem
u = Function(V)
solve(a == L, u, bcs)

# Evaluate integral of normal gradient over top boundary
n = FacetNormal(mesh)
result_flux = assemble(dot(grad(u), n)*ds(1))

# test conservation of water
expected_flux = -waterInFlux*assemble(1*dx)

print( "relative error of conservation of water %f" % ((result_flux-expected_flux)/expected_flux) )

output = File("/tmp/results.pvd")
output<<u
