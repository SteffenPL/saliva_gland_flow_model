#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 21:42:57 2018

@author: plunder
"""
from fenics import *
# import numpy as np
# import matplotlib.pyplot as plt

from video_data import VideoData

class DuctState:
    def __init__(self):

        self.mesh = None
        self.vid = None
        # parameters
        self.water_in_flux = 0.1

    def load_mesh(self, fn):
        self.mesh = Mesh(fn)

    def load_video_mesh(self, nx = 30, ny = 30):
        self.mesh = UnitSquareMesh(nx, ny)

    def load_video(self, fn):
        self.video_filename = fn

    def compute_potential_flow(self):

        class Hole(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 1) and 0.3 < x[1]  < 0.6

        #class ActiveArea(SubDomain):
        #    def inside(self, x, on_boundary):
        #        return between(x[0], (0., 0.1)) and between(x[1], (0, 0.2))

        # Initialize sub-domain instances
        hole = Hole()
        # self.activate = ActiveArea()

        # Initialize mesh function for interior domains
        self.domains = MeshFunction("size_t", self.mesh, 2)
        self.domains.set_all(0)
        # self.activate.mark(self.domains, 1)

        # Initialize mesh function for boundary domains
        self.boundaries = MeshFunction("size_t", self.mesh, 1)
        self.boundaries.set_all(0)
        hole.mark(self.boundaries, 1)

        # Define new measures associated with the interior domains and
        # exterior boundaries

        self.dx = Measure('dx', domain=self.mesh, subdomain_data=self.domains)
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        dx, ds = self.dx, self.ds

        # Define function space and basis functions
        V = FunctionSpace(self.mesh, "P", 1)
        self.V_vel = V
        phi = TrialFunction(V)  # potential
        v = TestFunction(V)

        # Define Dirichlet boundary conditions at top and bottom boundaries
        bcs = [DirichletBC(V, 0.0, self.boundaries, 1)]

        # Define input data
        a0 = Constant(1.0)
        a1 = Constant(0.01)
        f_water = Constant(self.water_in_flux)

        # Define the reaction equation

        # Define variational form
        F_pot = inner(a0 * grad(phi), grad(v)) * dx - f_water * v * dx

        # Separate left and right hand sides of equation
        a, L = lhs(F_pot), rhs(F_pot)

        # Solve problem
        phi = Function(V)
        solve(a == L, phi, bcs)

        # Evaluate integral of normal gradient over top boundary
        n = FacetNormal(self.mesh)
        result_flux = assemble(dot(grad(phi), n) * ds(1))

        # test conservation of water
        expected_flux = -self.water_in_flux * assemble(1 * dx)

        print("relative error of conservation of water %f" % ((result_flux - expected_flux) / expected_flux))

        self.phi = phi
        self.flow = -grad(phi)

        # Save result
        output = File("/tmp/potential.pvd")
        output << self.phi

    def F_diff_conv(self, u, v, n, flow, D, s, source):

        dx, ds = self.dx, self.ds
        return D * inner(grad(u), grad(v)) * dx \
               + s * inner(flow, grad(u)) * v * dx \
               - source * v * dx

    def compute_steady_state(self):

        names = {'Cl', 'Na', 'K'}

        P1 = FiniteElement('P', triangle, 1)
        element = MixedElement([P1, P1, P1])
        V = FunctionSpace(self.mesh, element)
        self.V_conc = V

        (u_cl, u_na, u_k) = TrialFunction(V)
        (v_cl, v_na, v_k) = TestFunction(V)

        assert (self.flow is not None)

        n = FacetNormal(self.mesh)

        # F = ( self.F_diff_conv(u_cl, v_cl, n, grad(self.phi), 1. ,1., 0.)
        #    + self.F_diff_conv(u_na, v_na, n, grad(self.phi), 1. ,1., 0.)
        #    + self.F_diff_conv(u_k , v_k , n, grad(self.phi), 1. ,1., 0.) )

        dx, ds = self.dx, self.ds
        flow = self.flow
        F = inner(grad(u_cl), grad(v_cl)) * dx \
            + inner(flow, grad(u_cl)) * v_cl * dx \
            + inner(grad(u_na), grad(v_na)) * dx \
            + inner(flow, grad(u_na)) * v_na * dx \
            + inner(grad(u_k), grad(v_k)) * dx \
            + inner(flow, grad(u_k)) * v_k * dx

        a, L = lhs(F), rhs(F)
        a_mat = assemble(a)
        L_vec = assemble(L)
        # solve
        u = Function(V)
        solve(a_mat, u.vector(), L_vec)

        u_cl, u_na, u_k = u.split()

        output1 = File('/tmp/steady_state_cl.pvd')
        output1 << u_cl
        output2 = File('/tmp/steady_state_na.pvd')
        output2 << u_na
        output3 = File('/tmp/steady_state_k.pvd')
        output3 << u_k

        self.u_cl = u_cl
        self.u_na = u_na
        self.u_k = u_k

    def compute_conv_diff_reac(self, initial_condition=None):

        names = {'Cl', 'Na', 'K'}

        dt = 0.1
        t = 0.
        t_end = 1.

        P1 = FiniteElement('P', triangle, 1)
        element = MixedElement([P1, P1, P1])
        V = FunctionSpace(self.mesh, element)
        self.V_conc = V

        u_init = Function(V)
        (u_cl, u_na, u_k) = TrialFunction(V)
        (v_cl, v_na, v_k) = TestFunction(V)

        if initial_condition is None:
            initial_condition = Expression(("exp(-((x[0]-0.1)*(x[0]-0.1)+x[1]*x[1])/0.01)",
                                            "exp(-((x[0]-0.12)*(x[0]-0.12)+x[1]*x[1])/0.01)",
                                            "0."), element=element)

        u_init = interpolate(initial_condition, V)
        u_init_cl = u_init[0]
        u_init_na = u_init[1]
        u_init_k = u_init[2]

        assert (self.flow is not None)

        n = FacetNormal(self.mesh)

        dx, ds = self.dx, self.ds
        flow = 10 * self.flow
        f_in = Constant(0.00)
        D = Constant(0.01)
        k1 = Constant(0.1)
        k_1 = Constant(0.001)
        F = (
                (u_cl - u_init_cl) * v_cl * dx
                + dt * D * inner(grad(u_cl), grad(v_cl)) * dx
                + dt * inner(flow, grad(u_cl)) * v_cl * dx
                + (u_na - u_init_na) * v_na * dx
                + dt * D * inner(grad(u_na), grad(v_na)) * dx
                + dt * inner(flow, grad(u_na)) * v_na * dx
                + (u_k - u_init_k) * v_k * dx
                + dt * D * inner(grad(u_k), grad(v_k)) * dx
                + dt * inner(flow, grad(u_k)) * v_k * dx
                + f_in * v_cl * dx
                + f_in * v_na * dx
                + f_in * v_k * dx
                + dt * k1 * u_init_cl * u_init_na * v_cl * dx
                + dt * k1 * u_init_cl * u_init_na * v_na * dx
                - dt * k1 * u_init_cl * u_init_na * v_k * dx
                - dt * k_1 * u_init_k * v_cl * dx
                - dt * k_1 * u_init_k * v_na * dx
                + dt * k_1 * u_init_k * v_k * dx
        )

        self.F = F

        a, L = lhs(F), rhs(F)
        a_mat = assemble(a)
        L_vec = assemble(L)

        output1 = File('/tmp/cl_dyn.pvd')
        output2 = File('/tmp/na_dyn.pvd')
        output3 = File('/tmp/k_dyn.pvd')
        output4 = File('/tmp/all_dyn.pvd')
        # solve

        self.sol = []

        while t < t_end:
            t = t + dt
            print(t)

            u = Function(V)

            a_mat = assemble(a)
            L_vec = assemble(L)
            solve(a_mat, u.vector(), L_vec)
            # NonlinearVariationalProblem(F,u)
            u_init.assign(u)

            u_cl, u_na, u_k = u.split()

            # u_init_cl.assign(u_cl)
            # u_init_na.assign(u_na)
            # u_init_k.assign(u_k)

            u_cl.rename("cl", "cl")
            u_na.rename("na", "na")
            u_k.rename("k", "k")
            output1 << u_cl, t
            output2 << u_na, t
            output3 << u_k, t
            self.sol.append((u_cl, u_na, u_k))

        plot(u_k)

        self.u_cl = u_cl
        self.u_na = u_na
        self.u_k = u_k


    def compute_conv_diff_reac_video(self, initial_condition=None):

        names = {'Cl', 'K'}


        P1 = FiniteElement('P', triangle, 1)
        element = MixedElement([P1, P1])
        V_single = FunctionSpace(self.mesh, P1)
        V = FunctionSpace(self.mesh, element)
        self.V_conc = V

        # load video
        video = VideoData(element=P1)
        video.load_video(self.video_filename)

        print(video)

        dt = 0.1
        t = 0.
        t_end = 10.

        u_init = Function(V)
        (u_cl, u_k) = TrialFunction(V)
        (v_cl, v_k) = TestFunction(V)

        #u_na = video


        if initial_condition is None:
            initial_condition = Expression(("f*exp(-0.5*((x[0]-a)*(x[0]-a)+(x[1]-b)*(x[1]-b))/var)/(sqrt(2*pi)*var)",
                                            "0."), a = 0.2, b= 0.5, var=0.1, f=0.01, pi=pi, element=element)

        u_init = interpolate(initial_condition, V)
        u_init_cl = u_init[0]

        u_init_k = u_init[1]
        u_init_na : VideoData= video

        assert (self.flow is not None)

        n = FacetNormal(self.mesh)

        dx, ds = self.dx, self.ds
        flow = 10 * self.flow
        f_in = Constant(0.00)
        D = Constant(0.01)

        C_na = Constant(0.01)
        k1 = Constant(0.1)
        k_1 = Constant(0.001)
        F = (
                (u_cl - u_init_cl) * v_cl * dx
                + dt * D * inner(grad(u_cl), grad(v_cl)) * dx
                + dt * inner(flow, grad(u_cl)) * v_cl * dx
                #+ (u_na - u_init_na) * v_na * dx
                #+ dt * D * inner(grad(u_na), grad(v_na)) * dx
                #+ dt * inner(flow, grad(u_na)) * v_na * dx
                + (u_k - u_init_k) * v_k * dx
                + dt * D * inner(grad(u_k), grad(v_k)) * dx
                + dt * inner(flow, grad(u_k)) * v_k * dx
                + f_in * v_cl * dx
                #+ f_in * v_na * dx
                + f_in * v_k * dx
                + dt * k1 * u_init_cl * C_na * u_init_na * v_cl * dx
                #+ dt * k1 * u_init_cl * u_init_na * v_na * dx
                - dt * k1 * u_init_cl * C_na * u_init_na * v_k * dx
                - dt * k_1 * u_init_k * v_cl * dx
                #- dt * k_1 * u_init_k * v_na * dx
                + dt * k_1 * u_init_k * v_k * dx
        )

        self.F = F

        a, L = lhs(F), rhs(F)
        a_mat = assemble(a)
        L_vec = assemble(L)

        output1 = File('/tmp/cl_dyn.pvd')
        output2 = File('/tmp/na_dyn.pvd')
        output3 = File('/tmp/k_dyn.pvd')
        output4 = File('/tmp/all_dyn.pvd')
        # solve

        self.sol = []

        u_na = Function(V_single)
        u_na = C_na*interpolate(u_init_na,V_single)

        na_inflow = 0

        while t < t_end:
            t = t + dt
            print(t)

            u = Function(V)

            u_init_na.set_time(10*t)

            a_mat = assemble(a)
            L_vec = assemble(L)
            solve(a_mat, u.vector(), L_vec)
            # NonlinearVariationalProblem(F,u)
            u_init.assign(u)

            u_cl, u_k = u.split()

            # u_init_cl.assign(u_cl)
            # u_init_na.assign(u_na)
            # u_init_k.assign(u_k)

            u_na = interpolate(u_init_na, V_single)

            u_cl.rename("cl", "cl")
            u_na.rename("na", "na")
            u_k.rename("k", "k")
            output1 << u_cl, t
            output2 << u_na, t
            output3 << u_k, t
            self.sol.append((u_cl, u_na, u_k))

        plot(u_k)

        self.u_cl = u_cl
        #self.u_na = u_na
        self.u_k = u_k


    def check_conservation_law(self, sol):
        dx = self.dx

        total = []
        out = []
        for s in sol:
            cl = s[0]
            na = s[1]
            k = s[2]

            # integrate interior concentration
            mass_cl = assemble(cl * dx)
            mass_na = assemble(na * dx)
            mass_k = assemble(k * dx)

            # integrate outflux
            n = FacetNormal(self.mesh)
            flow = self.flow
            outflux_cl = assemble(cl * inner(flow, n) * ds)
            outflux_na = assemble(na * inner(flow, n) * ds)
            outflux_k = assemble(k * inner(flow, n) * ds)

            total.append(mass_cl + mass_na + 2 * mass_k)
            out.append(outflux_cl + outflux_na + 2 * outflux_k)

        self.total_mass = total
        self.outflux = out


d = DuctState()


def main():
    d.load_mesh('nontrivial_mesh_med.xml')

    d.compute_potential_flow()
    d.compute_conv_diff_reac()


def demo_video():
    d.load_video_mesh(20, 20)
    d.compute_potential_flow()
    d.load_video("video_data/ach12.mp4")

    d.compute_conv_diff_reac_video()


def demo_video_complex():
    d.load_video_mesh(20, 20)
    d.compute_potential_flow()
    d.load_video("video_data/ach12.mp4")

    d.compute_conv_diff_reac_video()


if __name__ == "__main__":
    # main()

    demo_video_complex()
