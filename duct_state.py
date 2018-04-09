#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 21:42:57 2018

@author: plunder
"""
import fenics as fe
from fenics import Mesh, near, between, UnitSquareMesh, SubDomain
from fenics import div, grad, inner, dot
from fenics import FiniteElement, FunctionSpace, VectorFunctionSpace
from fenics import MeshFunction, TrialFunction, TestFunction, Function, Expression
from fenics import Measure
# import numpy as np
import matplotlib.pyplot as plt

from video_data import VideoData

class DuctState:
    def __init__(self):

        self.mesh = None
        self.vid = None
        # parameters
        self.water_in_flux = 0.1

        self.boundary_fnc = lambda x: near(x[0], 1) and 0.3 < x[1]  < 0.6

    def load_mesh(self, fn):
        self.mesh = Mesh(fn)

    def load_video_mesh(self, nx = 30, ny = 30):
        self.mesh = UnitSquareMesh(nx, ny)

    def load_video(self, fn):
        self.video_filename = fn



    def compute_potential_flow(self):

        class Hole(SubDomain):
            def inside(s, x, on_boundary):
                return on_boundary and self.boundary_fnc(x)

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
        bcs = [fe.DirichletBC(V, 0.0, self.boundaries, 1)]

        # Define input data
        a0 = fe.Constant(1.0)
        a1 = fe.Constant(0.01)
        f_water = fe.Constant(self.water_in_flux)

        # Define the reaction equation

        # Define variational form
        F_pot = inner(a0 * grad(phi), grad(v)) * dx - f_water * v * dx

        # Separate left and right hand sides of equation
        a, L = fe.lhs(F_pot), fe.rhs(F_pot)

        # Solve problem
        phi = fe.Function(V)
        fe.solve(a == L, phi, bcs)

        # Evaluate integral of normal gradient over top boundary
        n = fe.FacetNormal(self.mesh)
        result_flux = fe.assemble(dot(grad(phi), n) * ds(1))

        # test conservation of water
        expected_flux = -self.water_in_flux * fe.assemble(1 * dx)

        print("relative error of conservation of water %f" % ((result_flux - expected_flux) / expected_flux))

        self.phi = phi
        self.flow = -grad(phi)

        # Save result
        output = fe.File("/tmp/potential.pvd")
        output << self.phi




    def F_diff_conv(self, u, v, n, flow, D, s, source):

        dx, ds = self.dx, self.ds
        return D * inner(grad(u), grad(v)) * dx \
               + s * inner(flow, grad(u)) * v * dx \
               - source * v * dx

    def compute_steady_state(self):

        names = {'Cl', 'Na', 'K'}

        P1 = FiniteElement('P', fe.triangle, 1)
        element = MixedElement([P1, P1, P1])
        V = FunctionSpace(self.mesh, element)
        self.V_conc = V

        (u_cl, u_na, u_k) = TrialFunction(V)
        (v_cl, v_na, v_k) = TestFunction(V)

        assert (self.flow is not None)

        n = fe.FacetNormal(self.mesh)

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

        a, L = fe.lhs(F), fe.rhs(F)
        a_mat = fe.assemble(a)
        L_vec = fe.assemble(L)
        # solve
        u = Function(V)
        fe.solve(a_mat, u.vector(), L_vec)

        u_cl, u_na, u_k = u.split()

        output1 = fe.File('/tmp/steady_state_cl.pvd')
        output1 << u_cl
        output2 = fe.File('/tmp/steady_state_na.pvd')
        output2 << u_na
        output3 = fe.File('/tmp/steady_state_k.pvd')
        output3 << u_k

        self.u_cl = u_cl
        self.u_na = u_na
        self.u_k = u_k

    def compute_conv_diff_reac(self, initial_condition=None):

        names = {'Cl', 'Na', 'K'}

        dt = 0.1
        t = 0.
        t_end = 1.

        P1 = FiniteElement('P', fe.triangle, 3)
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

        u_init = fe.interpolate(initial_condition, V)
        u_init_cl = u_init[0]
        u_init_na = u_init[1]
        u_init_k = u_init[2]

        assert (self.flow is not None)

        n = fe.FacetNormal(self.mesh)

        dx, ds = self.dx, self.ds
        flow = 10 * self.flow
        f_in = fe.Constant(0.00)
        D = fe.Constant(0.01)
        k1 = fe.Constant(0.1)
        k_1 = fe.Constant(0.001)
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

        a, L = fe.lhs(F), fe.rhs(F)
        a_mat = fe.assemble(a)
        L_vec = fe.assemble(L)

        output1 = fe.File('/tmp/cl_dyn.pvd')
        output2 = fe.File('/tmp/na_dyn.pvd')
        output3 = fe.File('/tmp/k_dyn.pvd')
        output4 = fe.File('/tmp/all_dyn.pvd')
        # solve

        self.sol = []

        while t < t_end:
            t = t + dt
            print(t)

            u = Function(V)

            a_mat = fe.assemble(a)
            L_vec = fe.assemble(L)
            fe.solve(a_mat, u.vector(), L_vec)
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



        self.u_cl = u_cl
        self.u_na = u_na
        self.u_k = u_k


    def compute_conv_diff_reac_video(self, initial_condition=None, video_ref=None, video_size=None):

        names = {'Cl', 'K'}


        P1 = FiniteElement('P', fe.triangle, 1)
        element = fe.MixedElement([P1, P1])
        V_single = FunctionSpace(self.mesh, P1)
        V = FunctionSpace(self.mesh, element)
        self.V_conc = V

        # load video
        video = VideoData(element=P1)
        video.load_video(self.video_filename)

        if( video_ref is not None and video_size is not None):
            video.set_reference_frame(video_ref, video_size)


        print(video)

        dt = 0.05
        t = 0.
        t_end = 20.

        u_init = Function(V)
        (u_cl, u_k) = TrialFunction(V)
        (v_cl, v_k) = TestFunction(V)

        #u_na = video


        if initial_condition is None:
            #initial_condition = Expression(("f*exp(-0.5*((x[0]-a)*(x[0]-a)+(x[1]-b)*(x[1]-b))/var)/(sqrt(2*pi)*var)",
            #                                "0."), a = 80, b= 55, var=10, f=10, pi=fe.pi, element=element)

            initial_condition = Expression(("f",
                                            "0."), a=80, b=55, var=0.1, f=0.1, pi=fe.pi, element=element)

        u_init = fe.interpolate(initial_condition, V)
        u_init_cl = u_init[0]

        u_init_k = u_init[1]
        u_init_na : VideoData= video

        assert (self.flow is not None)

        n = fe.FacetNormal(self.mesh)

        dx, ds = self.dx, self.ds
        flow = 5. * self.flow
        f_in = fe.Constant(0.00)
        f_in_cl = fe.Constant(-0.05)
        D = fe.Constant(0.1)

        C_na = fe.Constant(0.1)
        k1 = fe.Constant(0.2)
        k_1 = fe.Constant(0.00001)
        # explicit
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
                + f_in_cl * v_cl * dx
                #+ f_in * v_na * dx
                + f_in * v_k * dx
                + dt * k1 * u_init_cl * C_na * u_init_na * v_cl * dx
                #+ dt * k1 * u_init_cl * u_init_na * v_na * dx
                - dt * k1 * u_init_cl * C_na * u_init_na * v_k * dx
                - dt * k_1 * u_init_k * v_cl * dx
                #- dt * k_1 * u_init_k * v_na * dx
                + dt * k_1 * u_init_k * v_k * dx
        )
        # implicit
        F = (
                (u_cl - u_init_cl) * v_cl * dx
                + dt * D * inner(grad(u_cl), grad(v_cl)) * dx
                + dt * inner(flow, grad(u_cl)) * v_cl * dx
                # + (u_na - u_init_na) * v_na * dx
                # + dt * D * inner(grad(u_na), grad(v_na)) * dx
                # + dt * inner(flow, grad(u_na)) * v_na * dx
                + (u_k - u_init_k) * v_k * dx
                + dt * D * inner(grad(u_k), grad(v_k)) * dx
                + dt * inner(flow, grad(u_k)) * v_k * dx
                + f_in_cl * v_cl * dx
                # + f_in * v_na * dx
                + f_in * v_k * dx
                + dt * k1 * u_cl * C_na * u_init_na * v_cl * dx
                # + dt * k1 * u_init_cl * u_init_na * v_na * dx
                - dt * k1 * u_cl * C_na * u_init_na * v_k * dx
                - dt * k_1 * u_k * v_cl * dx
                # - dt * k_1 * u_init_k * v_na * dx
                + dt * k_1 * u_k * v_k * dx
        )

        self.F = F

        a, L = fe.lhs(F), fe.rhs(F)
        a_mat = fe.assemble(a)
        L_vec = fe.assemble(L)

        output1 = fe.File('/tmp/cl_dyn.pvd')
        output2 = fe.File('/tmp/na_dyn.pvd')
        output3 = fe.File('/tmp/k_dyn.pvd')
        output4 = fe.File('/tmp/all_dyn.pvd')
        # solve

        self.sol = []

        u_na = Function(V_single)
        u_na = fe.interpolate(u_init_na,V_single)

        na_inflow = 0

        t_plot = 0.5
        t_last_plot = 0

        while t < t_end:
            t = t + dt
            t_last_plot += dt
            print(t)

            u = Function(V)

            u_init_na.set_time(5*t)

            a_mat = fe.assemble(a)
            L_vec = fe.assemble(L)
            fe.solve(a_mat, u.vector(), L_vec)
            # NonlinearVariationalProblem(F,u)
            u_init.assign(u)

            u_cl, u_k = u.split()

            # u_init_cl.assign(u_cl)
            # u_init_na.assign(u_na)
            # u_init_k.assign(u_k)

            u_na = fe.interpolate(u_init_na, V_single)

            u_cl.rename("cl", "cl")
            u_na.rename("na", "na")
            u_k.rename("k", "k")
            output1 << u_cl, t
            output2 << u_na, t
            output3 << u_k, t
            self.sol.append((u_cl, u_na, u_k))

            print( fe.assemble(u_cl*self.dx))

            if t_last_plot > t_plot:
                t_last_plot = 0
                plt.figure(figsize=(8,16))
                plt.subplot(211)
                fe.plot(u_cl)
                plt.subplot(212)
                fe.plot(u_k)
                plt.show()



        self.u_cl = u_cl
        #self.u_na = u_na
        self.u_k = u_k


    def check_conservation_law(self, sol):
        dx = self.dx
        ds = self.ds
        total = []
        out = []
        for s in sol:
            cl = s[0]
            na = s[1]
            k = s[2]

            # integrate interior concentration
            mass_cl = fe.assemble(cl * dx)
            mass_na = fe.assemble(na * dx)
            mass_k = fe.assemble(k * dx)

            # integrate outflux
            n = fe.FacetNormal(self.mesh)
            flow = self.flow
            outflux_cl = fe.assemble(cl * inner(flow, n) * ds)
            outflux_na = fe.assemble(na * inner(flow, n) * ds)
            outflux_k = fe.assemble(k * inner(flow, n) * ds)

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


from gmsh_interface import GmshInterface

def demo_video_complex():

    #comm_mpi4py = fe.mpi_comm_world().tompi4py()

    #rk = comm_mpi4py.Get_rank()
    #if rk == 0:
    g = GmshInterface("geo/movie_2.geo", dim=2)
    g.generate_xml("geo/movie_2.xml",lc=1.,replace=False)

    #comm_mpi4py.Barrier()

    d.load_mesh("geo/movie_2.xml")
    d.boundary_fnc = lambda x: x[0] > -5. and x[1] < 5.
    fe.plot(d.mesh)
    plt.show()


    d.compute_potential_flow()

    fe.plot(d.flow)
    plt.show()
    d.load_video("video_data/ach12.mp4")

    d.compute_conv_diff_reac_video(video_ref=[0,0],video_size=[99.,69.])


if __name__ == "__main__":

    # main()
    demo_video_complex()
