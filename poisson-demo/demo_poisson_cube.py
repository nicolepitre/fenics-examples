"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0        for x = 0 or x = 1
du/dn(x, y) = sin(5*x) for y = 0 or y = 1
"""

# Copyright (C) 2007-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2007-08-16
# Last changed: 2012-11-12

# Begin demo

from dolfin import *

# Create mesh and define function space

mesh = Mesh("cube_periodic.xml")
subdomains = MeshFunction("size_t", mesh, "cube_periodic_facet_region.xml")
V = FunctionSpace(mesh, "Lagrange", 1)

# Define Dirichlet boundary (x = 0 or x = 1)
#def boundary(x):
#    return x[0] < DOLFIN_EPS -1.5 or x[0] > 1.5 - DOLFIN_EPS

# Define boundary condition
u0 = Constant(0.0)
u1 = Constant(.0)
# u1 = Expression("x[0]")

#bc1 = DirichletBC(V, u1, 1)
#bc2 = DirichletBC(V, u2, 2)
#bc3 = DirichletBC(V, u1, 3)

bc0 = DirichletBC(V, u0, subdomains, 111)
bc1 = DirichletBC(V, u0, subdomains, 333)
bcs = [bc0, bc1]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Expression("10*exp(-(pow(x[0] - 1.0, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
f = Expression("exp(-(pow(x[0] - 1.0, 2) + pow(x[1] - 0.5, 2)) / 0.02) - exp(-(pow(x[0] + 1.0, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
#f = Expression("0.0")
g = Expression("sin(10*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx #+ g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("cube_periodic_poisson.pvd")
file << u

# Plot solution
plot(u, interactive=True)
