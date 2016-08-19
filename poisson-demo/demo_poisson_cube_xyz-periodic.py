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


class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        # return True if on left, front, or bottom boundary AND NOT on slave edges
        return bool((near(x[0], 0) or near(x[1], 0) or near(x[2], 0)) 
                    and (not ( (near(x[0], 1) and near(x[1], 0)) 
                        or (near(x[0], 1) and near(x[2], 0))
                        or (near(x[1], 1) and near(x[2], 0)) 
                        or (near(x[1], 1) and near(x[0], 0))
                        or (near(x[2], 1) and near(x[0], 0))
                        or (near(x[2], 1) and near(x[1], 0)))) 
                        and on_boundary)

    def map(self, x, y):
        if near(x[0], 1) and near(x[1], 1) and near(x[2], 1):
           y[0] = x[0] - 1.
           y[1] = x[1] - 1.
           y[2] = x[2] - 1.
        elif near(x[0], 1) and near(x[1], 1):
            y[0] = x[0] - 1.
            y[1] = x[1] - 1.
            y[2] = x[2]
        elif near(x[0], 1) and near(x[2], 1):
            y[0] = x[0] - 1.
            y[1] = x[1]
            y[2] = x[2] - 1.           
        elif near(x[1], 1) and near(x[2], 1):
            y[0] = x[0]
            y[1] = x[1] - 1.
            y[2] = x[2] - 1.           
        elif near(x[0], 1):
            y[0] = x[0] - 1.
            y[1] = x[1] 
            y[2] = x[2]
        elif near(x[1], 1):
            y[0] = x[0]
            y[1] = x[1] - 1.
            y[2] = x[2]            
        else: # near(x[2], 1):
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2] - 1.
        
pbc = PeriodicDomain()
mesh = UnitCubeMesh(10, 10, 10)
V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=pbc)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Expression("10*exp(-(pow(x[0] - 1.0, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
f = Expression("exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)+ pow(x[2] - 0.5, 2)) / 0.02)")
#f = Expression("0.0")
g = Expression("sin(10*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx #+ g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u)

# Save solution in VTK format
file = File("poisson_cube.pvd")
file << u

# Plot solution
plot(u, interactive=True)

