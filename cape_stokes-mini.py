"""This demo solves the Stokes equations, using linear elements
enriched with a bubble for the velocity and linear elements for the
pressure (Mini elements). The sub domains for the different boundary
conditions used in this simulation are computed by the demo program in
src/demo/mesh/subdomains."""

# Copyright (C) 2007 Kristian B. Oelgaard
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
# Modified by Anders Logg, 2008-2009.
#
# First added:  2007-11-16
# Last changed: 2010-04-01
# Begin demo

from __future__ import print_function
from dolfin import *

# Load mesh and subdomains
mesh = Mesh("cape.xml")
#sub_domains = MeshFunction("size_t", mesh, "cape_subdomains.xml")
sub_domains = MeshFunction("size_t", mesh, "cape_facet_region.xml")


# Define function spaces
P1 = VectorFunctionSpace(mesh, "Lagrange", 1)
B  = VectorFunctionSpace(mesh, "Bubble", 3)
Q  = FunctionSpace(mesh, "CG",  1)
V = P1 + B
Mini = V*Q

# Top boundary condition for velocity
# NOTE: Projection here is inefficient workaround of issue #489, FFC issue #69
top = project(Constant((0, 0)), V)
bc0 = DirichletBC(Mini.sub(0), top, sub_domains, 4)

# Right boundary condition for velocity
# NOTE: Projection here is inefficient workaround of issue #489, FFC issue #69
right = project(Expression(("-sin(x[1]*pi)", "0.0")), V)
bc1 = DirichletBC(Mini.sub(0), right, sub_domains, 1)

# Collect boundary conditions
bcs = [bc0, bc1]

# Define variational problem
(u, p) = TrialFunctions(Mini)
(v, q) = TestFunctions(Mini)
f = Constant((0, 0))
a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = inner(f, v)*dx

# Compute solution
w = Function(Mini)
solve(a == L, w, bcs)

# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)

print("Norm of velocity coefficient vector: %.15g" % u.vector().norm("l2"))
print("Norm of pressure coefficient vector: %.15g" % p.vector().norm("l2"))

# Split the mixed solution using a shallow copy
(u, p) = w.split()

# Save solution in VTK format
ufile_pvd = File("cape_velocity.pvd")
ufile_pvd << u
pfile_pvd = File("cape_pressure.pvd")
pfile_pvd << p

# Plot solution
plot(u)
plot(p)
interactive()
