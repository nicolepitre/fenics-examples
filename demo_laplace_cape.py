"""This demo program solves Laplace's equation

    - div grad u(x, y) = 0

on the cape mesh with four Dirichlet boundary conditions.
"""
# Begin demo

from dolfin import *

# Create mesh and define function space
mesh = Mesh("cape.xml")
subdomains = MeshFunction("size_t", mesh, "cape_facet_region.xml")
V = FunctionSpace(mesh, "Lagrange", 1)

# Define boundary condition
u1 = Expression("x[1] + 1")
u2 = Expression("x[0]")
u3 = Expression("x[1]")
u4 = Expression("x[0] + x[1]")

bc1 = DirichletBC(V, u1, subdomains, 1)
bc2 = DirichletBC(V, u2, subdomains, 2)
bc3 = DirichletBC(V, u3, subdomains, 3)
bc4 = DirichletBC(V, u4, subdomains, 4)

bcs = [bc1, bc2, bc3, bc4]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
#f = Expression("10*exp(-(pow(x[0] - 1.0, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
f = Expression("0.0")
g = Expression("sin(10*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx #+ g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("laplace_cape.pvd")
file << u

# Plot solution
plot(u, interactive=True)

