"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

and boundary conditions given by

    u(x, y) = 0        for x = 0 or x = 1
du/dn(x, y) = sin(5*x) for y = 0 or y = 1
"""

# Begin demo

from dolfin import *

class PeriodicBoundaryX(SubDomain):
    def __init__(self, tolerance=DOLFIN_EPS):
        SubDomain.__init__(self, tolerance)
        self.tol = tolerance
    def inside(self, x, on_boundary):
        return bool(x[0] < self.tol and x[0] > -self.tol and on_boundary)
    def map(self, x, y):
        y[0] = x[0] - 1.0
        y[1] = x[1]

# Create mesh and define function space
mesh = Mesh("cape.xml")
subdomains = MeshFunction("size_t", mesh, "cape_facet_region.xml")
pbc = PeriodicBoundaryX()
V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=pbc)

# boundary conditions
u0 = Constant(0.0)
bc0 = DirichletBC(V, u0, subdomains, 4)
bc1 = DirichletBC(V, u0, subdomains, 2)
bcs = [bc0, bc1]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2))/ 0.02)")
#f = Expression("0.0")
g = Expression("sin(5*x[0])")
a = inner(grad(u), grad(v))*dx
L = f*v*dx #+ g*v*ds

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("poisson_cape.pvd")
file << u

# Plot solution
plot(u, interactive=True)

