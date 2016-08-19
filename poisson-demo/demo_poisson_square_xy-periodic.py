"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)

"""

# Begin demo

from dolfin import *

class PeriodicBoundary(SubDomain):

    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 1) and (1, 0)
        return bool((near(x[0], 0) or near(x[1], 0)) and 
                (not ((near(x[0], 0) and near(x[1], 1)) 
                        or 
                        (near(x[0], 1) and near(x[1], 0)))
                        ) and on_boundary)

    def map(self, x, y):
        if near(x[0], 1) and near(x[1], 1):
            y[0] = x[0] - 1.
            y[1] = x[1] - 1.
        elif near(x[0], 1):
            y[0] = x[0] - 1.
            y[1] = x[1]
        else:   # near(x[1], 1)
            y[0] = x[0]
            y[1] = x[1] - 1.
            

# Create mesh and define function space
mesh = Mesh("square.xml")
subdomains = MeshFunction("size_t", mesh, "square_facet_region.xml")
pbc = PeriodicBoundary()
V = FunctionSpace(mesh, "Lagrange", 1, constrained_domain=pbc)

# boundary conditions
u0 = Constant(0.0)
bc0 = DirichletBC(V, u0, subdomains, 2)
bc1 = DirichletBC(V, u0, subdomains, 4)
bcs = [bc0, bc1]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2))/ 0.02)")
a = inner(grad(u), grad(v))*dx
L = f*v*dx

# Compute solution
u = Function(V)
solve(a == L, u, bcs)

# Save solution in VTK format
file = File("poisson_square.pvd")
file << u

# Plot solution
plot(u, interactive=True)
