from dolfin import *

set_log_level(1)

# Sub domain for top (mark whole boundary; bottom, right, left will overwrite)
class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Sub domain for bottom        
class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS and on_boundary        

# Sub domain for right wall
class Right(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1.0 - DOLFIN_EPS and on_boundary

# Sub domain for left wall
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary

# Read mesh
mesh = Mesh("cube.xml")

# Create mesh functions over the cell facets
sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

# Mark all facets as sub domain 4
sub_domains.set_all(4)

# Mark top as sub domain 0
top = Top()
top.mark(sub_domains, 0)

# Mark bottom as sub domain 1
bottom = Bottom()
bottom.mark(sub_domains, 1)

# Mark right wall as sub domain 2
right = Right()
right.mark(sub_domains, 2)

# Mark left wall as sub domain 3
left = Left()
left.mark(sub_domains, 3)

# Save sub domains to file
file = File("cube_subdomains.xml")
file << sub_domains

# Save sub domains to VTK files
file = File("cube_subdomains.pvd")
file << sub_domains


