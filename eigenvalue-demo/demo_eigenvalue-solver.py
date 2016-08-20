from dolfin import *                                                            
import scipy.sparse as sp                                                       
from scipy.sparse.linalg import eigs                                            

mesh = UnitSquareMesh(10,10)                                                    

V = FunctionSpace(mesh, "Lagrange", degree = 1)                                 
bc = DirichletBC(V, 0.0, 'on_boundary')                                         
u = TrialFunction(V)                                                            
v = TestFunction(V)                                                             

a = inner(grad(u), grad(v))*dx                                                  
m = u * v * dx                                                                  

L = inner(Constant(1), v)*dx                                                    

A = PETScMatrix()                                                               
assemble(a, tensor=A)                                                           
bc.apply(A)                                                                     

M = PETScMatrix()                                                               
assemble(m, tensor=M)                                                           

# This just converts PETSc to CSR                                               
A = sp.csr_matrix(A.mat().getValuesCSR()[::-1])                                 
M = sp.csr_matrix(M.mat().getValuesCSR()[::-1])                                 

v, V = eigs(A, 6, M, sigma=1.5)                                                 

print v
