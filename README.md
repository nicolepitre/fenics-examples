# fenics-examples
A collection of FEniCS examples that show simple ways to use FEniCS to solve PDEs. Written in Python.

This repository contains examples of Python code for FEniCS. It also contains geometry files written in the Gmsh scripting language. All files with the .geo extension can be displayed in Gmsh. The corresponding mesh files (.msh extension) are contained in the folder called meshes. The XML versions of the meshes are contained in the folder called xml-meshes. Any file with the .py extension should be run inside of FEniCS, including the dolfin-convert script.

The following tips and notes might be useful when learning FEniCS and Gmsh.

## FEniCS Notes
- List of FEniCS demos: http://fenicsproject.org/documentation/dolfin/1.6.0/python/demo/index.html
- Meshes used in FEniCS demos (XML format): http://fenicsproject.org/download/data.html#data
- What is Docker?
    -  https://www.dataquest.io/blog/docker-data-science/
- The dolfin-convert script can be found here: https://people.sc.fsu.edu/~jburkardt/py_src/dolfin-convert/dolfin-convert.html
- It only takes two commands to get a Jupyter notebook up and running in Docker with a FEniCS environment installed
  - The commands are:
    - fenicsproject notebook myproject
    - fenicsproject start myproject 
    - To close the notebook, type the following command: docker stop myproject
- Periodic boundary condition class in FEniCS: 
    - https://bitbucket.org/fenics-project/dolfin/pull-requests/4/adding-test-for/diff 
    - https://fenicsproject.org/documentation/dolfin/1.3.0/python/demo/documented/periodic/python/documentation.html
    - https://fenicsproject.org/qa/262/possible-specify-more-than-one-periodic-boundary-condition

## Gmsh Notes
- Tutorial on how to use Gmsh: https://openfoamwiki.net/index.php/2D_Mesh_Tutorial_using_GMSH
- Command for Mac OS X Terminal to set up the gmsh command line: 
  - ln -s //Applications/Gmsh.app/Contents/MacOS/gmsh /usr/local/bin/gmsh
- Command to make cape.geo file into a 2D mesh file:
  - gmsh -2 cape.geo -o cape.msh
- List of Gmsh command-line options: *Section 3.3* of http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options
- Gmsh Periodic Line/Surface commands: http://www.manpagez.com/info/gmsh/gmsh-2.7.1/gmsh_47.php

Once you are comfortable with FEniCS and Gmsh, you can create your own demos by following the steps below.

## Steps to create a FEniCS Demo for a particular geometry
1. Create a geometry file in the Gmsh script language.
2. Use Gmsh to generate a finite element mesh from the geometry file.
3. Convert the mesh into XML format using the dolfin-convert script. The script also produces an XML file called filename_facet_regions.xml if you have marked any physical lines/surfaces in Gmsh.
4. Write a Python demo script using the FEniCS DOLFIN module.
5. Add the XML mesh file into the Python demo script.
6. Open FEniCS in a virtual environment such as Docker and run the demo.
7. View the solution in VisIt.

