This is a quick implementation of a simple Finite Element (FE) solver using 
Eigen for the computations and Qt the for the graphical user interface (GUI).

This code solves the magnetostatic Poisson problem with the Finite Element 
Method on a two-dimensional triangular mesh. The mesh file is imported from 
Gmsh. The user defines the material parameters and excitation for each 
Physical region using the GUI. Zero Dirichlet condition is assumed on all 
Physical Lines.

The GUI visualizes the solutions with contour plots. Since the main purpose
of the code (for the author) is to play with visualizations, the solution is 
recomputed every time the material parameters is changed.

Technical details:

Mesh files generated with GMsh are imported with mesh.cc, mesh_element.cc,
mesh_file.cc and mesh.cc. Material parameters are specified with Region-
objects and assembled into a map according to the "physical numbers" 
(see region.cc and region.h). The element stiffness and mass matrices for 
first order basis functions are computed with Gaussian quadrature and 
assembled into global stiffness and mass matrices in element.cc and 
assembly.cc. 

Sparse LU solver of Eigen is used for the solution of the resulting system
of equations. This code is applicable to large problems even with ~1M 
degrees of freedom.

The GUI is implemented with Qt in meshplot.cpp, physlist.cpp, mainwindow.cpp 
and main.cpp:

- Files meshplot.cpp and meshplot.h implement a QWidget for visualizing 
  geometric primitives such as lines and polygons with QPainter.
- Files physlist.cpp and physlist.h implement a QWidget based on QTableWidget
  for the mainpulation of material parameters and source currents.

