//---------------------------------------------------
// Main page documentation
/**
   \mainpage libMesh - A C++ Finite Element Library

   The \p libMesh library is a C++ framework for the numerical
   simulation of partial differential equations on serial and parallel
   platforms.  Development began in March 2002 with the intent of
   providing a friendly interface to a number of high-quality software
   packages that are currently available.

   A major goal of the library is to provide support for adaptive mesh
   refinement (AMR) computations in parallel while allowing a research
   scientist to focus on the physics they are modeling.  The library
   currently offers:

   - Partitioning Algorithms
   - Metis K-Way weighted graph partitioning
   - Parmetis parallel graph partitioning
   - Hilbert and Morton-ordered space filling curves

   - Generic 2D Finite Elements
   - 3 and 6 noded triangles (\p Tri3, \p Tri6)
   - 4, 8, and 9 noded quadrilaterals (\p Quad4, \p Quad8, \p Quad9)
   - 4 and 6 noded infinite quadrilaterals (\p InfQuad4, \p InfQuad6)

   - Generic 3D Finite Elements
   - 4 and 10 noded tetrahedrals (\p Tet4, \p Tet10)
   - 8, 20, and 27 noded hexahedrals (\p Hex8, \p Hex20, \p Hex27)
   - 6, 15, and 18 noded prisms (\p Prism6, \p Prism15, \p Prism18)
   - 5 noded pyramids (\p Pyramid5)
   - 8, 16, and 18 noded infinite hexahedrals (\p InfHex8,
   \p InfHex16, \p InfHex18)
   - 6 and 12 noded infinite prisms (\p InfPrism6, \p InfPrism12)

   - Generic Finite Element Families
   - Lagrange
   - Hierarchic
   - Discontinuous Monomials

   - Dimension-independence
   - Operators are defined to allow the same code
   to run unmodified on 2D and 3D applications
   - The code you debug and verify on small 2D problems
   can immediately be applied to large, parallel 3D applications

   - Sparse Linear Algebra
   - \p PETSc provides a suite of iterative solvers and preconditioners
   for serial and parallel applications
   - Complex values are supported with \p PETSc
   - \p LASPACK provides iterative solvers and preconditioners for serial
   applications
   - The \p SparseMatrix, \p NumericVector, and \p LinearSolver
   allow for transparent switching between solver packages.  Adding
   a new solver interface is as simple as deriving from these classes

   - Mesh IO & Format Translation Utilities
   - Ideas Universal (UNV) format (.unv) with support through
   \p MeshData for arbitrary float data, like boundary conditions,
   associated with mesh entities
   - Sandia National Labs ExodusII format (.exd)
   - Amtec Engineering's Tecplot binary format (.plt)
   - Amtec Engineering's Tecplot ascii format (.dat)
   - Los Alamos National Labs GMV format (.gmv)
   - AVS Unstructured UCD format (.ucd)

   - Mesh Creation & Modification Utilities
   - refine or coarsen a mesh: prescribed, level-one-compatible, or uniform
   - build equispaced n-cubes out of \p Edge2, \p Tri3, \p Tri6,
   \p Quad4, \p Quad8, \p Quad9, \p Hex8, \p Hex20, \p Hex27
   - build circles/spheres out of \p Tri3, \p Tri6, \p Quad4,
   \p Quad8, \p Quad9, \p Hex8
   - add infinite elements to a volume-based mesh, handle symmetry planes
   - convert \p Quad4, \p Quad8, \p Quad9 to \p Tri3, \p Tri6
   - convert a mesh consisting of any of the fore-mentioned
   n-dimensional linear elements to their second-order
   counterparts
   - distort/translate/rotate/scale a mesh
   - determine bounding boxes/spheres
   - extract the mesh boundary for boundary condition handling or
   as a separate mesh
*/

// Local Variables:
// mode: html
// End:
