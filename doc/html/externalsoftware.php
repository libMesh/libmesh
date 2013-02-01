<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<!-- $Id: externalsoftware.php,v 1.7 2007-09-09 20:02:19 benkirk Exp $ -->

<html>
<head>
  <title>Packages used by libMesh</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("externalsoftware",$root)?>

<div class="content">
<h1>External Packages</h1>
<code>libMesh</code> interfaces to a number of high-quality software packages to provide certain functionality.  This page provides a list of packages and a description of their use in the library.

<h2>MPI</h2> The <a href="http://www-unix.mcs.anl.gov/mpi">Message Passing Interface</a> is a standard for parallel programming using the message passing model. PETSc requires MPI for its functionality.  <code>libMesh</code> makes use of MPI to when running in parallel for certain operations, including its <code>ParallelMesh</code> distributed-memory, fully unstructured mesh implementation. 

<h2>TBB</h2> Since February 2008 <code>libMesh</code> can be configured to use the <a href="http://threadingbuildingblocks.org/">Threading Building Blocks</a> for thread-based parallelism on shared memory machines.  Several key algorithms in the library have been refactored to be multithreaded, and this effort will continue as additional profiling reveals additional serial bottlenecks.  It is envisioned that for certain classes of problems multilevel parallelism (e.g. message passing between nodes and threading within nodes) will prove more scalable than message passing alone, especially with the introduction of commodity multi-core processors.  The reality is that for implicit problems this can only be achieved with a parallel linear algebra library that also uses multilevel parallelism.

<hr>
<h2>PETSc - Parallel Linear & Nonlinear Solvers</h2> The Portable, Extensible Toolkit for Scientific Computation <a href="http://www.mcs.anl.gov/petsc">(PETSc)</a> is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.

<h3>SLEPc</h3> The Scalable Library for Eigenvalue Computations <a href="http://www.grycap.upv.es/slepc">(SLEPc)</a> is a library for the solution of large scale sparse eigenvalue problems on parallel computers. It is an extension of PETSc and can be used for either standard or generalized eigenproblems, with real or complex arithmetic.

<h3>BLAS</h3> The <a href="http://www.netlib.org/blas">Basic Linear Algebra Subprograms</a> are routines that provide standard building blocks for performing basic vector and matrix operations. The Level 1 BLAS perform scalar, vector and vector-vector operations, the Level 2 BLAS perform matrix-vector operations, and the Level 3 BLAS perform matrix-matrix operations. High-performance implementations of the BLAS are generally provided by computer hardware manufacturers for a particular architecture.  PETSc makes extensive use of the BLAS hence a high-performance BLAS implementation is key to achieving high performance from the PETSc linear solvers.

<h3>LAPACK</h3>
<a href="http://www.netlib.org/lapack">LAPACK</a> s written in Fortran77 and provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems. The associated matrix factorizations (LU, Cholesky, QR, SVD, Schur, generalized Schur) are also provided, as are related computations such as reordering of the Schur factorizations and estimating condition numbers. Dense and banded matrices are handled, but not general sparse matrices.  PETSc makes use of LAPACK in several computational kernels within its linear solver framework.

<hr>
<h2>Trilinos</h2>
<a href="http://trilinos.sandia.gov">The Trilinos Project</a> is an effort to develop and implement robust algorithms and enabling technologies using modern object-oriented software design, while still leveraging the value of established libraries. It emphasizes abstract interfaces for maximum flexibility of component interchanging, and provides a full-featured set of concrete classes that implement all abstract interfaces. Current development efforts will permit <code>libMesh</code> can use the Epetra parallel matrix and vector data structures, as well as the AztecOO parallel linear solvers.

<hr>
<h2>LASPack - Serial Linear Solvers</h2>
 <a href="http://www.tu-freiberg.de/urz/anwendungen/sprodukte/soft/LASPACK/laspack.html">LASPack</a> is an object-oriented package of iterative methods, multigrid solvers, and auxiliary routines for the iterative solution of linear systems. It does not run in parallel. There are data structures for vectors, general and square matrices, and preconditioners; a large number of accessing and manipulating these objects is available.

<!--
Main features:
<ul>
  <li>The primary aim of LASPack is the implementation of efficient iterative methods for the solution of systems of linear equations. All routines and data structures are optimized for effective usage of resources especially with regard to large sparse matrices. The package can be accessed from an application through a straightforward interface defined in the form of procedure calls.</li>
  <li>Besides the obligatory Jacobi, successive over-relaxation, Chebyshev, and conjugate gradient solvers, LASPack contains selected state-of-the-art algorithms which are commonly used for large sparse systems:
  <ul>
   <li>CG-like methods for non-symmetric systems: CGN, GMRES, BiCG, QMR, CGS, and BiCGStab,</li>
   <li>multilevel methods such as the multigrid and conjugate gradient methods preconditioned by multigrid and BPX preconditioners</li>
  </ul>
  All the above solvers are applicable not only to the positive definite or non-symmetric matrices, but are also adopted for singular systems (e.g. arising from discretization of Neumann boundary value problems).</li>
  <li>The implementation is based on an object-oriented approach (although it is programmed in C). Vectors and matrices are defined as new data types in connection with the corresponding supporting routines. The basic operations are implemented in such a way that they allow the programming of linear algebra algorithms in a natural way.</li>
  <li>LASPack is extensible in a simple manner. An access to the internal representation of vectors and matrices is not necessary and is avoided, as required of the object-oriented programming. This allows an improvement of algorithms or a modification of data structures with no adjustment of application programs using the package</li>
</ul>
-->


<hr>
<h2>Mesh Generation</h2>
<h3>Triangle</h3> <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle</a> is the definitive two-dimensional delaunay triangulator written by Jonathan Richard Shewchuk. <code>libMesh</code> can use Triangle to produce Delaunay triangulizations for hybrid-element (not just triangles) input meshes. 

<h3>Tetgen</h3> <a href="http://tetgen.berlios.de"> Tetgen</a>generates the Delaunay tetrahedralization, Voronoi diagram, constrained Delaunay tetrahedralizations and quality tetrahedral meshes. The main goal of TetGen is to generate suitable meshes for solving partial differential equations by finite element or finite volume methods.


<h3>EXODUS II</h3> EXODUS II, available via <a href="http://sourceforge.net/projects/exodusii">sourceforge</a>, is a model developed to store and retrieve data for finite element analyses. It is used for preprocessing (problem definition), postprocessing (results visualization), as well as code to code data transfer. An EXODUS II data file is a random access, machine independent, binary file that is written and read via C, C++, or Fortran library routines which comprise the Application Programming Interface. <code>libMesh</code> contains source code for the EXODUS II library in <code>./contrib/exodusii</code> and can use it to read EXODUS II mesh files.  Additional information may be found on the <a href="http://endo.sandia.gov/SEACAS/Documentation/SEACAS.html">Sandia Engineering Analysis Code Access System</a>



<hr>
<h2>Visualization & Post-Processing</h2>

<h3>Tecplot</h3> <a href="http://www.tecplot.com">Tecplot</a> is a high-quality engineering and scientific visualization package.  <code>libMesh</code> can write simulation data in either ASCII or binary formatted Tecplot files.

<h3>GMV</h3> The <a href="http://laws.lanl.gov/XCM/gmv/GMVHome.html"> General Mesh Viewer</a> is "an easy to use, 3D scientific visualization tool designed to view simulation data from any type of structured or unstructured mesh."  GMV is developed at <a href="http://www.lanl.gov">Los Alamos National Laboratory</a> and is freely available for a wide range of platforms.  <code>libMesh</code> can write simulation data directly in the GMV file format.



<hr>
<h2>Utilities</h2>

<h3>XDR</h3> The <a href="http://www.faqs.org/rfcs/rfc1014.html">XDR: External Data Representation Standard</a> is a standard for the description and encoding of data.  It is useful for transferring data between different computer architectures, and as such provides a very simple, portable approach for writing platform independent binary files. <code>libMesh</code> uses the <a href="doxygen/classXdr.php">Xdr</a> class to provide a uniform interface to input/output operations for files in either XDR binary or ASCII text formats.

<h3>GetPot</h3> <a href="http://getpot.sourceforge.net">GetPot</a> is a powerful command line and configuration file parsing for C++, Python, Ruby and Java. This tool provides many features, such as separate treatment for options, variables, and flags, unrecognized object detection, prefixes and much more. <code>libMesh</code> uses GetPot to parse command line options upon initialization and for input file parsing in some of the <a href="examples.php">examples</a>.

<h3>libHilbert</h3> <a href="http://flame.cs.dal.ca/~chamilto/hilbert">libHilbert</a> is library written by <a href="http://flame.cs.dal.ca/~chamilto">Chris Hamilton</a> for producing compact Hilbert indices for multidimensional data.  <code>libMesh</code> uses <code>libHilbert</code> to assign unique, global identifiers for nodes an elements which are independent of a particular domain decomposition.



<br>
<br>
<br>


<?php make_footer() ?>

</body>
</html>

<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
