<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Examples</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("examples",$root)?>

<div class="content">

<h1>A Series of Example Programs</h1>
The following series of example programs have been
designed to get you started on the right foot.
For the most part they are arranged in order of
increasing complexity, and could be attempted in
that order.  Click the links below, or use the
menu on the left to navigate the examples.



<h2><a href="ex0.php">Example 0</a> - Solving 1D PDE Using Adaptive Mesh Refinement</h2>
This example demonstrates how to solve a simple 1D problem using
adaptive mesh refinement. The PDE that is solved is: -epsilon*u''(x) +
u(x) = 1, on the domain [0,1] with boundary conditions u(0) = u(1) = 0
and where epsilon << 1.

<br>
The approach used to solve 1D problems in libMesh is virtually
identical to solving 2D or 3D problems, so in this sense this example
represents a good starting point for new users. Note that many
concepts are used in this example which are explained more fully in
subsequent examples.




<h2><a href="ex1.php">Example 1</a> - Creation of a Mesh Object</h2>
This is the first example program.  It simply demonstrates
how to create a mesh object.  A mesh is read from file,
information is printed to the screen, and the mesh is then
written.




<h2><a href="ex2.php">Example 2</a> - Defining a Simple System</h2>
The second example program demonstrates how to
create an equation system for a simple scalar system.  This
example will also introduce some of the issues involved with using Petsc
in your application.




<h2><a href="ex3.php">Example 3</a> - Solving a 2D Poisson Problem</h2>
This is the third example program.  It builds on
the second example program by showing how to solve a simple
Poisson system.  This example also introduces the notion
of customized matrix assembly functions, working with an
exact solution, and using element iterators.
We will not comment on things that
were already explained in the second example.




<h2><a href="ex4.php">Example 4</a> - Solving a 2D or 3D Poisson Problem in Parallel</h2>
This is the fourth example program.  It builds on
the third example program by showing how to formulate
the code in a dimension-independent way.  Very minor
changes allow rhe example will allow the problem to be
solved in two or three dimensions and in parallel. 




<h2><a href="ex5.php">Example 5</a> - Run-time Selection of Quadrature Rules</h2>
This example changes the previous example by enabling
run-time selection of quadrature rules. 




<h2><a href="ex6.php">Example 6</a> - Infinite Elements for the Wave Equation</h2>
This example introduces "infinite elements" which may be used for
certain classes of applications.  The wave equation is solved in this
example. <i>For this example to work you must have configured the
library with the --enable-ifem option</i>




<h2><a href="ex7.php">Example 7</a> - Introduction to Complex Numbers and the "FrequencySystem"</h2>
This is the seventh example program.  It builds on
the previous example programs, introduces complex
numbers and the FrequencySystem class to solve a 
simple Helmholtz equation grad(p)*grad(p)+(omega/c)^2*p=0,
for multiple frequencies rather efficiently.
The FrequencySystem class offers two solution styles,
namely to solve large systems, or to solve
moderately-sized systems fast, for multiple frequencies.
The latter approach is implemented here.
For this example the library has to be compiled with
complex numbers enabled. 




<h2><a href="ex8.php">Example 8</a> - The Newmark System and the Wave Equation</h2>
This example solves the wave equation in a hybrid-mesh pipe.  The mesh
consists of <code>HEX8</code> and <code>PRISM6</code> element types.
The pressure at a point in the pipe is extracted and can be plotted as
a function of time.




<h2><a href="ex9.php">Example 9</a> - Solving a Transient Linear System in Parallel</h2>
This example shows how a simple, linear transient
system can be solved in parallel.  The system is simple
scalar convection-diffusion with a specified external
velocity.  The initial condition is given, and the
solution is advanced in time with a standard Crank-Nicholson
time-stepping strategy.




<h2><a href="ex10.php">Example 10</a> - Solving a Transient System with Adaptive Mesh Refinement</h2>
This example shows how a simple, linear transient
system can be solved in parallel.  The system is simple
scalar convection-diffusion with a specified external
velocity.  The initial condition is given, and the
solution is advanced in time with a standard Crank-Nicholson
time-stepping strategy.  This example differs from the previous
example by employing adaptive mesh refinement (AMR) and the
Kelly et. al. error indicator.




<h2><a href="ex11.php">Example 11</a> - Solving a System of Equations</h2>
This example shows how to solve a simple, linear system of equations.  The
familiar Stokes equations for incompressible fluid flow are solved.  To satisfy
the LBB criterion different approximation spaces are used for the velocity and
pressure fields.




<h2><a href="ex12.php">Example 12</a> - Using the <code>MeshData</code> class</h2>
This example describes the use of the <code>MeshData</code> class.
More on this later.



<h2><a href="ex13.php">Example 13</a> - Unsteady Navier-Stokes Equations - Unsteady Nonlinear System</h2>
This example shows how a simple, unsteady, nonlinear system of equations
can be solved in parallel.  The system of equations are the familiar
Navier-Stokes equations for low-speed incompressible fluid flow.  This
example introduces the concept of the inner nonlinear loop for each
timestep, and requires a good deal of linear algebra number-crunching
at each step.  If you have the General Mesh Viewer (GMV) installed,
the script movie.sh in this directory will also take appropriate screen
shots of each of the solution files in the time sequence.  These rgb files
can then be animated with the "animate" utility of ImageMagick if it is
installed on your system.  On a PIII 1GHz machine in debug mode, this
example takes a little over a minute to run.  If you would like to see
a more detailed time history, or compute more timesteps, that is certainly
possible by changing the n_timesteps and dt variables below.

<h2><a href="ex14.php">Example 14</a> - Laplace's Equation in an L-Shaped Domain</h2>
This example solves the Laplace equation on the classic "L-shaped"
domain with adaptive mesh refinement.  In this case, the exact
solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta), but
the standard Kelley error indicator is used to estimate the error.
The initial mesh contains three QUAD9 elements which represent the
standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
i.e.
Element 0: [-1,0]x[ 0,1]
Element 1: [ 0,1]x[ 0,1]
Element 2: [-1,0]x[-1,0]
The mesh is provided in the standard libMesh ASCII format file
named "lshaped.xda".  In addition, an input file named "ex14.in"
is provided which allows the user to set several parameters for
the solution so that the problem can be re-run without a
re-compile.  The solution technique employed is to have a
refinement loop with a linear solve inside followed by a
refinement of the grid and projection of the solution to the new grid
In the final loop iteration, there is no additional
refinement after the solve.  In the input file "ex14.in", the variable
"max_r_steps" controls the number of refinement steps,
"max_r_level" controls the maximum element refinement level, and
"refine_percentage" / "coarsen_percentage" determine the number of
elements which will be refined / coarsened at each step.

<h2><a href="ex15.php">Example 15</a> - Biharmonic Equation</h2>
This example solves the Biharmonic equation on a square domain using a
Galerkin formulation with C1 elements approximating the H^2_0 function
space.  The initial mesh contains two TRI6 elements.  The mesh is
provided in the standard libMesh ASCII format file named "domain.xda".
In addition, an input file named "ex15.in" is provided which allows
the user to set several parameters for the solution so that the
problem can be re-run without a re-compile.  The solution technique
employed is to have a refinement loop with a linear solve inside
followed by a refinement of the grid and projection of the solution to
the new grid In the final loop iteration, there is no additional
refinement after the solve.  In the input file "ex15.in", the variable
"max_r_steps" controls the number of refinement steps, and
"max_r_level" controls the maximum element refinement level.

<h2><a href="ex16.php">Example 16</a> - Solving an Eigen Problem</h2>
This example introduces the EigenSystem and shows how libMesh can be
used for eigenvalue analysis.  For solving eigen problems, libMesh
interfaces SLEPc
(<a href="http://www.grycap.upv.es/slepc/">www.grycap.upv.es/slepc/</a>)
which again is based on PETSc.  Hence, this example will only work if
the library is compiled with SLEPc support enabled.  In this example
some eigenvalues for a standard symmetric eigenvalue problem
A*x=lambda*x are computed, where the matrix A is assembled according
to a mass matrix.

<h2><a href="ex17.php">Example 17></a> - Solving a generalized Eigen Problem</h2>
This example shows how the previous EigenSolver example
can be adapted to solve generailzed eigenvalue problems.
For solving eigen problems, libMesh interfaces
SLEPc
(<a href="http://www.grycap.upv.es/slepc/">www.grycap.upv.es/slepc/</a>)
which again is based on PETSc.
Hence, this example will only work if the library is compiled
with SLEPc support enabled.

In this example some eigenvalues for a generalized symmetric
eigenvalue problem A*x=lambda*B*x are computed, where the
matrices A and B are assembled according to stiffness and
mass matrix, respectively.

<h2><a href="ex18.php">Example 18</a> - Unsteady Navier-Stokes Equations with DiffSystem</h2>
This example shows how the transient nonlinear problem from
example 13 can be solved using the DiffSystem class framework to
simplify the user-implemented equations.

</div>

<!--
<div id="navBeta">
</div>
-->

<?php make_footer() ?>

</body>
</html>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
