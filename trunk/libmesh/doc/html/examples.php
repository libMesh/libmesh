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
<br>
<br>



<h2><a href="ex1.php">Example 1</a> -- Creation of a Mesh Object</h2>
This is the first example program.  It simply demonstrates
how to create a mesh object.  A mesh is read from file,
information is printed to the screen, and the mesh is then
written.

<br>
<br>



<h2><a href="ex2.php">Example 2</a> -- Defining a Simple System</h2>
The second example program demonstrates how to
create an equation system for a simple scalar system.  This
example will also introduce some of the issues involved with using Petsc
in your application.

<br>
<br>



<h2><a href="ex3.php">Example 3</a> -- Solving a 2D Poisson Problem</h2>
This is the third example program.  It builds on
the second example program by showing how to solve a simple
Poisson system.  This example also introduces the notion
of customized matrix assembly functions, working with an
exact solution, and using element iterators.
We will not comment on things that
were already explained in the second example.

<br>
<br>



<h2><a href="ex4.php">Example 4</a> -- Solving a 2D or 3D Poisson Problem in Parallel</h2>
This is the fourth example program.  It builds on
the third example program by showing how to formulate
the code in a dimension-independent way.  Very minor
changes allow rhe example will allow the problem to be
solved in two or three dimensions and in parallel. 

<br>
<br>



<h2><a href="ex5.php">Example 5</a> -- Run-time Selection of Quadrature Rules</h2>
This example changes the previous example by enabling
run-time selection of quadrature rules. 

<br>
<br>



<h2><a href="ex6.php">Example 6</a> -- Infinite Elements for the Wave Equation</h2>
This example introduces "infinite elements" which may be used for certain classes
of applications.  The wave equation is solved in this example. <i>For this example to
work you must have configured the library with the --enable-ifem option</i> 

<br>
<br>

<h2><a href="ex10.php">Example 10</a> -- Solving a Transient System with Adaptive Mesh Refinement</h2>
This example shows how a simple, linear transient
system can be solved in parallel.  The system is simple
scalar convection-diffusion with a specified external
velocity.  The initial condition is given, and the
solution is advanced in time with a standard Crank-Nicholson
time-stepping strategy.
 


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