<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<!-- $Id: externalsoftware.php,v 1.5 2006-06-28 16:49:37 benkirk Exp $ -->

<html>
<head>
  <title>Packages used by libMesh</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("download",$root)?>

<div class="content">
<h1>External Packages</h1>
<code>libMesh</code> interfaces to a number of high-quality software packages to provide certain functionality.  This page provides a list of packages and a description of their use in the library.

<hr>
<h2>PETSc - Parallel Linear & Nonlinear Solvers</h2> <a href="http://www-unix.mcs.anl.gov/petsc/petsc-2">PETSc</a>
<h3>MPI</h3>
<h3>SLEPc</h3> <a href="http://www.grycap.upv.es/slepc">SLEPc</a>
<h3>BLAS</h3>
<h3>LAPACK</h3>

<hr>
<h2>LASPACK - Serial Linear Solvers</h2>



<hr>
<h2>Mesh Generation</h2>
<h3>Triangle</h3> <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle</a>
<h3>Tetgen</h3> <a href="http://tetgen.berlios.de">Tetgen</a>
<h3>EXODUS II</h3> EXODUS II, available via the <a href="http://endo.sandia.gov/SEACAS/Documentation/SEACAS.html">Sandia Engineering Analysis Code Access System</a>, is a model developed to store and retrieve data for finite element analyses. It is used for preprocessing (problem definition), postprocessing (results visualization), as well as code to code data transfer. An EXODUS II data file is a random access, machine independent, binary file that is written and read via C, C++, or Fortran library routines which comprise the Application Programming Interface. <code>libMesh</doc> can detect a valid EXODUS II installation and use it to read EXODUS II mesh files.  



<hr>
<h2>Visualization & Post-Processing</h2>

<h3>Tecplot</h3> <a href="http://www.tecplot.com">Tecplot</a> is a high-quality engineering and scientific visualization package.  <code>libMesh</code> can write simulation data in either ASCII or binary formatted Tecplot files.

<h3>GMV</h3> The <a href="http://laws.lanl.gov/XCM/gmv/GMVHome.html"> General Mesh Viewer</a> is "an easy to use, 3D scientific visualization tool designed to view simulation data from any type of structured or unstructured mesh."  GMV is developed at <a href="http://www.lanl.gov">Los Alamos National Laboratory</a> and is freely available for a wide range of platforms.  <code>libMesh</code> can write simulation data directly in the GMV file format.



<hr>
<h2>Utilities</h2>
<h3>XDR</h3>
<h3>GetPot</h3> <a href="http://getpot.sourceforge.net">GetPot</a>




<br>
<br>
<br>
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
