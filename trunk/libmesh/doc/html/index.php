<!-- $Id: index.php,v 1.20 2007-07-01 00:43:38 benkirk Exp $ -->
<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<HTML>
<head>
  <title>libMesh - A C++ Finite Element Library</title>
  <META name="description" content="A parallel adaptive C++ framework for the simulation of partial differential equations.">
  <META name="keywords" content="science math numerical analysis, science math numerical analysis software, numerical analysis research group, partial differential equations, applied partial differential equations, elliptic parabolic hyperbolic partial differential equations, finite element analysis modeling software, system of linear equation solver, simultaneous linear equation solver, computational fluid dynamics software">
  <META name="robots" content="index, follow">
  <?php load_style($root); ?>
</head>

<BODY>

<?php make_navigation("normal",$root)?>


<div class="content">

<!-- <img class="title_image" src="images/libmesh_mesh_small.png"> -->
<!-- <br><br><br> -->
<?php

// This php script makes the words wrap nicely around the image which
// we have split up into several sections.  This apparently only works
// for Gecko-enabled browsers, e.g. Firefox.
if (NW_IS_GECKO)
{
  // You can get these numbers using the 'file' command on most linux systems.
  $pic_widths  = array( 417, 417, 417, 417, 417, 400, 345, 270, 205, 145, 80 );
  $pic_heights = array( 26,   26,  26,  26,  26,  13,  13,  13,  13,  13, 13 );

  for ($i=0; $i<11; $i++)
    {
      if ($i < 10)
	{
	  echo "<img class=\"slant\" width=$pic_widths[$i] height=$pic_heights[$i] src=\"images/libmesh_mesh_layer_0$i.png\">";
	}
      
      else
	{
	  echo "<img class=\"slant\" width=$pic_widths[$i] height=$pic_heights[$i] src=\"images/libmesh_mesh_layer_10.png\">";
	}
    }
}

// Do not do slanting text
else
{
  // Use a slightly different CSS module for IE.
  if (NW_IS_IE > 5)
    {
  echo "<img class=\"title_image_ie\" src=\"images/libmesh_mesh_small.png\">";
    }
  else
    {
  echo "<img class=\"title_image\" src=\"images/libmesh_mesh_small.png\">";
    }
}
?>
The <code>libMesh</code> library is a C++ framework for the numerical
simulation of partial differential equations on serial and parallel
platforms. Development began in March 2002 with the intent of
providing a friendly interface to a number of high-quality software
packages that are publicly available. Currently the library supports
1D, 2D, and 3D steady and transient finite element and finite volume simulations.
<a href="http://www-unix.mcs.anl.gov/petsc/petsc-2">PETSc</a> is
currently used for the solution of linear systems on both serial and
parallel platforms, and
<a href="http://www.tu-dresden.de/mwism/skalicky/laspack/laspack.html">LASPack</a>
is included with the library to provide linear solver support on serial machines.
An optional interface to <a href="http://www.grycap.upv.es/slepc">SLEPc</a> is also
provided for solving both standard and generalized eigenvalue problems.


<br>
<br>
The <code>libMesh</code> library is primarily developed at The
University of Texas at Austin in the 
<a href="http://www.cfdlab.ae.utexas.edu">CFDLab</a>.  Major 
contributions have come from developers at the Technische 
Universit&auml;t Hamburg-Harburg
<a href="http://www.mub.tu-harburg.de/index_e.html">Institute of Modelling and Computation</a>, 
and recent contributions have been made by CFDLab graduates at
<a href="http://www.sandia.gov/">Sandia National Laboratories</a> and <a href="http://www.nasa.gov/">NASA</a> <a href="http://www.nasa.gov/centers/johnson/home/index.html">Lyndon B. Johnson Space Center</a>.
The <code>libMesh</code> <a href="http://libmesh.sf.net/developers.php">developers</a> welcome contributions
in the form of patches and bug reports (preferably with a minimal test case that reliably reproduces the error)
to the official <a href="http://sourceforge.net/mail/?group_id=71130">mailing lists</a>.
Many thanks to <a href="http://sourceforge.net">SourceForge</a> for
<a href="http://sourceforge.net/projects/libmesh">hosting the
project</a>.  You can find out what is currently happening in the
development branch by checking out the
<a href="http://libmesh.cvs.sourceforge.net/libmesh/libmesh">CVS
Logs</a> online, and you can see many people are downloading the library
on the <a href="http://sourceforge.net/project/stats/?group_id=71130&ugn=libmesh">statistics</a> page.

<br>
<br>
A major goal of the library is to provide support for adaptive mesh
refinement (AMR) computations in parallel while allowing a research
scientist to focus on the physics they are modeling.
The library makes use of high-quality, existing software whenever possible.
A complete list of external applications used in the library may be found <a href="externalsoftware.php">here</a>.

<br>
<br>

</div>
<br>

<!--
<?php if (1==1) {echo "PHP IS WORKING";} else { ?>
PHP Does not Work
<?php } ?>
-->

<?php make_footer() ?>

</BODY>
</HTML>


