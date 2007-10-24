<!-- $Id: index.php,v 1.21 2007-09-09 19:42:07 benkirk Exp $ -->
<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<HTML>
<head>
  <title>libMesh - C++ Finite Element Library</title>
  <META name="description" content="A parallel adaptive C++ finite element library for simulating partial differential equations.">
  <META name="keywords" content="finite element, finite element method, finite element modeling, finite element modelling, finite element software, finite element methods, science math numerical analysis, science math numerical analysis software, numerical analysis research group, partial differential equations, applied partial differential equations, elliptic parabolic hyperbolic partial differential equations, finite element analysis modeling software, system of linear equation solver, simultaneous linear equation solver, computational fluid dynamics software">
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
if (NW_IS_GECKO || NW_IS_MAC)
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

//  $pic_widths  = array( 56, 122, 180, 245, 310, 373, 410 );
//  $pic_heights = array( 13,  13,  12,  12,  12,  13, 128 );
//
//  for ($i=0; $i<7; $i++)
//    {
//      echo "<img class=\"slant_right\" width=$pic_widths[$i] height=$pic_heights[$i] src=\"images/libmesh_mesh_layer_rhs_0$i.png\">";
//    }
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
The <code>libMesh</code> library provides a framework for the
numerical simulation of partial differential equations using arbitrary
unstructured discretizations on serial and parallel platforms. A major
goal of the library is to provide support for adaptive mesh refinement
(AMR) computations in parallel while allowing a research scientist to
focus on the physics they are modeling.

<!-- Development began in March 2002 with the intent of providing a -->
<!-- friendly interface to a number of high-quality software packages that -->
<!-- are publicly available. -->

<br>
<br>
<code>libMesh</code> currently supports 1D, 2D, and 3D steady and transient finite element simulations.
The library makes use of high-quality, existing software whenever possible.
<a href="http://www-unix.mcs.anl.gov/petsc/petsc-2">PETSc</a> is
used for the solution of linear systems on both serial and
parallel platforms, and
<a href="http://www.mgnet.org/mgnet/Codes/laspack/html/laspack.html">LASPack</a>
is included with the library to provide linear solver support on serial machines.
An optional interface to <a href="http://www.grycap.upv.es/slepc">SLEPc</a> is also
provided for solving both standard and generalized eigenvalue problems.
A complete list of external applications used in the library may be found <a href="externalsoftware.php">here</a>.


<br>
<br>
The <code>libMesh</code> library has been actively developed at The
University of Texas at Austin in the 
<a href="http://www.cfdlab.ae.utexas.edu">CFDLab</a> since March 2002.  Major 
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
<a href="http://libmesh.svn.sourceforge.net/viewvc/libmesh/trunk/libmesh">SVN
Logs</a> online, and you can see many people are downloading the library
on the <a href="http://sourceforge.net/project/stats/?group_id=71130&ugn=libmesh">statistics</a> page.

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


