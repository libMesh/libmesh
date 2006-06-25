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
// we have split up into several sections.  Unfortunately, it doesnt
// work for internet explorer!
if (NW_IS_IE)
{
  echo "<img class=\"title_image\" src=\"images/libmesh_mesh_small.png\">";
}

else
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
?>
The <code>libMesh</code> library is a C++ framework for the numerical
simulation of partial differential equations on serial and parallel
platforms. Development began in March 2002 with the intent of
providing a friendly interface to a number of high-quality software
packages that are currently available. Currently the library supports
2D and 3D steady and transient finite element simulations.
<a href="http://www-unix.mcs.anl.gov/petsc/petsc-2">PETSc</a> is
currently used for the solution of linear systems on both serial and
parallel platforms, and
<a href="http://www.tu-dresden.de/mwism/skalicky/laspack/laspack.html">LASPack</a>
is included with the library to provide linear solver support on serial machines.

<br>
<br>
The <code>libMesh</code> library is actively developed at The
University of Texas at Austin in the <a href="http://www.cfdlab.ae.utexas.edu">CFDLab</a> and at
Technische Universit&auml;t Hamburg-Harburg,
<a href="http://www.mub.tu-harburg.de/index_e.html">Modelling and Computation</a> in Germany.
Many thanks to <a href="http://sourceforge.net">SourceForge</a> for
<a href="http://sourceforge.net/projects/libmesh">hosting the
project</a>.  You can find out what is currently happening in the
development branch by checking out the
<a href="http://libmesh.cvs.sourceforge.net/libmesh/libmesh">CVS
Logs</a> online.

<br>
<br>
A major goal of the library is to provide support for adaptive mesh
refinement (AMR) computations in parallel while allowing a research
scientist to focus on the physics they are modeling.

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


