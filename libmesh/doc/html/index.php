<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<HTML>
<head>
  <title>libMesh -- A C++ Finite Element Library</title>
  <?php load_style($root); ?>
</head>

<BODY>

<?php make_navigation("normal",$root)?>


<div class="content">
<!-- <h1>libMesh</h1> -->
   <img src="images/libmesh_mesh_small.png"><br>

 The libMesh library is a C++ framework for the numerical simulation of
partial differential equations on serial and parallel
platforms. Development began in March 2002 with the intent of
providing a friendly interface to a number of high-quality software
packages that are currently available. Currently the library supports
2D and 3D steady and transient finite element simulations.  PETSc is
currently used for the solution of linear systems on both serial and
parallel platforms, and LASPACK is included with the library to
provide linear solver support on serial machines.

<br>
<br>
The libMesh library is actively developed at The
University of Texas at Austin in the <a href="http://www.cfdlab.ae.utexas.edu">CFDLab</a> and at
Technische Universit&auml;t Hamburg-Harburg, <a
href="http://www.mum.tu-harburg.de/english"> Mechanics and Ocean
Engineering</a> in Germany. Many thanks to <a
href="http://sourceforge.net"> SourceForge </a> for <a
href="http://sourceforge.net/projects/libmesh"> hosting the
project</a>.  You can find out what is currently happening in the
development branch by checking out the <a
href="http://cvs.sourceforge.net/cgi-bin/viewcvs.cgi/libmesh"> CVS
Logs</a> online.

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


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>