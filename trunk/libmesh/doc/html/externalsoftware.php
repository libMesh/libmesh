<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>Packages used by libMesh</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("externalsoftware",$root)?>

<div class="content">
<code>libMesh</code> interfaces to a number of high-quality software packages to provide certain functionality.  This page provides a list of packages and a description of their use in the library.

<ul>
 <li> <a href="http://www-unix.mcs.anl.gov/petsc/petsc-2">PETSc</a> </li>
 <li> <a href="http://www.grycap.upv.es/slepc">SLEPc</a> </li>
 <li> <a href="http://www.cs.cmu.edu/~quake/triangle.html">Triangle</a> </li>
 <li> <a href="http://tetgen.berlios.de">Tetgen</a> </li>
 <li> <a href="http://www.tecplot.com">Tecplot</a> </li>
 <li> <a href="http://getpot.sourceforge.net">GetPot</a> </li>
</ul>




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
