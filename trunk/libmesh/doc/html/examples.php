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

<h2><a href="ex2.php">Example 2</a> -- Using PETSc to Solve a Simple System</h2>
The second example program demonstrates how to
create an equation system for a simple scalar system.  This
example will also introduce some of the issues involved with using Petsc
in your application.

<br>
<br>

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