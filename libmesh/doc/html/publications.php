<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Publications</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("publications",$root)?>

<div class="content">
<h1>Publications</h1>

<ul>

<?php echo
"<li> John's <a href=\"http://www.cfdlab.ae.utexas.edu/~peterson/masters.pdf\">
Masters Report</a>. </li>";
?>

<?php echo
"<li> Ben's <a href=\"http://www.cfdlab.ae.utexas.edu/~peterson/ben_dissertation.pdf\">
PhD Dissertation</a>. </li>";
?>

<?php echo
"<li>  A general <a href=\"howto/howto.pdf\">HOWTO</a>
document by M. Luthi containing some hints
and programming tips for writing effective LibMesh programs. </li>";
?>

<?php echo
"<li>  A <a href=\"xda_format/xda_format.pdf\">description</a>
of the XDA file format used by LibMesh. </li>";
?>

<?php echo
"<li> Texas Advanced Computing Center
<a href=\"http://www.tacc.utexas.edu/general/press/announcements/20040112_01.php\">
press release</a> commemorating the launch of the Lonestar cluster. </li>";
?>

<?php echo
"<li>A <a href=\"http://ondrej.certik.cz/libmesh/fem.ps\">description</a> of the Newmark
System class by Ondrej Certik.</li>";
?>

<?php echo
"<li>libMesh: A C++ Library for Parallel Adaptive
 Mesh Refinement/Coarsening Simulations - Engineering with Computers.
      <ul>
	<li>
	  <a href=\"http://cfdlab.ae.utexas.edu/~benkirk/libmesh-ewc-preprint.pdf\">Preprint</a>
	</li>
	<li>
	  <a href=\"http://cfdlab.ae.utexas.edu/~peterson/ewc_journal_version.pdf\">Journal Version</a>
	</li>
      </ul>
</li>";
?>


</ul>
</div>

<br>
<br>
<br>
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
