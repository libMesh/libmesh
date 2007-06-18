<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>LibMesh Developers</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("developers",$root)?>

<?php function dev_table_entry($pic, $name, $homepage, $titles, $institutions) { ?>
<table>
<tr>
  <td>
    <?php echo "<img src=\"$pic\">"; ?>
  </td>
  <td>
  <?php echo "<b>Name</b>: <a href=\"$homepage\">$name</a> <br>"; ?>
  <?php echo "<b>Title(s)</b>: $titles <br>"; ?>
  <?php echo "<b>Institution(s)</b>: $institutions"; ?>
  </td>
</tr>
</table>
<?php } ?>

<div class="content">
<h1>Active LibMesh Developers</h1>

<?php
dev_table_entry("images/benkirk.jpg",
                "Dr. Benjamin S. Kirk",
                "http://www.cfdlab.ae.utexas.edu/~benkirk",
                "libMesh Project Manager",
	        "<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<?php
dev_table_entry("images/spetersen.jpg",
		"Steffen Petersen",
		"http://www.mub.tu-harburg.de/deutsch/mitarbeiter/petersen.html",
		"PhD Student, libMesh Developer",
		"<a href=\"http://www.tu-harburg.de\">Hamburg University of Technology</a>");
?>

<?php
dev_table_entry("images/jwpeterson.jpg",
                "John W. Peterson",
                "http://www.cfdlab.ae.utexas.edu/~peterson",
		"PhD Student, libMesh Developer",
		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<?php
dev_table_entry("images/roystgnr.jpg",
                "Roy Stogner",
                "http://www.cfdlab.ae.utexas.edu/~roystgnr",
		"PhD Student, libMesh Developer",
		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<?php
 dev_table_entry("images/drgasto_s.jpg",
 		"Derek Gaston",
 		"http://www.cfdlab.ae.utexas.edu/~drgasto",
 		"libMesh Developer",
 		"<a href=\"http://www.sandia.gov\">UT-Austin, Sandia Nat'l Lab</a>");
?>

<p>

<h1>Past LibMesh Developers</h1>

<?php
 dev_table_entry("images/mikeando.jpg",
                 "Michael L. Anderson",
                 "http://www.cfdlab.ae.utexas.edu/~michaela",
 		"Post-doctoral Fellow, libMesh Developer",
 		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>,
 		<br><a href=\"http://www.uwa.edu.au\">University of Western Australia</a>");
?>


<?php
 dev_table_entry("images/bbarth2.jpg",
                 "William L. Barth",
                 "http://www.cfdlab.ae.utexas.edu/~bbarth",
 		"PhD Student, libMesh Developer",
 		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>


<?php
 dev_table_entry("images/blank.png",
 		"Daniel Dreyer",
 		"http://www.audi.com",
 		"PhD Student, libMesh Developer",
 		"<a href=\"http://www.audi.com\">Audi AG</a>,
// 		<br><a href=\"http://www.tu-harburg.de\">Technical University Hamburg-Harburg</a>");
?>

<?php
 dev_table_entry("images/knezed01.jpg",
 		"David Knezevic",
 		"http://web.comlab.ox.ac.uk/oucl/work/david.knezevic/",
 		"PhD Student, libMesh Developer",
 		"<a href=\"http://web.comlab.ox.ac.uk/oucl/work/david.knezevic/\">Oxford University Computing Laboratory</a>");
?>

<?php
 dev_table_entry("images/blank.png",
 		"Hendrik van der Heijden",
 		"http://www.tu-harburg.de",
 		"Student, libMesh Developer",
 		"<a href=\"http://www.tu-harburg.de\">Technical University Hamburg-Harburg</a>");
?>

<?php
 dev_table_entry("images/blank.png",
		"Florian Prill",
 		"http://www.tu-harburg.de",
 		"Student, libMesh Developer",
 		"<a href=\"http://www.tu-harburg.de\">Technical University Hamburg-Harburg</a>");
?>


</div>


<?php make_footer() ?>

</body>
</html>

