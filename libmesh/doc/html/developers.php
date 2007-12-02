<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Developers</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("developers",$root)?>

<?php function dev_table_entry($pic, $name, $homepage, $other, $titles, $institutions) { ?>
<table>
<tr>
  <td>
    <?php echo "<img src=\"$pic\">"; ?>
  </td>
  <td>
  <?php echo "<b>Name</b>: <a href=\"$homepage\">$name</a> $other <br>"; ?>
  <?php echo "<b>Title(s)</b>: $titles <br>"; ?>
  <?php echo "<b>Institution(s)</b>: $institutions"; ?>
  </td>
</tr>
</table>
<?php } ?>

<div class="content">
<h1>Active libMesh Developers</h1>

<?php
dev_table_entry("images/benkirk.jpg",
                "Dr. Benjamin S. Kirk",
                "http://www.cfdlab.ae.utexas.edu/~benkirk",
		"<a href=\"http://www.cfdlab.ae.utexas.edu/~benkirk/vitae.pdf\">(curriculum vitae)</a>",
                "Project Manager",
	        "<a href=\"http://www.nasa.gov/centers/johnson/home/index.html\">NASA Lyndon B. Johnson Space Center</a>");
?>

<?php
dev_table_entry("images/spetersen.jpg",
		"Steffen Petersen",
		"http://www.mub.tu-harburg.de/deutsch/mitarbeiter/petersen.html",
		"",
		"PhD Student, Developer",
		"<a href=\"http://www.tu-harburg.de\">Hamburg University of Technology</a>");
?>

<?php
dev_table_entry("images/jwpeterson.jpg",
                "John W. Peterson",
                "http://www.cfdlab.ae.utexas.edu/~peterson",
                "<a href=\"http://www.cfdlab.ae.utexas.edu/~peterson/resume.pdf\">(curriculum vitae)</a>",
		"PhD Student, Developer",
		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<?php
dev_table_entry("images/roystgnr.jpg",
                "Roy Stogner",
                "http://www.cfdlab.ae.utexas.edu/~roystgnr",
		"",
		"PhD Student, Chief Software Architect",
		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<?php
 dev_table_entry("images/drgasto_s.jpg",
 		"Derek Gaston",
 		"http://www.cfdlab.ae.utexas.edu/~drgasto",
		"",
 		"Developer",
 		"<a href=\"http://www.sandia.gov\">UT-Austin, Sandia Nat'l Lab</a>");
?>

<br>
<br>

<h1>Past libMesh Developers</h1>

<?php
 dev_table_entry("images/mikeando.jpg",
                 "Dr. Michael L. Anderson",
                 "http://www.cfdlab.ae.utexas.edu/~michaela",
		 "",
		 "Post-doctoral Fellow, Developer",
		 "<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>,
 		 <br><a href=\"http://www.uwa.edu.au\">University of Western Australia</a>");
?>


<?php
 dev_table_entry("images/bbarth2.jpg",
                 "Dr. William L. Barth",
                 "http://www.tacc.utexas.edu/general/staff/barth",
		 "",
		 "Consultant, Developer",
		 "<a href=\"http://www.tacc.utexas.edu\">Texas Advanced Computing Center</a>");
?>


<?php
 dev_table_entry("images/blank.png",
		 "Dr. Daniel Dreyer",
		 "http://www.audi.com",
		 "",
		 "Developer",
		 "<a href=\"http://www.audi.com\">Audi AG</a>");
?>

<?php
 dev_table_entry("images/knezed01.jpg",
		 "David Knezevic",
		 "http://web.comlab.ox.ac.uk/oucl/work/david.knezevic/",
		 "",
		 "PhD Student, Developer",
		 "<a href=\"http://web.comlab.ox.ac.uk/oucl/work/david.knezevic/\">Oxford University Computing Laboratory</a>");
?>

<?php
 dev_table_entry("images/blank.png",
		 "Hendrik van der Heijden",
		 "http://www.tu-harburg.de",
		 "",
		 "Student, Developer",
		 "<a href=\"http://www.tu-harburg.de\">Technical University Hamburg-Harburg</a>");
?>

<?php
 dev_table_entry("images/blank.png",
		 "Florian Prill",
		 "http://www.tu-harburg.de",
		 "",
		 "Student, Developer",
		 "<a href=\"http://www.tu-harburg.de\">Technical University Hamburg-Harburg</a>");
?>


</div>


<?php make_footer() ?>

</body>
</html>

