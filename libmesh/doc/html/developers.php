<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Developers</title>
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
<h1>libMesh Developers</h1>

<?php
dev_table_entry("images/benkirk.jpg",
                "Benjamin S. Kirk",
                "http://www.cfdlab.ae.utexas.edu/~benkirk",
                "PhD Student, libMesh Project Manager",
	        "<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

<br>

<?php
dev_table_entry("images/jwpeterson.jpg",
                "John W. Peterson",
                "http://www.cfdlab.ae.utexas.edu/~peterson",
		"PhD Student, libMesh Developer",
		"<a href=\"http://www.utexas.edu\">University of Texas at Austin</a>");
?>

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

</div>


<?php make_footer() ?>

</body>
</html>

<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>