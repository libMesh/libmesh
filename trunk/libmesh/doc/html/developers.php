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
  <?php echo "Name: <a href=\"$homepage\">$name</a> <br>"; ?>
  <?php echo "Title(s): $titles <br>"; ?>
  <?php echo "Institution(s): $institutions"; ?>
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
	        "University of Texas at Austin");
?>

<br>

<?php
dev_table_entry("images/jwpeterson.jpg",
                "John W. Peterson",
                "http://www.cfdlab.ae.utexas.edu/~peterson",
		"PhD Student, libMesh Developer",
		"University of Texas at Austin");
?>

<?php
dev_table_entry("images/mikeando.jpg",
                "Michael L. Anderson",
                "http://www.cfdlab.ae.utexas.edu/~michaela",
		"Post-doctoral Fellow, libMesh Developer",
		"University of Texas at Austin, University of Western Australia");
?>


<?php
dev_table_entry("images/bbarth.jpg",
                "William L. Barth",
                "http://www.cfdlab.ae.utexas.edu/~bbarth",
		"PhD Student, libMesh Developer",
		"University of Texas at Austin");
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