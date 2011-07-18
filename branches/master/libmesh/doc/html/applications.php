<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Applications</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("applications",$root)?>

<div class="content">
<h2>Here you will find presentations, organized by application area, that use libMesh</h2>
<ul>
  <li>    Foo </li>
  <li>    Bar </li>
</ul>

<h2>Natural Convection</h2>
<ul>
  <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/adv_flows_project.pdf">Ben's Adavnced Flows Class Project</a> </li>
</ul>

<h2>Exterior Acoustics</h2>
<ul>
  <li> <a href="http://www.mum.tu-harburg.de/ovestorff/projects/project_dd.html">Acoustic Radiation of Ocean Engineering Structures</a> </li>
  <li> <a href="http://www.mum.tu-harburg.de/ovestorff/projects/project_sp.html">Tire Rolling Noise</a> </li>
</ul>

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