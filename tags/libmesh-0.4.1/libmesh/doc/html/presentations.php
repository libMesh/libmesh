<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Presentations</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("presentations",$root)?>

<div class="content">
<h1>Presentations</h1>
<ul>
  <li> <h2>Ben's Presentations</h2> </li>
  <ul>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/adv_flows_project.pdf">Adavnced Flows Class Project</a> (natural convection, 2002)</li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/siam_geo.pdf">SIAM Geosciences 2003</a></li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/libmesh_uwa03.pdf">libMesh Presentation</a> (from the June 2003 UWA Shortcourse)</li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/applications_uwa03.pdf">AMR Applications</a> (from the June 2003 UWA Shortcourse)</li>
  </ul>
  <li> USNCCM VII </li>
  <li> ICIAM </li>
  <li> SAVIAC </li>
</ul>
</div>

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