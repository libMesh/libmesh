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
<h1>This is the starting point for applications</h1>
<ul>
  <li>    Foo </li>
  <li>    Bar </li>
</ul>

<h2>Natural Convection</h2>
<ul>
  <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/adv_flows_project.pdf">Ben's Adavnced Flows Class Project</a> </li>
  <li> Bar </li>
</ul>

<h2>Sub-Heading 2 </h2>
<ul>
  <li> Foo </li>
  <li> Bar </li>
</ul>

</div>

<!--
<div id="navBeta">
</div>
-->

<?php make_footer() ?>

</body>
</html>
