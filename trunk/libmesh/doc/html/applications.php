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

<h2>Sub-Heading 1</h2>
<ul>
  <li> Foo </li>
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
