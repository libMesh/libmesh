<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>Application 1</title>
  <?php load_style($root); ?>
</head>

<body>

<?php make_navigation("applications",$root)?>

<div class="content">
<h1>Application 1: Why You Are A Panty-Waist</h1>
</div>

<div class="content">
<h2>Sub-heading 1</h2>
<ul>
<li> Foo </li>
<li> Bar </li>
</ul>
</div>

<div class="content">
<h2>Sub-heading 2</h2>
<ul>
<li> Foo </li>
<li> Bar </li>
</ul>
</div>

<?php make_footer() ?>
</body>
</html>
