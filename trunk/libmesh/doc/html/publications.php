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
  <li>Theses and dissertations</li>
  John's <a href="http://www.cfdlab.ae.utexas.edu/~peterson/masters.pdf">Masters Report</a>.


  <li>Web Links</li>
  
  <a href="howto/howto.pdf">HOWTO</a><br>
  XDA file <a href="xda_format/xda_format.pdf">description</a><br>
  <a href="http://www.tacc.utexas.edu/general/press/announcements/20040112_01.php">TACC</a> press release.<br>
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