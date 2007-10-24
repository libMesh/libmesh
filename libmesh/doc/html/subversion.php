<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>Download and Installation</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("download",$root)?>

<div class="content">
<h1>Subversion Repository Information</h1>
You may  access the SVN source tree for the latest code. You can get access
to the SVN repository via:
<br>

<div class="fragment">
  <pre>svn checkout https://libmesh.svn.sourceforge.net/svnroot/libmesh/trunk/libmesh </pre>
</div>

<br>
If you would like to contribute to the project you will need a SourceForge
developer account, or you can contribute patches. To create a patch from a
modified SVN tree simply do:
<br>
<div class="fragment">
  <pre>svn diff &gt; patch </pre>
</div>

<br>
in the top-level directory. You can then submit the file <code>patch</code>.


</div>

<?php make_footer() ?>

</body>
</html>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
