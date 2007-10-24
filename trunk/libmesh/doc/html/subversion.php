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


<h2>Repository Organization</h2>
  The base repository for the <code>libMesh</code> project is <pre>https://libmesh.svn.sourceforge.net/svnroot/libmesh </pre>
  which holds the active development branch and previously tagged versions.  The top-level repostory is organized as follows:
<div class="fragment">
  <pre>
  ./branches --> Branches for major feature additions.
  ./tags     --> Tagged versions of the project.
  ./trunk    --> The active development source tree.
  </pre>
</div>

<br>  
  

<h2>Tagged Versions</h2>
<h3>Creating a new Tag</h3>
We often create tagged versions of the repository.  These tags often correspond to a released version of the code, or perhaps a stable snapshot in time we might want to easily revert to.  Creating a tagged version is simple:
<div class="fragment">
  <pre>
  $ svn copy https://libmesh.svn.sourceforge.net/svnroot/libmesh/trunk \
             https://libmesh.svn.sourceforge.net/svnroot/libmesh/tags/libmesh-x.y.z
  </pre>
</div>
<br>
will create a tagged version <code>libmesh-x.y.z</code>.

<h3>Cheking out a Specified Tag</h3>
It is equally simple to check out a tagged version of the code:
<div class="fragment">
  <pre>
  $ svn checkout https://libmesh.svn.sourceforge.net/svnroot/libmesh/tags/libmesh-x.y.z/libmesh libmesh-x.y.z
  </pre>
</div>
<br>
will checkout the tagged version <code>libmesh-x.y.z</code> and place it in the local directory <code>libmesh-x.y.z</code>
</div>
<br>

<?php make_footer() ?>

</body>
</html>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
