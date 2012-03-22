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
<div class="fragment">
  <pre>svn checkout https://libmesh.svn.sourceforge.net/svnroot/libmesh/trunk/libmesh </pre>
</div>
If you would like to contribute to the project you will need a SourceForge
developer account, or you can contribute patches. To create a patch from a
modified SVN tree simply do:
<div class="fragment">
  <pre>svn diff &gt; patch </pre>
</div>
in the top-level directory. You can then submit the file <code>patch</code>.


<h2>Repository Organization</h2>
  The base repository for the <code>libMesh</code> project is <code>https://libmesh.svn.sourceforge.net/svnroot/libmesh </code>
  which holds the active development branch and previously tagged versions.  The top-level repository is organized as follows:
<div class="fragment">
  <pre>
  ./branches --> Branches for major feature additions.
  ./tags     --> Tagged versions of the project.
  ./trunk    --> The active development source tree.
  </pre>
</div>


<h2>Tagged Versions</h2>
<h3>Creating a new Tag</h3>
We often create tagged versions of the repository.  These tags often correspond to a released version of the code, or perhaps a stable snapshot in time we might want to easily revert to.  Creating a tagged version is simple:
<div class="fragment">
  <pre>
  $ svn copy https://libmesh.svn.sourceforge.net/svnroot/libmesh/trunk \
             https://libmesh.svn.sourceforge.net/svnroot/libmesh/tags/libmesh-x.y.z
  </pre>
</div>
will create a tagged version <code>libmesh-x.y.z</code>.

<h3>Checking out a Specified Tag</h3>
It is equally simple to check out a tagged version of the code:
<div class="fragment">
  <pre>
  $ svn checkout https://libmesh.svn.sourceforge.net/svnroot/libmesh/tags/libmesh-x.y.z/libmesh \
                 libmesh-x.y.z
  </pre>
</div>
will checkout the tagged version <code>libmesh-x.y.z</code> and place it in the local directory <code>libmesh-x.y.z</code>

<h2>Reference</h2>
For a comprehensive discussion of subversion, its features, and examples of its usage read <a href="http://svnbook.red-bean.com"><i>Version Control with</i> Subversion</a>.

</div>

<?php make_footer() ?>

</body>
</html>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
