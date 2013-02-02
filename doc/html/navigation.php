<?php
function load_style($root)
{ 
  // Do it the php way
  require($root."detect.php");
  
  // Call the detection script to set definitions
  detect();
  
  // Now check for CSS compatibility
  if (NW_MOZ_VERSION || NW_IS_GECKO || NW_IS_KONQ || NW_IS_OPERA || NW_IS_MAC)
    {
      echo "<link href='gecko.css' type='text/css' rel='stylesheet' media='screen'>\n";
    }

  else if (NW_IS_IE > 5)
    {
      echo "<link href='layout.css' type='text/css' rel='stylesheet' media='screen'>\n";
    }

  else if (NW_IS_NN <= 4)
    {
      echo "<link href='orig.css' type='text/css' rel='stylesheet' media='screen'>\n";
    }
}
?>



<?php function make_navigation($mode,$root) { ?>
 <div id="navAlpha">
    <?php echo "<a class=\"L1\" href = \"", $root, "index.php\" title=\"Main Page\">Home</a><BR>"; ?>    

    <?php echo "<a class=\"L1\" href = \"", $root, "publications.php\">Publications</a><BR>"; ?>
	<?php if ($mode=="publications") {?>
              <?php echo "<a class=\"L2\" href = \"", $root, "publications.php#Dissertations\" title=\"Dissetations\">Dissertations</a><BR>"; ?>
              <?php echo "<a class=\"L2\" href = \"", $root, "publications.php#Preprints\" title=\"Dissetations\">Preprints</a><BR>"; ?>
              <?php echo "<a class=\"L2\" href = \"", $root, "publications.php#Articles\" title=\"Dissetations\">Articles</a><BR>"; ?>
	<?php } ?>
 
    <?php echo "<a class=\"L1\" href = \"", $root, "presentations.php\">Presentations</a><BR>"; ?>

    <?php echo "<a class=\"L1\" href = \"", $root, "support.php\">Support</a><BR>"; ?>
	
    <?php echo "<a class=\"L1\" href = \"", $root, "installation.php\" title=\"Download & Installation Instructions\">Download</a><BR>"; ?>
        <?php if ($mode=="download") { ?>
              <?php echo "<a class=\"L2\" href = \"", $root, "installation.php\" title=\"Installation\">Installation</a><BR>"; ?>
              <!-- 
	        <?php echo "<a class=\"L2\" href = \"", $root, "subversion.php\" title=\"SVN Repository\">SVN Repository</a><BR>"; ?>
              -->
              <?php echo "<a class=\"L2\" href = \"", $root, "https://sourceforge.net/projects/libmesh/files/libmesh\" title=\"Source Files\">Source Files</a><BR>"; ?>
        <?php } ?>
	
    <?php echo "<a class=\"L1\" href = \"", $root, "examples.php\">Examples</a><BR>"; ?>
    <!--
    There are too many examples now to expand them all in the menu...
    <?php if (ereg("^ex[0-9]+|examples",$mode))
            {
              for ($i=0; $i<30; $i++)
                {
                  make_example_subs($i, $root, $mode);
                }
            } ?>
    -->

    <?php echo "<a class=\"L1\" href = \"", $root, "developers.php\">Developers</a><BR>"; ?>

    <?php echo "<a class=\"L1\" href = \"", $root, "doxygen/index.php\" title=\"C++ Class Documentation\">Class Docs</a><BR>"; ?>
	<?php if ($mode=="documentation") {?>
              <?php echo "<a class=\"L2\" href = \"", $root, "doxygen/namespaces.php\" title=\"Namespaces\">Namespaces</a><BR>"; ?>
              <?php echo "<a class=\"L2\" href = \"", $root, "doxygen/classes.php\" title=\"Classes\">Classes</a><BR>"; ?>
              <?php echo "<a class=\"L2\" href = \"", $root, "doxygen/files.php\" title=\"Files\">Files</a><BR>"; ?>
	<?php } ?>
 
    <?php echo "<a class=\"L1\" href = \"http://sourceforge.net/mail/?group_id=71130\" title=\"Sourceforge's Mailing List Page\">Mailing Lists</a><BR>" ?>

    <?php echo "<a class=\"L1\" href = \"", $root, "wiki/index.php/Gallery\" title=\"Results From Applications Which Use libMesh\">Gallery</a><BR>"; ?>

    <?php echo "<a class=\"L1\" href = \"", $root, "wiki\" title=\"External Wiki Page\">Wiki</a><BR>"; ?>
</div>


<?php } ?>



<?php function make_example_subs($num,$root,$mode) { ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex", $num, ".php\" title=\"Example $num\">Example $num</a><BR>\n"; ?>
    <?php if ($mode=="ex$num") { ?>
            <?php echo "<a class=\"L3\" href = \"", $root, "ex", $num, ".php#comments   \" title=\"Commented Code\">Comments</a><BR>\n"; ?>
	    <?php echo "<a class=\"L3\" href = \"", $root, "ex", $num, ".php#nocomments \" title=\"Source Code\">Source</a><BR>\n"; ?>
	    <?php echo "<a class=\"L3\" href = \"", $root, "ex", $num, ".php#output     \" title=\"Console Output\">Console Output</a><BR>\n"; ?>
    <?php } ?>
<?php } ?>

<?php function make_footer() { ?>

<div class="search">
  <!-- Google CSE Search Box Begins  -->
  <form action="http://www.google.com/cse" id="cse-search-box">
    <input type="hidden" name="cx" value="012399064992308798676:e_bxmf5nqrq">
    <input type="text" name="q" size="25">
    <input type="submit" name="sa" value="Search">
  </form>
  <script type="text/javascript" src="http://www.google.com/coop/cse/brand?form=cse-search-box&lang=en"></script>
  <!-- Google CSE Search Box Ends -->
</div>

<hr>
<i>Site Created By:</i>
<!-- <a href="mailto:libmesh-users@lists.sourceforge.net">libMesh Developers</a><br> -->
<a href="http://libmesh.sourceforge.net/developers.php">libMesh Developers</a><br>
<?php 
   // Newer versions of PHP require you to set the default timezone or a warning message is printed
   date_default_timezone_set('UTC');
?>
<i><?php echo "Last modified: ".date( "F d Y H:i:s", getlastmod() ); echo " UTC";?></i>

<br>
<br>
<i>Hosted By:</i>
<br>
<A href="http://sourceforge.net">
<IMG src="http://sourceforge.net/sflogo.php?group_id=71130&amp;type=1" width="88" height="31" border="0" alt="SourceForge.net Logo"></A>

<!-- Google Analytics stuff -->
<script type="text/javascript">

  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-24978333-1']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();

</script>
									    
<?php } ?>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>
