<?php function load_style($root) { ?>
   <script language="javascript" src="detect.js"></script>
   <script>
    if (is_gecko)
      {
        // Use special Gecko CSS style for all CSS-enabled types
        document.write("<link href='gecko.css' type='text/css' rel='stylesheet' media='screen'>\n")
      }
    else if (is_ie5up)
      {
        document.write("<link href='layout.css' type='text/css' rel='stylesheet' media='screen'>\n")
      }
    else if (is_nav && (is_major <= 4))
      {
        document.write("<link href='orig.css' type='text/css' rel='stylesheet' media='screen'>\n")
      }
    else
      {
        // document.write() statements to create static HTML markup for earlier versions
      }
     </script>

<!-- <style type="text/css" media="screen"> -->
<!-- <?php echo "@import \"".$root."layout.css\";"; ?> -->
<!-- </style> -->

<?php } ?>



<?php function make_navigation($mode,$root) { ?>
 <div id="navAlpha">
    <?php echo "<a class=\"L1\" href = \"", $root, "index.php\" title=\"Main Page\">Home</a><BR>"; ?>
    
    <?php echo "<a class=\"L1\" href = \"", $root, "examples.php\">Examples</a><BR>"; ?>

    <?php if (ereg("^ex[1-9]+|examples",$mode)) { ?>

    <?php for ($i=1; $i<7; $i++)
          {
            make_example_subs($i, $root, $mode);
          } ?>
      

    <?php make_example_subs(9, $root,$mode) ?>
    <?php make_example_subs(10,$root,$mode) ?>
    
    <?php } ?>

    <?php echo "<a class=\"L1\" href = \"", $root, "applications.php\">Applications</a><BR>"; ?>
    <?php if ($mode=="applications") { ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "application1.php\" title=\"Application 1\">Application 1</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "application2.php\" title=\"Application 2\">Application 2</a><BR>"; ?>
    <?php } ?>
	
    <?php echo "<a class=\"L1\" href = \"", $root, "presentations.php\">Presentations</a><BR>"; ?>
    <?php echo "<a class=\"L1\" href = \"", $root, "publications.php\">Publications</a><BR>"; ?>
    <?php echo "<a class=\"L1\" href = \"", $root, "developers.php\">Developers</a><BR>"; ?>
    <?php echo "<a class=\"L1\" href = \"", $root, "doxygen/index.html\">Class Docs</a><BR>"; ?>
    <?php echo "<a class=\"L1\" href = \"", $root, "installation.php\">Download</a><BR>" ?>
    <?php if ($mode=="download") { ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "installation.php\" title=\"Installation\">Installation</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "http://sourceforge.net/project/showfiles.php?group_id=71130\" title=\"Source Files\">Source Files</a><BR>"; ?>
    <?php } ?>
	
    <?php echo "<a class=\"L1\" href = \"", $root, "http://sourceforge.net/mail/?group_id=71130\">Mailing Lists</a><BR>" ?>
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
<hr>
<i>Site Created By:</i>
<a href="mailto:libmesh-users@lists.sourceforge.net">libMesh Developers</a><br>
<i><?php echo "Last modified: ".date( "F d Y H:i:s.", getlastmod() ); ?></i>

<br>
<br>
<i>Hosted By:</i>
<br>
<A href="http://sourceforge.net">
<IMG src="http://sourceforge.net/sflogo.php?group_id=71130&amp;type=1" width="88" height="31" border="0" alt="SourceForge.net Logo"></A>
									    
<?php } ?>


<?php if (0) { ?>
# Local Variables:
# mode: html
# End:
<?php } ?>