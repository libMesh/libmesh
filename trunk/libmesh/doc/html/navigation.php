<?php function load_style($root) { ?>
<style type="text/css" media="screen"><?php echo "@import \"".$root."layout.css\";"; ?></style>
<?php } ?>



<?php function make_navigation($mode,$root) { ?>
 <div id="navAlpha">
    <?php echo "<a class=\"L1\" href = \"", $root, "index.php\" title=\"Main Page\">Home</a><BR>"; ?>
    
    <?php echo "<a class=\"L1\" href = \"", $root, "examples.php\">Examples</a><BR>"; ?>
    <?php if ($mode=="examples") { ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex1.php\" title=\"Example 1\">Example 1</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex2.php\" title=\"Example 2\">Example 2</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex3.php\" title=\"Example 3\">Example 3</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex4.php\" title=\"Example 4\">Example 4</a><BR>"; ?>
    <?php echo "<a class=\"L2\" href = \"", $root, "ex5.php\" title=\"Example 5\">Example 5</a><BR>"; ?>
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
    <?php echo "<a class=\"L1\" href = \"", $root, "http://sourceforge.net/project/showfiles.php?group_id=71130\">Download</a><BR>" ?>
    <?php echo "<a class=\"L1\" href = \"", $root, "http://sourceforge.net/mail/?group_id=71130\">Mailing Lists</a><BR>" ?>
</div>
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