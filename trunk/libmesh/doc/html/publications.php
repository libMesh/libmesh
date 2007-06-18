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
<h1>Publications by LibMesh Developers</h1>

<ul>

<?php echo
"<li>libMesh: A C++ Library for Parallel Adaptive
 Mesh Refinement/Coarsening Simulations - Engineering with Computers.
      <ul>
	<li>
	  <a href=\"http://cfdlab.ae.utexas.edu/~benkirk/libmesh-ewc-preprint.pdf\">Preprint</a>
	</li>
	<li>
	  <a href=\"http://cfdlab.ae.utexas.edu/~peterson/ewc_journal_version.pdf\">Journal Version</a>
	</li>
      </ul>
</li>";
?>

<?php echo
"<li> Ben Kirk's <a href=\"http://www.cfdlab.ae.utexas.edu/~peterson/ben_dissertation.pdf\">
PhD Dissertation</a>. </li>";
?>

<?php echo
"<li> John Peterson's <a href=\"http://www.cfdlab.ae.utexas.edu/~peterson/masters.pdf\">
Masters Report</a>. </li>";
?>
</ul>

<h1>Publications and Preprints Using LibMesh</h1>

<ul>
<li>
Graham F. Carey, William Barth, Juliette A. Woods, Benjamin S. Kirk, Michael L. Anderson,
Sum Chow, and Wolfgang Bangerth, 
<i>Modelling error and constitutive relations in simulation of flow and transport.</i>
<a href="http://www3.interscience.wiley.com/cgi-bin/abstract/109747121/ABSTRACT?CRETRY=1&SRETRY=0">
IJNMF</a> Volume 46, Issue 12, Pages 1211 - 1236, 2004.
</li>

<li>
Michael Schindler, Peter Talkner, and Peter Hanggi,					     
<i>Computing stationary free-surface shapes in microfluidics.</i>
<a href="http://arxiv.org/pdf/physics/0511217">Physics Preprint arXiv</a>
</li>

<li>
Paul Simedrea, Luca Antiga, and David A. Steinman,
<i>Towards a New Framework for Simulating Magnetic Resonance Imaging.</i>
<a href="http://cscbc2006.cs.queensu.ca/assets/documents/Papers/paper108.pdf">
First Canadian Student Conference on Biomedical Computing (CSCBC)</a>, 2005
</li>

<li>
Jose Camata, Alvaro Coutinho, and Graham Carey,
<i>Numerical Evaluation of the LCD Method Implemented in the LibMesh Library.</i>
<a href="http://www.inf.ufes.br/~avalli/papers/2005/CIL-0692.pdf.gz">CILAMCE</a> 2005.
</li>
</ul>


<h1>Misc. Publications</h1>

<ul>
<?php echo
"<li>  A general <a href=\"howto/howto.pdf\">HOWTO</a>
document by M. Luthi containing some hints
and programming tips for writing effective LibMesh programs. </li>";
?>

<?php echo
"<li>  A <a href=\"xda_format/xda_format.pdf\">description</a>
of the XDA file format used by LibMesh. </li>";
?>

<?php echo
"<li> Texas Advanced Computing Center
<a href=\"http://www.tacc.utexas.edu/general/press/announcements/20040112_01.php\">
press release</a> commemorating the launch of the Lonestar cluster. </li>";
?>

<?php echo
"<li>A <a href=\"http://ondrej.certik.cz/libmesh/fem.ps\">description</a> of the Newmark
System class by Ondrej Certik.</li>";
?>

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
