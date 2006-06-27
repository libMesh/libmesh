<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Presentations</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("presentations",$root)?>

<div class="content">
<h1>Presentations</h1>
<ul>

  <li> <h2>Ben's Presentations</h2> </li>
  <ul>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/adv_flows_project.pdf">Adavnced Flows Class Project</a> (natural convection, 2002)</li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/siam_geo.pdf">SIAM Geosciences 2003</a></li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/libmesh_uwa03.pdf">libMesh Presentation</a> (from the June 2003 UWA Shortcourse)</li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/applications_uwa03.pdf">AMR Applications</a> (from the June 2003 UWA Shortcourse)</li>
    <li> <a href="http://cfdlab.ae.utexas.edu/~benkirk/seminar/talk.pdf">Finite Elements - Introduction and Applications</a> (April 2004 NASA/JSC technical seminar)</li>
  </ul>
  
  <li> <h2>USNCCM VII</h2> </li>
  <ul>
    <li> <a href="http://usnccm.sandia.gov/mslist/upload/IncomprssblCFD/2230_Final_abstract.pdf">Abstract </a> </li>
    <li> <a href="http://www.cfdlab.ae.utexas.edu/~peterson/usnccm.pdf">Presentation.</a> </li>
  </ul>

  <li> <h2>USNCCM VIII</h2> </li>
  <ul>
    <li> <a href="http://www.cfdlab.ae.utexas.edu/~benkirk/carter-USNCCMVIII.pdf">Compressible Flow Studies Using Parallel Adaptive Mesh Refinement</a> </li>
    <li> <a href="http://www.cfdlab.ae.utexas.edu/~benkirk/carey-USNCCMVIII.pdf"> Algorithms for Compressible Flows with Adaptive Mesh Refinement</a> </li>
    <li> <a href="http://www.cfdlab.ae.utexas.edu/~peterson/USNCCM8.pdf">         Adaptive Finite Element Simulations of Thermosolutal Convection in Porous Media</a> </li>
  </ul>
  
  <li> <h2>Finite Element Rodeo</h2> </li>
  <ul>
    <li> <a href="http://www.cfdlab.ae.utexas.edu/~peterson/Rodeo_2006.pdf">A Stabilized h-Adaptive Continuation Method for Double-Diffusive Convection in Porous Media</a> </li>
  </ul>

  
  <!--   add some blank space -->
  <h2> </h2> q
  <li> <a href="http://www.cfdlab.ae.utexas.edu/~peterson/scraper.pdf">Scraper Flow</a> Presentation with I. Schoegl.</li>
  <li> <a href="http://www.challenge.nm.org/FinalReports/08.pdf">Albuquerque Academy</a> Presentation on Stress Analysis in a Torus.</li>
  <li> <a href="http://www.saviac.org/74th_Symposium/abstracts/U040.htm">SAVIAC</a> conference abstract.</li>

</ul>
</div>

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
