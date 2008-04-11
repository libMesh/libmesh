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
<h1>Please use the following citation to reference libMesh</h1>
<ul>
<li>
  B. Kirk, J. W. Peterson, R. H. Stogner, and G. F. Carey,
<code>libMesh</code><i>: A
  C++ Library for Parallel Adaptive Mesh Refinement/Coarsening Simulations.</i>
<a href="http://dx.doi.org/10.1007/s00366-006-0049-3">Engineering with Computers</a>,
vol. 22, no. 3--4, pp. 237--254, 2006. 	 (<a href="http://cfdlab.ae.utexas.edu/~benkirk/libmesh-ewc-preprint.pdf">preprint</a>)
<pre>
@Article{libMeshPaper,
  author = {B.~Kirk and J.~W.~Peterson and R.~H.~Stogner and G.~F.~Carey},
  title = {{\texttt{libMesh}: A C++ Library for Parallel Adaptive Mesh
           Refinement/Coarsening Simulations}},
  journal = {Engineering with Computers},
  volume = {22},
  number = {3--4},
  pages = {237--254},
  year = {2006},
  note = {\url{http://dx.doi.org/10.1007/s00366-006-0049-3}}
}
</pre>  
</li>
</ul>

<br>
<h1>Publications by libMesh Developers</h1>

<ul>


<li> Benjamin S. Kirk, 
<i>Adaptive Finite Element Simulation of Flow and Transport Applications on Parallel Computers.</i>
<a href="http://www.cfdlab.ae.utexas.edu/~benkirk/dissertation.pdf">PhD Dissertation</a>, 
The University of Texas at Austin, May 2007. 
</li>

<br>
<li> John W. Peterson,
<i>A Numerical Investigation of Benard Convection in Small Aspect Ratio Containers.</i>
<a href="http://www.cfdlab.ae.utexas.edu/~peterson/masters.pdf">Masters Report</a>, 
The University of Texas at Austin, August 2004. 
</li>

<br>
<li>
  Benjamin S. Kirk, Graham F. Carey, 
  <i>A Parallel, Adaptive Finite Element Scheme for Modeling Chemotactic Biological Systems</i>, submitted. (<a href="http://cfdlab.ae.utexas.edu/~benkirk/2008_chemotaxis_CNME.pdf">preprint</a>)
</li>

<br>
<li>
  Benjamin S. Kirk, Graham F. Carey, 
  <i>Development and Validation of a SUPG Finite Element Scheme for the Compressible Navier-Stokes Equations using a Modified Inviscid Flux Discretization.</i> <a href="http://dx.doi.org/10.1002/fld.1635">IJNMF</a>, vol. 57, no. 3, pp 265--293, May 2008. (<a href="http://cfdlab.ae.utexas.edu/~benkirk/2007_fins_IJNMF.pdf">preprint</a>)
</li>

<br>
<li>
J.W. Peterson, G.F. Carey, D.J. Knezevic, and B.T. Murray, <i>Adaptive
  finite element methodology for tumor angiogenesis modeling.</i>
  <a href="http://dx.doi.org/10.1002/nme.1802">
  IJNME</a>, vol. 69, no. 6, pp. 1212--1238, 2007.
<pre>
@Article{Tumor,
  author = {J.~W.~Peterson and G.~F.~Carey and D.~J.~Knezevic and B.~T.~Murray},
  title = {{Adaptive finite element methodology for tumor angiogenesis modeling}},
  journal = {Int.~J.~Numer.~Meth.~Eng.},
  volume = {69},
  number  = {6},
  pages = {1212--1238},
  year = {2007},
  note = {\url{http://dx.doi.org/10.1002/nme.1802}}
}
</pre>
</li>

<br>



<li>
G. F. Carey, W. Barth, B. Kirk, and J. W. Peterson,
<i>Parallel CFD for Flow and Transport Applications Including Unstructured and Adaptive Grids.</i>
  In <a href="http://ceani.ulpgc.es/pcfd04/">
  Proceedings of Parallel CFD 2004: Multidisciplinary Applications</a>,
  G. Winter, A. Ecer, J. Periaux, N. Satofuka and P. Fox (Eds), (Amsterdam,
  The Netherlands), Elsevier Science B.V., Oct 2005. ISBN: 0444520244.
</li>

<br>
<li>
R. H. Stogner and G. F. Carey,
<i>C<sup>1</sup> macroelements in adaptive finite element methods.</i>
<a href="http://doi.wiley.com/10.1002/nme.1912">IJNME</a> vol 70, issue 9, Pages 1076 - 1095, 2006.
</li>

<br>
<li>
Carey GF, Anderson M, Carnes B, and Kirk B, <i>Some aspects of adaptive grid
  technology related to boundary and interior layers.</i>
  <a href="http://dx.doi.org/10.1016/j.cam.2003.09.036"> J. Comput. Appl. Math.</a>
    166(1):55--86, ISSN 0377-0427, 2004.
</li>

<br>
<li>
Dreyer D, Petersen S, and O. von Estorff,
<i> Effectiveness and robustness of
improved infinite elements for exterior acoustics.</i>
<a href="http://dx.doi.org/10.1016/j.cma.2005.01.019">
CMAME</a>
 195(29-32):3591--3607, 2006.
</li>

<br>
<li>
Petersen S, Dreyer D, and O. von Estorff,
<i>Assessment of finite and spectral
element shape functions for efficient iterative simulations of interior
acoustics.</i>
<a href="http://dx.doi.org/10.1016/j.cma.2006.01.008">
CMAME</a>
Volume 195, Issues 44-47, Pages 6463-6478, 2006.
</li>
        
</ul>

<h1>Publications and Preprints Using libMesh</h1>

<ul>
<li>
Graham F. Carey, William Barth, Juliette A. Woods, Benjamin S. Kirk, Michael L. Anderson,
Sum Chow, and Wolfgang Bangerth, 
<i>Modelling error and constitutive relations in simulation of flow and transport.</i>
<a href="http://dx.doi.org/10.1002/fld.797">
IJNMF</a> Volume 46, Issue 12, Pages 1211 - 1236, 2004.
</li>


<br>
<li>
Jeremiah J. Marichalar, William C. Rochelle, Benjamin S. Kirk, and Charles H. Campbell, 
<i>Boundary Layer/Streamline Surface Catalytic Heating Predictions on Space Shuttle Orbiter.</i>
<a href="http://www.aiaa.org/content.cfm?pageid=322&lupubid=25&sItem=6">Journal of Spacecraft and Rockets</a>, 
Volume 43,  Number 6, Pages 1202--1215, 2006.
</li>


<br>
<li>
Michael Schindler, Peter Talkner, and Peter Hänggi,
<i>Computing stationary free-surface shapes in microfluidics.</i>
<a href="http://link.aip.org/link/?PHFLE6/18/103303/1">Physics of Fluids</a>, 
Volume 18, 2006.
</li>

<br>
<li>
Paul Simedrea, Luca Antiga, and David A. Steinman,
<i>Towards a New Framework for Simulating Magnetic Resonance Imaging.</i>
<a href="http://cscbc2006.cs.queensu.ca/assets/documents/Papers/paper108.pdf">
First Canadian Student Conference on Biomedical Computing (CSCBC)</a>, 2005
</li>

<br>
<li>
Jose Camata, Alvaro Coutinho, and Graham Carey,
<i>Numerical Evaluation of the LCD Method Implemented in the libMesh Library.</i>
<a href="http://www.inf.ufes.br/~avalli/papers/2005/CIL-0692.pdf.gz">CILAMCE</a> 2005.
</li>


<br>
<li>
Matthew Anderson and Jung-Han Kimn,
<i>A Numerical Approach to Space-Time Finite Elements for the Wave Equation.</i>
<a href="http://arxiv.org/pdf/gr-qc/0601099">arXiv.org</a>, 2006.
</li>


<br>
<li>
Stefano Berrone and Endre Suliz,
<i>Two-sided a posteriori error bounds for incompressible quasi-Newtonian flows.</i>
Dipartimento di Matematica Politecnico di Torino, <a href="http://calvino.polito.it/ricerca/2006/2006.html">calvino.polito.it/ricerca</a>, 2006. 
</li>


<br>
<li>
M.P. Luthi,
<i>A Full Ice Stream Model for Jakobshavn Isbrae.</i>
<a href="http://www.cosis.net/abstracts/EGU2007/02503/EGU2007-J-02503.pdf?PHPSESSID=f00857a95461bcae8e598e9edf0f1ba8">Geophysical Research Abstracts</a>, Vol 9., 2007.
</li>



</ul>




<h1>Miscellaneous</h1>

<ul>
<li>  A general <a href="howto/howto.pdf">HOWTO</a> document by M. Luthi containing some hints
and programming tips for writing effective libMesh programs. </li>

<li>  A <a href="xda_format/xda_format.pdf">description</a> of the XDA file format used by libMesh. </li>

<li> Texas Advanced Computing Center <a href="http://www.tacc.utexas.edu/general/news/archive/20040112_01.php">press release</a> commemorating the launch of the Lonestar cluster. </li>

<li>A <a href="http://ondrej.certik.cz/libmesh/fem.ps">description</a> of the Newmark System class by Ondrej Certik.</li>

</ul>
</div>

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
