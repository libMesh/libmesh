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
vol. 22, no. 3--4, pp. 237--254, 2006.
(<a href="http://cfdlab.ae.utexas.edu/~benkirk/libmesh-ewc-preprint.pdf">preprint</a>)
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
<h1>Dissertations and Theses</h1>
<ul>
<br>
<li> Roy H. Stogner,
<i>Parallel Adaptive C1 Macro-Elements for Nonlinear Thin Film and Non-Newtonian Flow Problems.</i>
<a href="http://www.cfdlab.ae.utexas.edu/~roystgnr/dissertation/dissertation-stogner.pdf">PhD Dissertation</a>, 
The University of Texas at Austin, August 2008. 
</li>

<br>
<li> John W. Peterson, 
<i>Parallel Adaptive Finite Element Methods for Problems in Natural Convection.</i>
<a href="http://www.cfdlab.ae.utexas.edu/~peterson/diss.pdf">PhD Dissertation</a>, 
The University of Texas at Austin, May 2008. 
</li>

<br>
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

</ul>

<br>
<h1>Publications utilizing libMesh</h1>

<ol>

<li>
Benjamin S. Kirk,
<i>A Multidimensional Thermal Analysis to Assess Modeling Error
  in High-Speed Wind Tunnel Heat Transfer Data Reduction Schemes.</i>
<a href="http://www.aiaa.org/content.cfm?pageid=318">AIAA Journal of
  Thermophysics and Heat Transfer</a>,
vol. 23, no. 1, pp 186--191, January, 2009.
<pre>
@Article{Kirk_2009,
   author = {B.~S.~Kirk},
    title = {{Multidimensional Assessment of Modeling Error in Typical
              High-Speed Wind-Tunnel Heat-Transfer Data-Reduction Schemes}},
  journal = {AIAA Journal of Thermophysics and Heat Transfer},
     year = 2009,
    month = jan,
   volume = 23,
   number = 1,
    pages = {186--191}
}
</pre>
</li>


<br>
<li>
  Benjamin S. Kirk, Graham F. Carey, 
  <i>A Parallel, Adaptive Finite Element Scheme for
    Modeling Chemotactic Biological Systems.</i>,
  <a href="http://dx.doi.org/10.1002/cnm.1173">CNME</a>,
  accepted. (<a href="http://cfdlab.ae.utexas.edu/~benkirk/2008_chemotaxis_CNME.pdf">preprint</a>)
<pre>
@Article{Kirk_2008b,
   author = {Benjamin S. Kirk and Graham F. Carey},
    title = {{A parallel, adaptive finite element scheme for
              modeling chemotactic biological systems}},
  journal = {Communications in Numerical Methods in Engineering},
     year = 2008,
     note = {\url{http://dx.doi.org/10.1002/cnm.1173}}
}
</pre>  
</li>

<br>
<li>
  Benjamin S. Kirk, Graham F. Carey, 
  <i>Development and Validation of a SUPG Finite Element Scheme
    for the Compressible Navier-Stokes Equations using a Modified
    Inviscid Flux Discretization.</i> <a href="http://dx.doi.org/10.1002/fld.1635">IJNMF</a>,
    vol. 57, no. 3, pp 265--293, May 2008.
    (<a href="http://cfdlab.ae.utexas.edu/~benkirk/2007_fins_IJNMF.pdf">preprint</a>)
<pre>
@Article{Kirk_2008,
   author = {Benjamin S. Kirk and Graham F. Carey},
    title = {{Development and validation of a SUPG finite element scheme for the compressible
              Navier-Stokes equations using a modified inviscid flux discretization}},
  journal = {International Journal for Numerical Methods in Fluids},
     year = 2008,
   volume = 57,
   number = 3,
    pages = {265--293},
     note = {\url{http://dx.doi.org/10.1002/fld.1635}}
}
</pre>  
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
R. H. Stogner and G. F. Carey,
<i>C<sup>1</sup> macroelements in adaptive finite element methods.</i>
<a href="http://doi.wiley.com/10.1002/nme.1912">IJNME</a> vol 70, issue 9, Pages 1076 - 1095, 2006.
<pre>
@Article{Stogner_2007,
   author = {R.~H.~Stogner and G.~F.~Carey},
    title = {{$C^1$} macroelements in adaptive finite element methods},
  journal = {Int. J. Numerical Methods in Engineering},
   volume = 70,
   number = 9,
    pages = "1076--1095",
    month = {May},
     year = 2007
}
</pre>  
</li>

<br>
<li>
Dreyer D, Petersen S, and O. von Estorff,
<i> Effectiveness and robustness of
improved infinite elements for exterior acoustics.</i>
<a href="http://dx.doi.org/10.1016/j.cma.2005.01.019">
CMAME</a>
 195(29-32):3591--3607, 2006.
<pre>
@Article{Dreyer_2006,
   author = "D. Dreyer and S. Petersen and O. von Estorff",
    title = "Effectiveness and robustness of improved infinite elements for exterior acoustics",
  journal = "Computer Methods in Applied Mechanics and Engineering",
   volume = "195",
   number = "29--32",
    pages = "3591--3607",
     year = "2006",
     issn = "0045-7825",
     note = {\url{"http://dx.doi.org/10.1016/j.cma.2005.01.019"}}
}
</pre>
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
<pre>
@Article{Petersen_2006,
    author = {Petersen, S. and Dreyer, D. and {von Estorff}, O.},
    title  = {Assessment of finite and spectral element shape functions
              for efficient iterative simulations of interior acoustics},
   journal = {Computer Methods in Applied Mechanics and Engineering},
      year = {2006},
    volume = {195},
    number = {44--47},
     pages = {6463--6478}
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
<pre>
@InProceedings{Carey_2005,
     author = {G.~F.~Carey and W.~Barth and B.~Kirk and J.~W.~Peterson},
      title = {{Parallel CFD for Flow and Transport Applications
                Including Unstructured and Adaptive Grids}},
  booktitle = {Proceedings of Parallel CFD 2004:  Multidisciplinary Applications, G.~Winter,
               A.~Ecer, J.~Periaux, N.~Satofuka and P.~Fox (Eds)},
  publisher = {Elsevier Science B.V.},
    address = {Amsterdam, The Netherlands},
      month = {Oct},
       year = {2005},
       note = {ISBN: 0444520244}
}
</pre>
</li>

<br>
<li>
Carey GF, Anderson M, Carnes B, and Kirk B, <i>Some aspects of adaptive grid
  technology related to boundary and interior layers.</i>
  <a href="http://dx.doi.org/10.1016/j.cam.2003.09.036"> J. Comput. Appl. Math.</a>
    166(1):55--86, ISSN 0377-0427, 2004.
<pre>
@Article{Carey_2004a,
    author = {G.~F. Carey and M. Anderson and B. Carnes and B. Kirk},
     title = {Some aspects of adaptive grid technology related to boundary and interior layers},
   journal = {J. Comput. Appl. Math.},
    volume = {166},
    number = {1},
      year = {2004},
      issn = {0377-0427},
     pages = {55--86},
      note = {\url{http://dx.doi.org/10.1016/j.cam.2003.09.036}}
}
</pre>
</li>

        
<li>
Graham F. Carey, William Barth, Juliette A. Woods, Benjamin S. Kirk, Michael L. Anderson,
Sum Chow, and Wolfgang Bangerth, 
<i>Modelling error and constitutive relations in simulation of flow and transport.</i>
<a href="http://dx.doi.org/10.1002/fld.797">
IJNMF</a> Volume 46, Issue 12, Pages 1211 - 1236, 2004.
<pre>
@Article{Carey_2004,
  author  = {Graham F. Carey and William Barth and Juliette A. Woods
             and Ben Kirk and Michael L. Anderson and Sum Chow 
             and Wolfgang Bangerth},
  title   = {Modeling error and constitutive relations in simulation
             of flow and transport},
  journal = {Int. J. Numerical Methods in Fluid},
  year    = {2004},
  volume  = {46},
  number  = {12},
  pages   = {1211--1236}
}
</pre>  
</li>


<br>
<li>
Jeremiah J. Marichalar, William C. Rochelle, Benjamin S. Kirk, and Charles H. Campbell, 
<i>Boundary Layer/Streamline Surface Catalytic Heating Predictions on Space Shuttle Orbiter.</i>
<a href="http://www.aiaa.org/content.cfm?pageid=322&lupubid=25&sItem=6">Journal of Spacecraft and Rockets</a>, 
Volume 43,  Number 6, Pages 1202--1215, 2006.
<pre>
@Article{Marichalar_2006,
   author = {J.~J.~Marichalar and W.~C.~Rochelle and B.~S.~Kirk and C.~H.~Campbell},
    title = {{Boundary Layer/Streamline Surface Catalytic Heating Predictions on Space Shuttle Orbiter}},
  journal = {Journal of Spacecraft and Rockets},
     year = 2006,
   volume = 43,
   number = 6,
    pages = {1202--1215}
}
</pre>

</li>


<br>
<li>
Michael Schindler, Peter Talkner, and Peter Hanggi,
<i>Computing stationary free-surface shapes in microfluidics.</i>
<a href="http://link.aip.org/link/?PHFLE6/18/103303/1">Physics of Fluids</a>, 
Volume 18, 2006.
<pre>
@Article{Schindler_2006,
   author = {{Schindler}, M. and {Talkner}, P. and {H{\"a}nggi}, P.},
    title = "{Computing stationary free-surface shapes in microfluidics}",
  journal = {Physics of Fluids},
   eprint = {arXiv:physics/0511217},
     year = 2006,
    month = oct,
   volume = 18,
   number = 10,
    pages = {103303},
     note = {\url{http://dx.doi.org/10.1063/1.2361291}}
}
</pre>
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
Journal of Computational Physics, Volume 226, Issue 1, 10 September
2007, Pages 466-476, <a href="http://arxiv.org/pdf/gr-qc/0601099">preprint at arXiv.org</a>
<pre>
@Article{Anderson_2007,
   author = {{Anderson}, M. and {Kimn}, {J.-H.}},
    title = "{A numerical approach to space-time finite elements for the wave equation}",
  journal = {Journal of Computational Physics},
   eprint = {arXiv:gr-qc/0601099},
     year = 2007,
    month = sep,
   volume = 226,
    pages = {466--476},
     note = {\url{http://dx.doi.org/10.1016/j.jcp.2007.04.021}}
}
</pre>
</li>


<br>
<li>
Stefano Berrone and Endre Suli,
<i>Two-sided a posteriori error bounds for incompressible quasi-Newtonian flows.</i>
Dipartimento di Matematica Politecnico di Torino,
<a href="http://calvino.polito.it/ricerca/2006/2006.html">calvino.polito.it/ricerca</a>, 2006.
<pre>
@Article{Berrone_2008,
   author = {S.~Berrone and E.~S\"{u}li},
    title = {{Two-sided a posteriori error bounds for incompressible quasi-Newtonian flows}},
  journal = {IMA Journal of Numerical Analysis},
     year = 2008,
   volume = 28,
   number = 2,
    pages = {382--421},
     note = {\url{http://dx.doi.org/10.1093/imanum/drm017}}
}
</pre>  
</li>


<br>
<li>
M.P. Luthi,
<i>A Full Ice Stream Model for Jakobshavn Isbrae.</i>
<a href="http://www.cosis.net/abstracts/EGU2007/02503/EGU2007-J-02503.pdf?PHPSESSID=f00857a95461bcae8e598e9edf0f1ba8">Geophysical Research Abstracts</a>, Vol 9., 2007.
</li>


<br>
<li>
M. Auf der Maur, M. Povolotskyi, F. Sacconia and A. Di Carlo,
<i>TiberCAD: A new multiscale simulator for electronic and optoelectronic
  devices</i>. <a href="http://dx.doi.org/doi:10.1016/j.spmi.2007.03.011">Superlattices
  and Microstructures</a>, Volume 41, Issues 5-6, May-June 2007, Pages 381-385.
<pre>
@Article{AufderMaur_2007,
   author = {{Auf der Maur}, M. and {Povolotskyi}, M. and {Sacconi}, F. and {di Carlo}, A.},
    title = "{TIBERCAD: A new multiscale simulator for electronic and optoelectronic devices}",
  journal = {Superlattices and Microstructures},
     year = 2007,
    month = may,
   volume = 41,
    pages = {381--385},
     note = {\url{http://dx.doi.org/10.1016/j.spmi.2007.03.011}}
}
</pre>
</li>

<br>
<li>
Fabio Galbusera, Margherita Cioffi, and Manuela T. Raimondi,
<i>An in silico bioreactor for simulating laboratory experiments in tissue engineering</i>.
<a href="http://dx.doi.org/10.1007/s10544-008-9164-9">Biomedical Microdevices</a>,
Volume 10, Number 4, Pages 547-554, August, 2008.
<pre>
@Article{Galbusera_2008,
   author = {Fabio Galbusera and Margherita Cioffi and Manuela T. Raimondi},
    title = {{An in silico bioreactor for simulating laboratory experiments in tissue engineering}},
  journal = {Biomedical Microdevices},
     year = 2008,
   volume = 10,
   number = 4,
    month = aug,
    pages = {547--554}
}
</pre>
</li>


<br>
<li>
Hyung Taek Ahn, Mikhail Shashkov, Mark A. Christon,
<i>The moment-of-fluid method in action</i>.
Communications in Numerical Methods in Engineering, Vol. 25, Issue 10, pp. 1009--1018.
Published <a href="http://dx.doi.org/10.1002/cnm.1135">online</a>, July 2008.
<pre>
@Article{Ahn_2008,
   author = {{Ahn}, H.~T. and {Shashkov}, M. and M.~A.~Christon},
    title = "{The moment-of-fluid method in action}",
  journal = {Communications in Numerical Methods in Engineering},
     year = 2008,
    month = jul,
   volume = {25},
   number = {10},
    pages = {1009--1018}
}
</pre>
</li>

<br>
<li>
Marina Piccinelli, Lorenzo Botti, Bogdan Ene-Iordache, Andrea Remuzzi,
Alessandro Veneziani, Luca Antiga
<i>Link between vortex structures and voronoi diagram in cerebral aneurysms</i>
Journal of Biomechanics, Volume 41, Supplement 1, July 2008, Page S12
<pre>
@Article{Piccinelli_2008,
   author = {Marina Piccinelli and Lorenzo Botti and Bogdan Ene-Iordache
             and Andrea Remuzzi and Alessandro Veneziani and Luca Antiga},
    title = {{Link between vortex structures and voronoi diagram in cerebral aneurysms}},
  journal = {Journal of Biomechanics},
     year = 2008,
    month = jul,
   volume = {41, Supplement 1},
    pages = {S12}
}
</pre>
</li>


<br>
<li>
Maik Brinkmeier, Udo Nackenhorst, Steffen Petersen, Otto von Estorff
<i>A finite element approach for the simulation of tire rolling noise</i>
Journal of Sound and Vibration, Volume 309, Issues 1-2, 8 January
2008, Pages 20-39
<pre>
@Article{Brinkmeier_2008,
    author = {Brinkmeier, M. and Nackenhorst, U. and Petersen, S. and {von Estorff}, O.},
     title = {A Numerical Model for the Simulation of Tire Rolling Noise},
   journal = {Journal of Sound and Vibration},
      year = {2008},
     month = jan,
    volume = {309},
    number = {1--2},
     pages = {20--39},
      note = {\url{http://dx.doi.org/10.1016/j.jsv.2006.11.040}}
}
</pre>
</li>


<br>
<li>
Tim Kroger and Tobias Preusser
<i>Stability of the 8-tetrahedra shortest-interior-edge partitioning method</i>
Numerische Mathematik
Volume 109, Number 3, May 2008, pp. 435-457.
<a href="http://dx.doi.org/10.1007/s00211-008-0148-8">DOI</a>
<pre>
@Article{Kroger_2008,
   author = {Tim Kr\"{o}ger and Tobias Preusser},
    title = {{Stability of the 8-tetrahedra shortest-interior-edge partitioning method}},
  journal = {Numerische Mathematik},
     year = 2008,
    month = may,
   volume = 109,
   number = 3,
    pages = {435--457},
     note = {\url{http://dx.doi.org/10.1007/s00211-008-0148-8}}
}
</pre>
</li>


<br>
<li>
Koen Van Canneyt, Radoslav Kaminsky, Lorenzo Botti, Luca Antiga, Jan
Tordoir, Pascal Verdonck, Sunny Eloot
<i>Can a kinked arterio-venous graft cause instable flow? A
patient-specific validated computational study</i>
Journal of Biomechanics, Volume 41, Supplement 1, July 2008, Page S210
<pre>
@Article{VanCanneyt_2008,
   author = "Koen Van Canneyt and Radoslav Kaminsky and Lorenzo Botti and
             Luca Antiga and Jan Tordoir and Pascal Verdonck and Sunny Eloot",
    title = "CAN A KINKED ARTERIO-VENOUS GRAFT CAUSE INSTABLE FLOW?
             A PATIENT-SPECIFIC VALIDATED COMPUTATIONAL STUDY",
  journal = "Journal of Biomechanics",
   volume = "41",
   number = "Supplement 1",
    pages = "S210",
     year = "2008",
    month = jul,
     note = "Abstracts of the 16th Congress, European Society of Biomechanics",
     issn = "0021-9290",
     note = {\url{http://dx.doi.org/10.1016/S0021-9290(08)70210-0}}
}
</pre>
</li>

<br>
<li>
Lorenzo Botti, Marina Piccinelli, Bogdan Ene-Iordache, Andrea Remuzzi,
Luca Antiga
<i>An open source parallel AMR FE solver for biological flows based on
the libMesh C++ library</i>
Journal of Biomechanics, Volume 41, Supplement 1, July 2008, Page S211
<pre>
@Article{Botti_2008,
   author = "Lorenzo Botti and Marina Piccinelli and Bogdan Ene-Iordache and Andrea Remuzzi and Luca Antiga",
    title = "AN OPEN SOURCE PARALLEL {AMR} {FE} SOLVER FOR BIOLOGICAL FLOWS BASED ON THE {LibMesh C++} LIBRARY",
  journal = "Journal of Biomechanics",
   volume = "41",
   number = "Supplement 1",
    pages = "S211",
     year = "2008",
    month = jul,
     note = "Abstracts of the 16th Congress, European Society of Biomechanics",
     issn = "0021-9290",
     note = {\url{http://dx.doi.org/10.1016/S0021-9290(08)70211-2}}
}
</pre>
</li>

<br>
<li>
Jan Biermann, Otto von Estorff, Steffen Petersen, Christina Wenterodt
<i>Higher order finite and infinite elements for the solution of Helmholtz problems</i>
Computer Methods in Applied Mechanics and Engineering, Volume 198,
Issues 13-14, 1 March 2009, Pages 1171-1188
<pre>
@Article{Biermann_2009,
   author = "Jan Biermann and Otto von Estorff and Steffen Petersen and Christina Wenterodt",
    title = {{Higher order finite and infinite elements for the solution of Helmholtz problems}},
  journal = "Computer Methods in Applied Mechanics and Engineering",
   volume = "198",
   number = "13--14",
    pages = "1171--1188",
     year = "2009",
     note = "HOFEM07 - International Workshop on High-Order Finite Element Methods, 2007",
     issn = "0045-7825",
     note = {\url{http://dx.doi.org/10.1016/j.cma.2008.11.009}}
}
</pre>  
</li>

<br>
<li>
Matthew F. Barone, Irina Kalashnikova, Daniel J. Segalman, Heidi K. Thornquist,
<i>Stable Galerkin reduced order models for linearized compressible flow</i>
Journal of Computational Physics, Volume 228, Issue 6, 1 April 2009,
Pages 1932-1946
<pre>
@Article{Barone_2009,
   author = "Matthew F. Barone and Irina Kalashnikova and Daniel J. Segalman and Heidi K. Thornquist",
    title = {{Stable Galerkin reduced order models for linearized compressible flow}},
  journal = "Journal of Computational Physics",
   volume = "228",
   number = "6",
    pages = "1932--1946",
     year = "2009",
    month = apr,
     issn = "0021-9991",
     note = {\url{http://dx.doi.org/10.1016/j.jcp.2008.11.015}}
}
</pre>
</li>


<br>
<li>
A.M.P. Valli, R.N. Elias, G.F. Carey and A.L.G.A.Coutinho,
<i>PID Adaptive Control of Incremental and Arclength Continuation in Nonlinear Applications</i>,
International Journal for Numerical Methods in Fluids,
Published
<a href="http://www.inf.ufes.br/~avalli/papers/2009/avalli_ijnmf09.pdf">Online</a>, Jan 20 2009,
<a href="http://www3.interscience.wiley.com/journal/121638902/abstract">abstract page</a>.
<pre>
@Article{Valli_2009,
   author = {A. M. P. Valli and R. N. Elias and G. F. Carey and A. L. G. A. Coutinho},
    title = {{PID adaptive control of incremental and arclength continuation in nonlinear applications}},
  journal = {International Journal for Numerical Methods in Fluids},
     year = 2009,
   volume = {Early View},
    note  = {\url{http://dx.doi.org/10.1002/fld.1998}}
}
</pre>
</li>


<br>
<li>
Adam C. Powell and Raymundo Arroyave
<i>Open source software for materials and process modeling</i>
JOM Journal of the Minerals, Metals and Materials Society,
Volume 60, Number 5 / May, 2008, <a href="http://www.tms.org/pubs/journals/jom/0805/powell-0805.html">link</a>.
<pre>
@Article{Powell_2008,
   author = {{Powell}, IV, A.~C. and {Arroyave}, R.},
    title = "{Open source software for materials and process modeling}",
  journal = {Journal of the Minerals, Metals and Materials Society},
     year = 2008,
    month = may,
   volume = 60,
   number = 5,
    pages = {32--39},
     note = {\url{http://dx.doi.org/10.1007/s11837-008-0057-4}}
}
</pre>
</li>

</ol>




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
