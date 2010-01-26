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
J.&nbsp;W. Peterson, G.&nbsp;F. Carey, and B.&nbsp;T. Murray.
  Multi-Resolution Simulation of Double-Diffusive Convection in Porous
  Media.
  <em>International Journal of Numerical Methods for Heat &amp; Fluid
  Flow</em>, 20(1):37-65, 2010.
  <a href="http://dx.doi.org/10.1108/09615531011008118">http://dx.doi.org/10.1108/09615531011008118</a>.<br />
<pre>
@Article{Peterson_2010,
  author  = {J.~W.~Peterson and G.~F.~Carey and B.~T.~Murray},
  title   = {{Multi-Resolution Simulation of Double-Diffusive 
              Convection in Porous Media}},
  journal = {International Journal of Numerical Methods for Heat \& Fluid Flow},
  volume  = 20,
  number  = 1,
  year    = {2010},
  pages   = {37--65},
  note    = {\url{http://dx.doi.org/10.1108/09615531011008118}}
}
</pre>
</li>


<li>
V.&nbsp;S. Mahadevan, J.&nbsp;C. Ragusa, and V.&nbsp;Mousseau.
  Verification of multiphysics software: Space and time convergence
  studies for nonlinearly coupled applications.
  In <em>International Conference on Mathematics, Computaional Methods
  &amp; Reactor Physics (M&amp;C 2009)</em>, Saratoga Springs, New York, 2009.
  <a href="http://www.inl.gov/technicalpublications/Documents/4247162.pdf">http://www.inl.gov/technicalpublications/Documents/4247162.pdf</a>.<br />
<pre>
@Conference{Mahadevan_2009,
  author    = {V.~S.~Mahadevan and J.~C.~Ragusa and V.~Mousseau},
  title     = {Verification of Multiphysics software: Space and time
               convergence studies for nonlinearly coupled applications},
  booktitle = {International Conference on Mathematics, Computaional
               Methods \& Reactor Physics (M\&C 2009)},
  address   = {Saratoga Springs, New York},
  year      = {2009},
  note      = {\url{http://www.inl.gov/technicalpublications/Documents/4247162.pdf}}
}
</pre>
</li>


<li>
D.&nbsp;Gaston, G.&nbsp;Hansen, S.&nbsp;Kadioglu, D.&nbsp;A. Knoll, C.&nbsp;Newman, H.&nbsp;Park, C.&nbsp;Permann,
  and W.&nbsp;Taitano.
  Parallel multiphysics algorithms and software for computational
  nuclear engineering.
  <em>Journal of Physics: Conference Series</em>, 180:012012 (10pp), 2009.
  <a href="http://stacks.iop.org/1742-6596/180/012012">http://stacks.iop.org/1742-6596/180/012012</a>.<br />
<pre>
@Article{Gaston_2009,
  author = {D.~Gaston and G.~Hansen and S.~Kadioglu and D.~A.~Knoll
            and C.~Newman and H.~Park and C.~Permann and W.~Taitano},
  title  = {Parallel multiphysics algorithms and software for computational nuclear engineering},
  journal= {Journal of Physics: Conference Series},
  volume = {180},
  pages  = {012012 (10pp)},
  year   = {2009},
  note   = {\url{http://stacks.iop.org/1742-6596/180/012012}}
}
</pre>
</li>

<li>
M.&nbsp;P. L&uuml;thi.
  Transient response of idealized glaciers to climate variations.
  <em>Journal of Glaciology</em>, 55(193):918-930, 2009.<br />  
<pre>
@Article{Luethi2009,
 author  = {M.~P.~L\"{u}thi},
 title   = {Transient response of idealized glaciers to climate variations},
 journal = {Journal of Glaciology},
 year    = 2009,
 number  = 193,
 volume  = 55,
 pages   = {918--930}
}
</pre>
</li>


<li>
B.&nbsp;S. Kirk.
  Multidimensional Assessment of Modeling Error in Typical High-Speed
  Wind-Tunnel Heat-Transfer Data-Reduction Schemes.
  <em>AIAA Journal of Thermophysics and Heat Transfer</em>,
    23(1):186-191, January 2009.<br />
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
Benjamin&nbsp;S. Kirk and Graham&nbsp;F. Carey.
  A parallel, adaptive finite element scheme for modeling chemotactic
  biological systems.
  <em>Communications in Numerical Methods in Engineering</em>, 2008.
    <a href="http://dx.doi.org/10.1002/cnm.1173">http://dx.doi.org/10.1002/cnm.1173</a>.<br />
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
Benjamin&nbsp;S. Kirk and Graham&nbsp;F. Carey.
  Development and validation of a SUPG finite element scheme for the
  compressible Navier-Stokes equations using a modified inviscid flux
  discretization.
  <em>International Journal for Numerical Methods in Fluids</em>,
    57(3):265-293, 2008.
    <a href="http://dx.doi.org/10.1002/fld.1635">http://dx.doi.org/10.1002/fld.1635</a>.<br />
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

<li>
V.&nbsp;S. Mahadevan and J.&nbsp;C. Ragusa.
  High-order spatio-temporal schemes for coupled, multi-physics reactor
  simulations.
  Technical report, Idaho National Laboratory, September 2008.
  <a href="http://www.inl.gov/technicalpublications/Documents/4120524.pdf">http://www.inl.gov/technicalpublications/Documents/4120524.pdf</a>.<br />
<pre>
@TechReport{Mahadevan_2008,
  author      = {V.~S.~Mahadevan and J.~C.~Ragusa},
  title       = {High-order spatio-temporal schemes for coupled, multi-physics reactor simulations},
  institution = {Idaho National Laboratory},
  month       = {September},
  year        = {2008},
  note        = {\url{http://www.inl.gov/technicalpublications/Documents/4120524.pdf}}
}
</pre>
</li>

<br>
<li>
J.&nbsp;W. Peterson, G.&nbsp;F. Carey, D.&nbsp;J. Knezevic, and B.&nbsp;T. Murray.
  Adaptive finite element methodology for tumor angiogenesis
  modeling.
  <em>Int.&nbsp;J.&nbsp;Numer.&nbsp;Meth.&nbsp;Eng.</em>, 69(6):1212-1238, 2007.
    <a href="http://dx.doi.org/10.1002/nme.1802">http://dx.doi.org/10.1002/nme.1802</a>.<br />
<pre>
@Article{Peterson_2007,
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
R.&nbsp;H. Stogner and G.&nbsp;F. Carey.
  <em>C</em><sup>1</sup> macroelements in adaptive finite element methods.
    <em>Int. J. Numerical Methods in Engineering</em>, 70(9):1076-1095, May
      2007.<br />
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
D.&nbsp;Dreyer, S.&nbsp;Petersen, and O.&nbsp;von Estorff.
  Effectiveness and robustness of improved infinite elements for
  exterior acoustics.
  <em>Computer Methods in Applied Mechanics and Engineering</em>,
    195(29-32):3591-3607, 2006.
    <a href="http://dx.doi.org/10.1016/j.cma.2005.01.019">http://dx.doi.org/10.1016/j.cma.2005.01.019</a>.<br />
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
G.&nbsp;F. Carey, W.&nbsp;Barth, B.&nbsp;Kirk, and J.&nbsp;W. Peterson.
  Parallel CFD for Flow and Transport Applications Including
  Unstructured and Adaptive Grids.
  In <em>Proceedings of Parallel CFD 2004: Multidisciplinary
      Applications, G.&nbsp;Winter, A.&nbsp;Ecer, J.&nbsp;Periaux, N.&nbsp;Satofuka and P.&nbsp;Fox (Eds)</em>,
  Amsterdam, The Netherlands, October 2005. Elsevier Science B.V.
  ISBN: 0444520244.<br />
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
G.&nbsp;F. Carey, M.&nbsp;Anderson, B.&nbsp;Carnes, and B.&nbsp;Kirk.
  Some aspects of adaptive grid technology related to boundary and
  interior layers.
  <em>J. Comput. Appl. Math.</em>, 166(1):55-86, 2004.
    <a href="http://dx.doi.org/10.1016/j.cam.2003.09.036">http://dx.doi.org/10.1016/j.cam.2003.09.036</a>.<br />
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
Graham&nbsp;F. Carey, William Barth, Juliette&nbsp;A. Woods, Ben Kirk, Michael&nbsp;L.
  Anderson, Sum Chow, and Wolfgang Bangerth.
  Modeling error and constitutive relations in simulation of flow and
  transport.
  <em>Int. J. Numerical Methods in Fluid</em>, 46(12):1211-1236, 2004.<br />
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
J.&nbsp;J. Marichalar, W.&nbsp;C. Rochelle, B.&nbsp;S. Kirk, and C.&nbsp;H. Campbell.
  Boundary Layer/Streamline Surface Catalytic Heating Predictions on
  Space Shuttle Orbiter.
  <em>Journal of Spacecraft and Rockets</em>, 43(6):1202-1215, 2006.<br />
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
M.&nbsp;Schindler, P.&nbsp;Talkner, and P.&nbsp;H&auml;nggi.
  Computing stationary free-surface shapes in microfluidics.
  <em>Physics of Fluids</em>, 18(10):103303, October 2006.
    <a href="http://dx.doi.org/10.1063/1.2361291">http://dx.doi.org/10.1063/1.2361291</a>.<br />
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
M.&nbsp;Anderson and J.-H. Kimn.
  A numerical approach to space-time finite elements for the wave
  equation.
  <em>Journal of Computational Physics</em>, 226:466-476, September 2007.
    <a href="http://dx.doi.org/10.1016/j.jcp.2007.04.021">http://dx.doi.org/10.1016/j.jcp.2007.04.021</a>.<br />
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
S.&nbsp;Berrone and E.&nbsp;S&uuml;li.
  Two-sided a posteriori error bounds for incompressible
  quasi-Newtonian flows.
  <em>IMA Journal of Numerical Analysis</em>, 28(2):382-421, 2008.
    <a href="http://dx.doi.org/10.1093/imanum/drm017">http://dx.doi.org/10.1093/imanum/drm017</a>.<br />
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
M.P. L&uuml;thi,
<i>A Full Ice Stream Model for Jakobshavn Isbrae.</i>
<a href="http://www.cosis.net/abstracts/EGU2007/02503/EGU2007-J-02503.pdf?PHPSESSID=f00857a95461bcae8e598e9edf0f1ba8">Geophysical Research Abstracts</a>, Vol 9., 2007.
</li>


<br>
<li>
M.&nbsp;Auf der Maur, M.&nbsp;Povolotskyi, F.&nbsp;Sacconi, and A.&nbsp;di Carlo.
  TIBERCAD: A new multiscale simulator for electronic and
  optoelectronic devices.
  <em>Superlattices and Microstructures</em>, 41:381-385, May 2007.
    <a href="http://dx.doi.org/10.1016/j.spmi.2007.03.011">http://dx.doi.org/10.1016/j.spmi.2007.03.011</a>.<br />
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
Fabio Galbusera, Margherita Cioffi, and Manuela&nbsp;T. Raimondi.
  An in silico bioreactor for simulating laboratory experiments in
  tissue engineering.
  <em>Biomedical Microdevices</em>, 10(4):547-554, August 2008.<br />
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
H.&nbsp;T. Ahn, M.&nbsp;Shashkov, and M.&nbsp;A. Christon.
  The moment-of-fluid method in action.
  <em>Communications in Numerical Methods in Engineering</em>,
    25(10):1009-1018, July 2008.<br />
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
  Alessandro Veneziani, and Luca Antiga.
  Link between vortex structures and voronoi diagram in cerebral
  aneurysms.
  <em>Journal of Biomechanics</em>, 41, Supplement 1:S12, July 2008.<br />
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
M.&nbsp;Brinkmeier, U.&nbsp;Nackenhorst, S.&nbsp;Petersen, and O.&nbsp;von Estorff.
  A numerical model for the simulation of tire rolling noise.
  <em>Journal of Sound and Vibration</em>, 309(1-2):20-39, January 2008.
    <a href="http://dx.doi.org/10.1016/j.jsv.2006.11.040">http://dx.doi.org/10.1016/j.jsv.2006.11.040</a>.<br />
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
Tim Kr&ouml;ger and Tobias Preusser.
  Stability of the 8-tetrahedra shortest-interior-edge partitioning
  method.
  <em>Numerische Mathematik</em>, 109(3):435-457, May 2008.
    <a href="http://dx.doi.org/10.1007/s00211-008-0148-8">http://dx.doi.org/10.1007/s00211-008-0148-8</a>.<br />
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
Koen&nbsp;Van Canneyt, Radoslav Kaminsky, Lorenzo Botti, Luca Antiga, Jan Tordoir,
  Pascal Verdonck, and Sunny Eloot.
  Can a kinked arterio-venous graft cause instable flow? a
  patient-specific validated computational study.
  <em>Journal of Biomechanics</em>, 41(Supplement 1):S210, July 2008.
    <a href="http://dx.doi.org/10.1016/S0021-9290(08)70210-0">http://dx.doi.org/10.1016/S0021-9290(08)70210-0</a>.<br />
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
Lorenzo Botti, Marina Piccinelli, Bogdan Ene-Iordache, Andrea Remuzzi, and Luca
  Antiga.
  An open source parallel AMR FE solver for biological flows based
  on the LibMesh C++ library.
  <em>Journal of Biomechanics</em>, 41(Supplement 1):S211, July 2008.
    <a href="http://dx.doi.org/10.1016/S0021-9290(08)70211-2">http://dx.doi.org/10.1016/S0021-9290(08)70211-2</a>.<br />
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
Jan Biermann, Otto von Estorff, Steffen Petersen, and Christina Wenterodt.
  Higher order finite and infinite elements for the solution of
  Helmholtz problems.
  <em>Computer Methods in Applied Mechanics and Engineering</em>,
    198(13-14):1171-1188, 2009.
    <a href="http://dx.doi.org/10.1016/j.cma.2008.11.009">http://dx.doi.org/10.1016/j.cma.2008.11.009</a>.<br />
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
Matthew&nbsp;F. Barone, Irina Kalashnikova, Daniel&nbsp;J. Segalman, and Heidi&nbsp;K.
Thornquist.
Stable Galerkin reduced order models for linearized compressible
flow.
<em>Journal of Computational Physics</em>, 228(6):1932-1946, April
2009.
<a href="http://dx.doi.org/10.1016/j.jcp.2008.11.015">http://dx.doi.org/10.1016/j.jcp.2008.11.015</a>.<br />
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
A.&nbsp;M.&nbsp;P. Valli, R.&nbsp;N. Elias, G.&nbsp;F. Carey, and A.&nbsp;L. G.&nbsp;A. Coutinho.
  PID adaptive control of incremental and arclength continuation in
  nonlinear applications.
  <em>International Journal for Numerical Methods in Fluids</em>,
  61(11):1181-1200, 2009.
  <a href="http://dx.doi.org/10.1002/fld.1998">http://dx.doi.org/10.1002/fld.1998</a>.<br />
<pre>
@Article{Valli_2009,
   author = {A. M. P. Valli and R. N. Elias and G. F. Carey and A. L. G. A. Coutinho},
    title = {{PID adaptive control of incremental and arclength continuation in nonlinear applications}},
  journal = {International Journal for Numerical Methods in Fluids},
     year = 2009,
   volume = {61},
   number = {11},
    pages = {1181--1200},
    note  = {\url{http://dx.doi.org/10.1002/fld.1998}}
}
</pre>
</li>


<br>
<li>
A.&nbsp;C. Powell, IV and R.&nbsp;Arroyave.
  Open source software for materials and process modeling.
  <em>Journal of the Minerals, Metals and Materials Society</em>,
    60(5):32-39, May 2008.
    <a href="http://dx.doi.org/10.1007/s11837-008-0057-4">http://dx.doi.org/10.1007/s11837-008-0057-4</a>.<br />
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
