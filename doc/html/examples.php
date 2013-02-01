<?php $root=""; ?>
<?php require($root."navigation.php"); ?>

<html>
<head>
  <title>libMesh Examples</title>
  <?php load_style($root); ?>
</head>

<body>
<?php make_navigation("examples",$root)?>

<div class="content">

<h1>A Series of Example Programs</h1>
The following series of example programs have been
designed to get you started on the right foot.
They are arranged into categories based on which
library features they demonstrate.


<ul>

<li> Introduction <ol>

<li><L1><a href="introduction_ex1.php">Creation of a Mesh Object</a></L1></li>

<li><L1><a href="introduction_ex2.php">Defining a Simple System</a></L1></li>

<li><L1><a href="introduction_ex3.php">Solving a 2D Poisson Problem</a></L1></li>

<li><L1><a href="introduction_ex4.php">Solving a 2D or 3D Poisson Problem in Parallel</a></L1></li>

<li><L1><a href="introduction_ex5.php">Run-time Quadrature Rule Selection</a></L1></li>

</ol> </li>


<li> Systems of Equations <ol>

<li><L1><a href="systems_of_equations_ex1.php">Stokes Equations</a></L1></li>

<li><L1><a href="systems_of_equations_ex2.php">Unsteady Nonlinear Navier-Stokes</a></L1></li>

<li><L1><a href="systems_of_equations_ex3.php">Navier-Stokes with SCALAR Lagrange Multiplier</a></L1></li>

<li><L1><a href="systems_of_equations_ex4.php">Linear Elastic Cantilever</a></L1></li>

<li><L1><a href="systems_of_equations_ex5.php">Linear Elastic Cantilever with Constraint</a></L1></li>

<li><L1><a href="systems_of_equations_ex6.php">3D Linear Elastic Cantilever</a></L1></li>

</ol> </li>


<li> Transient Systems <ol>

<li><L1><a href="transient_ex1.php">Solving a Transient Linear System in Parallel</a></L1></li>

<li><L1><a href="transient_ex2.php">The Newmark System and the Wave Equation</a></L1></li>

</ol> </li>


<li> Adaptivity <ol>

<li><L1><a href="adaptivity_ex1.php">Solving 1D PDE Using Adaptive Mesh Refinement</a></L1></li>

<li><L1><a href="adaptivity_ex2.php">Solving a Transient System with Adaptive Mesh Refinement</a></L1></li>

<li><L1><a href="adaptivity_ex3.php">Laplace's Equation in an L-Shaped Domain</a></L1></li>

<li><L1><a href="adaptivity_ex4.php">Solving the Biharmonic Equation</a></L1></li>

<li><L1><a href="adaptivity_ex5.php">Periodic Boundary Conditions</a></L1></li>

</ol> </li>


<li> Eigenproblems <ol>

<li><L1><a href="eigenproblems_ex1.php">Solving an Eigen Problem</a></L1></li>

<li><L1><a href="eigenproblems_ex2.php">Solving a generalized Eigen Problem</a></L1></li>

<li><L1><a href="eigenproblems_ex3.php">Can you "hear the shape" of a drum?</a></L1></li>

</ol> </li>


<li> Subdomains <ol>

<li><L1><a href="subdomains_ex1.php">Solving on a Subdomain</a></L1></li>

<li><L1><a href="subdomains_ex2.php">Subdomain-Restricted Variables</a></L1></li>

</ol> </li>


<li> FEMSystem Framework <ol>

<li><L1><a href="fem_system_ex1.php">Unsteady Navier-Stokes Equations with FEMSystem</a></L1></li>

<li><L1><a href="fem_system_ex2.php">Nonlinear Elasticity with FEMSystem</a></L1></li>

</ol> </li>


<li> Reduced Basis <ol>

<li><L1><a href="reduced_basis_ex1.php">Certified Reduced Basis Method</a></L1></li>

<li><L1><a href="reduced_basis_ex2.php">Successive Constraint Method</a></L1></li>

<li><L1><a href="reduced_basis_ex3.php">Transient Reduced Basis Problem</a></L1></li>

<li><L1><a href="reduced_basis_ex4.php">Empirical Interpolation Method</a></L1></li>

<li><L1><a href="reduced_basis_ex5.php">Reduced Cantilever Problem</a></L1></li>

<li><L1><a href="reduced_basis_ex6.php">Heat Transfer on a Curved Domain in 3D</a></L1></li>

</ol> </li>


<li> Adjoints <ol>

<li><L1><a href="adjoints_ex1.php">Adjoint-based Goal Oriented Refinement</a></L1></li>

<li><L1><a href="adjoints_ex2.php">Adjoint-based Parameter Sensitivities</a></L1></li>

<li><L1><a href="adjoints_ex3.php">Adjoint-based coupled coupled Stokes + Convection Diffusion Goal Oriented Refinement</a></L1></li>

<li><L1><a href="adjoints_ex4.php">Adjoint-based Goal Oriented Refinement</a></L1></li>

<li><L1><a href="adjoints_ex5.php">SolutionHistory, General Localized Vectors and Unsteady Adjoints</a></L1></li>

</ol> </li>


<li> Vector-Valued Finite Elements <ol>

<li><L1><a href="vector_fe_ex1.php">Uncoupled Poisson Problem</a></L1></li>

<li><L1><a href="vector_fe_ex2.php">Unsteady Navier-Stokes with FEMSystem</a></L1></li>

</ol> </li>


<li> Miscellaneous <ol>

<li><L1><a href="miscellaneous_ex1.php">Infinite Elements for the Wave Equation</a></L1>

<li><L1><a href="miscellaneous_ex2.php">Complex Numbers and the "FrequencySystem"</a></L1></li>

<li><L1><a href="miscellaneous_ex3.php">2D Laplace-Young Problem Using Nonlinear Solvers</a></L1></li>

<li><L1><a href="miscellaneous_ex4.php">Using a Shell Matrix</a></L1></li>

<li><L1><a href="miscellaneous_ex5.php">Interior Penalty Discontinuous Galerkin</a></L1></li>

<li><L1><a href="miscellaneous_ex6.php">Meshing with Triangle and Tetgen</a></L1></li>

<li><L1><a href="miscellaneous_ex7.php">Variational Inequalities and PetscDMNonlinearSolver</a></L1></li>

<li><L1><a href="miscellaneous_ex8.php">Pointclould-based Meshfree Interpolation</a></L1></li>

</ol> </li>

</ul>

</div>

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
