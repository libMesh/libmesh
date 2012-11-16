<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex7",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "biharmonic.h"
        
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Print usage information if requested on command line
</div>

<div class ="fragment">
<pre>
        void print_help(int argc, char** argv);
        
        int main(int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
          if (on_command_line("--help"))
            print_help(argc, argv);
          else
            {
        #if !defined(LIBMESH_ENABLE_SECOND_DERIVATIVES)
              libmesh_example_assert(false, "--enable-second");
        #elif !defined(LIBMESH_ENABLE_PERIODIC)
              libmesh_example_assert(false, "--enable-periodic");
        #endif
        
</pre>
</div>
<div class = "comment">
This is a PETSc-specific solver
</div>

<div class ="fragment">
<pre>
              libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
              const int dim = command_line_value("dim",1);
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
              libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D/3D support");
        
              Biharmonic* biharmonic;
              Biharmonic::Create(&biharmonic);
              biharmonic-&gt;viewParameters();
              biharmonic-&gt;init();
              biharmonic-&gt;run();
              Biharmonic::Destroy(&biharmonic);
            }
          return 0;
        }
        
        
        
        
        
        void print_help(int, char** argv)
        {
          libMesh::out &lt;&lt; "This example solves the Cahn-Hillard equation with chemical potential f:\n"
        	       &lt;&lt; "    u_t = \\div(M(u)\\grad f(u))\n"
        	       &lt;&lt; "Here we have\n"
        	       &lt;&lt; "    u, -1 &lt;= u &lt;= 1        -- relative concentration (difference of two concentrations in a binary mixture) \n"
        	       &lt;&lt; "    M, M &gt;= 0              -- mobility of the mixture\n"
        	       &lt;&lt; "    f = \\delta E/\\delta u  -- variational derivative of the free energy functional E\n"
        	       &lt;&lt; "    E = \\int[\\kappa/2 |\\grac u|^ + g(u)]\n"
        	       &lt;&lt; "where the gradient term is the interfacial energy density with \\kappa quantifying the energy of the interface,\n"
        	       &lt;&lt; "and g(u) is the bulk energy density\n"
        	       &lt;&lt; "    g(u) = \\theta L(u) + \\theta_c W(u),\n"
        	       &lt;&lt; "L(u) is the (optional, in this model) logarithmic term corresponding to the entropy of the mixture:\n"
        	       &lt;&lt; "    L(u) = (\\theta/2)[(1+u)\\ln((1+u)/2) + (1-u)\\ln((1-u)/2)],\n"
        	       &lt;&lt; "where \\theta is related to the Boltzmann factor k_B T - a proxy for the absolute temperature T.\n"
        	       &lt;&lt; "L can be optionally approximated ('truncated') using a quadratic or a cubic polynomial on [-1,1]\n"
        	       &lt;&lt; "W(u) is the (optional, in this model) potential promoting demixing.  It can take the form of \n"
        	       &lt;&lt; "a 'double well' potential\n"
        	       &lt;&lt; "    W(u) = \\theta_c (u^4/4 - u^2/2),\n"
        	       &lt;&lt; "         or \n"
        	       &lt;&lt; "a 'double obstacle' potential\n"
        	       &lt;&lt; "    W(u) = (\\theta_c/2)(1-u^2),\n"
        	       &lt;&lt; "where \\theta_c is the critical 'temperature'.\n"
        	       &lt;&lt; "Finally, mobility M can be constant of 'degenerate', by which we mean that M is varying with u and \n"
        	       &lt;&lt; "vanishing (degenerating) whenever u reaches critical values +1 or -1:\n"
        	       &lt;&lt; "    M(u) = 1.0\n"
        	       &lt;&lt; "      or\n"
        	       &lt;&lt; "    M(u) = (1.0 - u^2)\n"
        	       &lt;&lt; "Degenerate mobility should generally be used only in conjunction with logarithmic free energy terms.\n\n"
        	       &lt;&lt; "The equation is solved on a periodic domain (in 1D, 2D or 3D)\n"
        	       &lt;&lt; "using a Galerkin formulation with C^1 elements approximating the H^2_{per} function space.\n\n"
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "COMPILING: "
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "Compile as follows (assuming libmesh has been built): \n"
        	       &lt;&lt; "METHOD=&lt;method&gt; make \n"
        	       &lt;&lt; "where &lt;method&gt; is dbg or opt.\n"
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "HELP:        "
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "Print this help message:\n"
        	       &lt;&lt; argv[0] &lt;&lt; " --help\n"
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "RUNNING:     "
        	       &lt;&lt; "\n-----------\n"
        	       &lt;&lt; "Run in serial with build METHOD &lt;method&gt; as follows:\n"
        	       &lt;&lt; "\n"
        	       &lt;&lt; argv[0] &lt;&lt; "\n"
        	       &lt;&lt; "               [--verbose] dim=&lt;1|2|3&gt; N=&lt;number_of_linear_elements&gt; \n"
        	       &lt;&lt; "               kappa=&lt;kappa_value&gt; growth=&lt;yes|no&gt; degenerate=&lt;yes|no&gt; [--cahn-hillard]                                           \n"
        	       &lt;&lt; "               [--netforce]  energy=&lt;double_well|double_obstacle|log_double_well|log_double_obstacle&gt;  log_truncation_order=&lt;2|3&gt; \n"
        	       &lt;&lt; "               theta=&lt;theta_value&gt; theta_c=&lt;theta_c_value&gt;                                                                        \n"
        	       &lt;&lt; "               initial_state=&lt;ball|rod|strip&gt; initial_center='x [y [z]]' initial_width=&lt;width&gt;                                    \n"
        	       &lt;&lt; "               min_time=&lt;initial_time&gt; max_time=&lt;final_time&gt; dt=&lt;timestep_size&gt; crank_nicholson_weight=&lt;between_0_and_1&gt;          \n"
        	       &lt;&lt; "               output_base=&lt;base_filename&gt; output_dt=&lt;output_timestep_size&gt; [--use-petsc-dm -snes_type virs]                      \n"
        	       &lt;&lt; "\n"
        	       &lt;&lt; argv[0] &lt;&lt; " --verbose \n"
        	       &lt;&lt; "is a pretty good start.\n"
        	       &lt;&lt; "\nModeling a 1D system with 2D or 3D (for a strip the second and third components of the center are immaterial):\n"
        	       &lt;&lt; argv[0]&lt;&lt; " --verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6\n"
        	       &lt;&lt; argv[0]&lt;&lt; " --verbose dim=2 N=64   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
        	       &lt;&lt; argv[0]&lt;&lt; " --verbose dim=3 N=32   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
        	       &lt;&lt; "\n"
        	       &lt;&lt; "Modeling a 2D system with 3D (for a rod the third component of the center is immaterial) \n"
        	       &lt;&lt; argv[0]&lt;&lt; " --verbose dim=2 N=64   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
        	       &lt;&lt; argv[0]&lt;&lt; " --verbose dim=3 N=32   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
        	       &lt;&lt; "\n"
        	       &lt;&lt; "A 3D system with an initial ball in the center\n"
        	       &lt;&lt; argv[0] &lt;&lt; " --verbose dim=3 N=32   initial_state=ball initial_center='0.5 0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n"
        	       &lt;&lt; "\n"
        	       &lt;&lt; "Add --use-petsc-dm -snes_type virs to run the variational inequality version that ensures the solution is between -1.0 and 1.0 at all times.\n\n"
        	       &lt;&lt; std::endl;
        }
        
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;biharmonic.h&quot;</FONT></B>
  
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> print_help(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
    <B><FONT COLOR="#A020F0">if</FONT></B> (on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--help&quot;</FONT></B>))
      print_help(argc, argv);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_ENABLE_SECOND_DERIVATIVES)
        libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-second&quot;</FONT></B>);
  #elif !defined(LIBMESH_ENABLE_PERIODIC)
        libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-periodic&quot;</FONT></B>);
  #endif
  
        libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;dim&quot;</FONT></B>,1);
  
        libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
        Biharmonic* biharmonic;
        <B><FONT COLOR="#5F9EA0">Biharmonic</FONT></B>::Create(&amp;biharmonic);
        biharmonic-&gt;viewParameters();
        biharmonic-&gt;init();
        biharmonic-&gt;run();
        <B><FONT COLOR="#5F9EA0">Biharmonic</FONT></B>::Destroy(&amp;biharmonic);
      }
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> print_help(<B><FONT COLOR="#228B22">int</FONT></B>, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;This example solves the Cahn-Hillard equation with chemical potential f:\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    u_t = \\div(M(u)\\grad f(u))\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Here we have\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    u, -1 &lt;= u &lt;= 1        -- relative concentration (difference of two concentrations in a binary mixture) \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    M, M &gt;= 0              -- mobility of the mixture\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    f = \\delta E/\\delta u  -- variational derivative of the free energy functional E\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    E = \\int[\\kappa/2 |\\grac u|^ + g(u)]\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;where the gradient term is the interfacial energy density with \\kappa quantifying the energy of the interface,\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;and g(u) is the bulk energy density\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    g(u) = \\theta L(u) + \\theta_c W(u),\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L(u) is the (optional, in this model) logarithmic term corresponding to the entropy of the mixture:\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    L(u) = (\\theta/2)[(1+u)\\ln((1+u)/2) + (1-u)\\ln((1-u)/2)],\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;where \\theta is related to the Boltzmann factor k_B T - a proxy for the absolute temperature T.\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L can be optionally approximated ('truncated') using a quadratic or a cubic polynomial on [-1,1]\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;W(u) is the (optional, in this model) potential promoting demixing.  It can take the form of \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;a 'double well' potential\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    W(u) = \\theta_c (u^4/4 - u^2/2),\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;         or \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;a 'double obstacle' potential\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    W(u) = (\\theta_c/2)(1-u^2),\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;where \\theta_c is the critical 'temperature'.\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Finally, mobility M can be constant of 'degenerate', by which we mean that M is varying with u and \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;vanishing (degenerating) whenever u reaches critical values +1 or -1:\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    M(u) = 1.0\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;      or\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    M(u) = (1.0 - u^2)\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Degenerate mobility should generally be used only in conjunction with logarithmic free energy terms.\n\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;The equation is solved on a periodic domain (in 1D, 2D or 3D)\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;using a Galerkin formulation with C^1 elements approximating the H^2_{per} function space.\n\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;COMPILING: &quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Compile as follows (assuming libmesh has been built): \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;METHOD=&lt;method&gt; make \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;where &lt;method&gt; is dbg or opt.\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;HELP:        &quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Print this help message:\n&quot;</FONT></B>
  	       &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --help\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;RUNNING:     &quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n-----------\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Run in serial with build METHOD &lt;method&gt; as follows:\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               [--verbose] dim=&lt;1|2|3&gt; N=&lt;number_of_linear_elements&gt; \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               kappa=&lt;kappa_value&gt; growth=&lt;yes|no&gt; degenerate=&lt;yes|no&gt; [--cahn-hillard]                                           \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               [--netforce]  energy=&lt;double_well|double_obstacle|log_double_well|log_double_obstacle&gt;  log_truncation_order=&lt;2|3&gt; \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               theta=&lt;theta_value&gt; theta_c=&lt;theta_c_value&gt;                                                                        \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               initial_state=&lt;ball|rod|strip&gt; initial_center='x [y [z]]' initial_width=&lt;width&gt;                                    \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               min_time=&lt;initial_time&gt; max_time=&lt;final_time&gt; dt=&lt;timestep_size&gt; crank_nicholson_weight=&lt;between_0_and_1&gt;          \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;               output_base=&lt;base_filename&gt; output_dt=&lt;output_timestep_size&gt; [--use-petsc-dm -snes_type virs]                      \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;is a pretty good start.\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\nModeling a 1D system with 2D or 3D (for a strip the second and third components of the center are immaterial):\n&quot;</FONT></B>
  	       &lt;&lt; argv[0]&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6\n&quot;</FONT></B>
  	       &lt;&lt; argv[0]&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=2 N=64   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n&quot;</FONT></B>
  	       &lt;&lt; argv[0]&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=3 N=32   initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-6 \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Modeling a 2D system with 3D (for a rod the third component of the center is immaterial) \n&quot;</FONT></B>
  	       &lt;&lt; argv[0]&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=2 N=64   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n&quot;</FONT></B>
  	       &lt;&lt; argv[0]&lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=3 N=32   initial_state=rod initial_center='0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;A 3D system with an initial ball in the center\n&quot;</FONT></B>
  	       &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; --verbose dim=3 N=32   initial_state=ball initial_center='0.5 0.5 0.5' initial_width=0.1 dt=1e-10 max_time=1e-6 \n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
  	       &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Add --use-petsc-dm -snes_type virs to run the variational inequality version that ensures the solution is between -1.0 and 1.0 at all times.\n\n&quot;</FONT></B>
  	       &lt;&lt; std::endl;
  }
  
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking miscellaneous_ex7-opt...
***************************************************************
* Running Example  mpirun -np 6 ./miscellaneous_ex7-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Biharmonic parameters:
verbose mode is on
mesh dimension           = 1
initial linear mesh size = 1024
kappa                    = 1
growth                   = 0
degenerate               = 0
Cahn-Hillard             = 0
netforce                 = 0
energy                   = 1
tol                      = 1e-08
theta                    = 0.001
theta_c                  = 1
log truncation           = 2
initial timestep size    = 1e-10
initial state:             strip
initial state center     = 0.5
initial state width      = 0.1
initial time (min_time)  = 0
integration time         = 1e-08
final time   (max_time)  = 1e-08
Crank-Nicholson weight   = 1
Output timestep          = 0
Output filename base:      bih.1
>>> Initializing Biharmonic
>>> Initializing Biharmonic::JR
<<< Initializing Biharmonic::JR
<<< Initializing Biharmonic
Writing state 1 at time 0 to file bih.1.e; output a total of 1 states so far
Solving for state 2, time 0
  NL step  0, |residual|_2 = 3.865471e+00
  NL step  1, |residual|_2 = 1.455971e-10
Writing state 2 at time 1.000000e-10 to file bih.1.e; output a total of 2 states so far
Solving for state 3, time 1.000000e-10
  NL step  0, |residual|_2 = 3.861870e+00
  NL step  1, |residual|_2 = 1.422037e-09
Writing state 3 at time 2.000000e-10 to file bih.1.e; output a total of 3 states so far
Solving for state 4, time 2.000000e-10
  NL step  0, |residual|_2 = 3.858634e+00
  NL step  1, |residual|_2 = 1.745300e-08
Writing state 4 at time 3.000000e-10 to file bih.1.e; output a total of 4 states so far
Solving for state 5, time 3.000000e-10
  NL step  0, |residual|_2 = 3.855483e+00
  NL step  1, |residual|_2 = 5.816269e-08
  NL step  2, |residual|_2 = 1.499954e-13
Writing state 5 at time 4.000000e-10 to file bih.1.e; output a total of 5 states so far
Solving for state 6, time 4.000000e-10
  NL step  0, |residual|_2 = 3.852464e+00
  NL step  1, |residual|_2 = 1.231523e-08
Writing state 6 at time 5.000000e-10 to file bih.1.e; output a total of 6 states so far
Solving for state 7, time 5.000000e-10
  NL step  0, |residual|_2 = 3.849501e+00
  NL step  1, |residual|_2 = 2.710877e-07
  NL step  2, |residual|_2 = 1.107187e-12
Writing state 7 at time 6.000000e-10 to file bih.1.e; output a total of 7 states so far
Solving for state 8, time 6.000000e-10
  NL step  0, |residual|_2 = 3.846621e+00
  NL step  1, |residual|_2 = 6.486216e-07
  NL step  2, |residual|_2 = 4.553978e-12
Writing state 8 at time 7.000000e-10 to file bih.1.e; output a total of 8 states so far
Solving for state 9, time 7.000000e-10
  NL step  0, |residual|_2 = 3.843784e+00
  NL step  1, |residual|_2 = 4.820632e-07
  NL step  2, |residual|_2 = 2.341594e-12
Writing state 9 at time 8.000000e-10 to file bih.1.e; output a total of 9 states so far
Solving for state 10, time 8.000000e-10
  NL step  0, |residual|_2 = 3.841009e+00
  NL step  1, |residual|_2 = 1.510397e-07
  NL step  2, |residual|_2 = 2.140622e-13
Writing state 10 at time 9.000000e-10 to file bih.1.e; output a total of 10 states so far
Solving for state 11, time 9.000000e-10
  NL step  0, |residual|_2 = 3.838268e+00
  NL step  1, |residual|_2 = 5.952732e-07
  NL step  2, |residual|_2 = 4.767653e-12
Writing state 11 at time 1.000000e-09 to file bih.1.e; output a total of 11 states so far
Solving for state 12, time 1.000000e-09
  NL step  0, |residual|_2 = 3.835577e+00
  NL step  1, |residual|_2 = 5.796059e-07
  NL step  2, |residual|_2 = 3.889766e-12
Writing state 12 at time 1.100000e-09 to file bih.1.e; output a total of 12 states so far
Solving for state 13, time 1.100000e-09
  NL step  0, |residual|_2 = 3.832915e+00
  NL step  1, |residual|_2 = 1.041354e-06
  NL step  2, |residual|_2 = 2.441486e-15
Writing state 13 at time 1.200000e-09 to file bih.1.e; output a total of 13 states so far
Solving for state 14, time 1.200000e-09
  NL step  0, |residual|_2 = 3.830295e+00
  NL step  1, |residual|_2 = 9.869192e-07
  NL step  2, |residual|_2 = 2.584606e-15
Writing state 14 at time 1.300000e-09 to file bih.1.e; output a total of 14 states so far
Solving for state 15, time 1.300000e-09
  NL step  0, |residual|_2 = 3.827700e+00
  NL step  1, |residual|_2 = 1.217420e-06
  NL step  2, |residual|_2 = 2.667277e-15
Writing state 15 at time 1.400000e-09 to file bih.1.e; output a total of 15 states so far
Solving for state 16, time 1.400000e-09
  NL step  0, |residual|_2 = 3.825140e+00
  NL step  1, |residual|_2 = 1.118917e-06
  NL step  2, |residual|_2 = 2.431359e-15
Writing state 16 at time 1.500000e-09 to file bih.1.e; output a total of 16 states so far
Solving for state 17, time 1.500000e-09
  NL step  0, |residual|_2 = 3.822604e+00
  NL step  1, |residual|_2 = 1.140332e-06
  NL step  2, |residual|_2 = 2.603906e-15
Writing state 17 at time 1.600000e-09 to file bih.1.e; output a total of 17 states so far
Solving for state 18, time 1.600000e-09
  NL step  0, |residual|_2 = 3.820098e+00
  NL step  1, |residual|_2 = 1.041043e-06
  NL step  2, |residual|_2 = 2.569107e-15
Writing state 18 at time 1.700000e-09 to file bih.1.e; output a total of 18 states so far
Solving for state 19, time 1.700000e-09
  NL step  0, |residual|_2 = 3.817613e+00
  NL step  1, |residual|_2 = 9.006636e-07
  NL step  2, |residual|_2 = 2.635896e-15
Writing state 19 at time 1.800000e-09 to file bih.1.e; output a total of 19 states so far
Solving for state 20, time 1.800000e-09
  NL step  0, |residual|_2 = 3.815156e+00
  NL step  1, |residual|_2 = 8.652038e-07
  NL step  2, |residual|_2 = 2.504349e-15
Writing state 20 at time 1.900000e-09 to file bih.1.e; output a total of 20 states so far
Solving for state 21, time 1.900000e-09
  NL step  0, |residual|_2 = 3.812718e+00
  NL step  1, |residual|_2 = 5.765407e-07
  NL step  2, |residual|_2 = 2.386284e-15
Writing state 21 at time 2.000000e-09 to file bih.1.e; output a total of 21 states so far
Solving for state 22, time 2.000000e-09
  NL step  0, |residual|_2 = 3.810305e+00
  NL step  1, |residual|_2 = 7.135885e-07
  NL step  2, |residual|_2 = 2.493672e-15
Writing state 22 at time 2.100000e-09 to file bih.1.e; output a total of 22 states so far
Solving for state 23, time 2.100000e-09
  NL step  0, |residual|_2 = 3.807909e+00
  NL step  1, |residual|_2 = 2.418441e-07
  NL step  2, |residual|_2 = 2.811973e-13
Writing state 23 at time 2.200000e-09 to file bih.1.e; output a total of 23 states so far
Solving for state 24, time 2.200000e-09
  NL step  0, |residual|_2 = 3.805536e+00
  NL step  1, |residual|_2 = 6.967777e-07
  NL step  2, |residual|_2 = 2.364348e-15
Writing state 24 at time 2.300000e-09 to file bih.1.e; output a total of 24 states so far
Solving for state 25, time 2.300000e-09
  NL step  0, |residual|_2 = 3.803179e+00
  NL step  1, |residual|_2 = 2.225682e-07
  NL step  2, |residual|_2 = 2.384040e-15
Writing state 25 at time 2.400000e-09 to file bih.1.e; output a total of 25 states so far
Solving for state 26, time 2.400000e-09
  NL step  0, |residual|_2 = 3.800842e+00
  NL step  1, |residual|_2 = 8.246108e-07
  NL step  2, |residual|_2 = 2.249607e-15
Writing state 26 at time 2.500000e-09 to file bih.1.e; output a total of 26 states so far
Solving for state 27, time 2.500000e-09
  NL step  0, |residual|_2 = 3.798522e+00
  NL step  1, |residual|_2 = 5.234375e-07
  NL step  2, |residual|_2 = 2.592257e-15
Writing state 27 at time 2.600000e-09 to file bih.1.e; output a total of 27 states so far
Solving for state 28, time 2.600000e-09
  NL step  0, |residual|_2 = 3.796220e+00
  NL step  1, |residual|_2 = 1.017306e-06
  NL step  2, |residual|_2 = 2.518285e-15
Writing state 28 at time 2.700000e-09 to file bih.1.e; output a total of 28 states so far
Solving for state 29, time 2.700000e-09
  NL step  0, |residual|_2 = 3.793932e+00
  NL step  1, |residual|_2 = 8.108586e-07
  NL step  2, |residual|_2 = 2.568864e-15
Writing state 29 at time 2.800000e-09 to file bih.1.e; output a total of 29 states so far
Solving for state 30, time 2.800000e-09
  NL step  0, |residual|_2 = 3.791663e+00
  NL step  1, |residual|_2 = 1.214613e-06
  NL step  2, |residual|_2 = 2.628195e-15
Writing state 30 at time 2.900000e-09 to file bih.1.e; output a total of 30 states so far
Solving for state 31, time 2.900000e-09
  NL step  0, |residual|_2 = 3.789407e+00
  NL step  1, |residual|_2 = 1.057661e-06
  NL step  2, |residual|_2 = 2.418453e-15
Writing state 31 at time 3.000000e-09 to file bih.1.e; output a total of 31 states so far
Solving for state 32, time 3.000000e-09
  NL step  0, |residual|_2 = 3.787167e+00
  NL step  1, |residual|_2 = 1.390170e-06
  NL step  2, |residual|_2 = 2.301078e-15
Writing state 32 at time 3.100000e-09 to file bih.1.e; output a total of 32 states so far
Solving for state 33, time 3.100000e-09
  NL step  0, |residual|_2 = 3.784941e+00
  NL step  1, |residual|_2 = 1.260362e-06
  NL step  2, |residual|_2 = 2.490627e-15
Writing state 33 at time 3.200000e-09 to file bih.1.e; output a total of 33 states so far
Solving for state 34, time 3.200000e-09
  NL step  0, |residual|_2 = 3.782730e+00
  NL step  1, |residual|_2 = 1.534471e-06
  NL step  2, |residual|_2 = 2.592762e-15
Writing state 34 at time 3.300000e-09 to file bih.1.e; output a total of 34 states so far
Solving for state 35, time 3.300000e-09
  NL step  0, |residual|_2 = 3.780531e+00
  NL step  1, |residual|_2 = 1.420471e-06
  NL step  2, |residual|_2 = 2.581077e-15
Writing state 35 at time 3.400000e-09 to file bih.1.e; output a total of 35 states so far
Solving for state 36, time 3.400000e-09
  NL step  0, |residual|_2 = 3.778347e+00
  NL step  1, |residual|_2 = 1.645441e-06
  NL step  2, |residual|_2 = 2.597340e-15
Writing state 36 at time 3.500000e-09 to file bih.1.e; output a total of 36 states so far
Solving for state 37, time 3.500000e-09
  NL step  0, |residual|_2 = 3.776174e+00
  NL step  1, |residual|_2 = 1.541545e-06
  NL step  2, |residual|_2 = 2.651339e-15
Writing state 37 at time 3.600000e-09 to file bih.1.e; output a total of 37 states so far
Solving for state 38, time 3.600000e-09
  NL step  0, |residual|_2 = 3.774015e+00
  NL step  1, |residual|_2 = 1.724389e-06
  NL step  2, |residual|_2 = 2.390707e-15
Writing state 38 at time 3.700000e-09 to file bih.1.e; output a total of 38 states so far
Solving for state 39, time 3.700000e-09
  NL step  0, |residual|_2 = 3.771868e+00
  NL step  1, |residual|_2 = 1.627988e-06
  NL step  2, |residual|_2 = 2.736122e-15
Writing state 39 at time 3.800000e-09 to file bih.1.e; output a total of 39 states so far
Solving for state 40, time 3.800000e-09
  NL step  0, |residual|_2 = 3.769733e+00
  NL step  1, |residual|_2 = 1.774159e-06
  NL step  2, |residual|_2 = 2.447678e-15
Writing state 40 at time 3.900000e-09 to file bih.1.e; output a total of 40 states so far
Solving for state 41, time 3.900000e-09
  NL step  0, |residual|_2 = 3.767609e+00
  NL step  1, |residual|_2 = 1.684397e-06
  NL step  2, |residual|_2 = 2.498992e-15
Writing state 41 at time 4.000000e-09 to file bih.1.e; output a total of 41 states so far
Solving for state 42, time 4.000000e-09
  NL step  0, |residual|_2 = 3.765498e+00
  NL step  1, |residual|_2 = 1.798211e-06
  NL step  2, |residual|_2 = 2.626303e-15
Writing state 42 at time 4.100000e-09 to file bih.1.e; output a total of 42 states so far
Solving for state 43, time 4.100000e-09
  NL step  0, |residual|_2 = 3.763397e+00
  NL step  1, |residual|_2 = 1.715223e-06
  NL step  2, |residual|_2 = 2.469608e-15
Writing state 43 at time 4.200000e-09 to file bih.1.e; output a total of 43 states so far
Solving for state 44, time 4.200000e-09
  NL step  0, |residual|_2 = 3.761307e+00
  NL step  1, |residual|_2 = 1.800130e-06
  NL step  2, |residual|_2 = 2.456156e-15
Writing state 44 at time 4.300000e-09 to file bih.1.e; output a total of 44 states so far
Solving for state 45, time 4.300000e-09
  NL step  0, |residual|_2 = 3.759228e+00
  NL step  1, |residual|_2 = 1.724580e-06
  NL step  2, |residual|_2 = 2.484353e-15
Writing state 45 at time 4.400000e-09 to file bih.1.e; output a total of 45 states so far
Solving for state 46, time 4.400000e-09
  NL step  0, |residual|_2 = 3.757160e+00
  NL step  1, |residual|_2 = 1.783358e-06
  NL step  2, |residual|_2 = 2.379383e-15
Writing state 46 at time 4.500000e-09 to file bih.1.e; output a total of 46 states so far
Solving for state 47, time 4.500000e-09
  NL step  0, |residual|_2 = 3.755101e+00
  NL step  1, |residual|_2 = 1.716175e-06
  NL step  2, |residual|_2 = 2.426828e-15
Writing state 47 at time 4.600000e-09 to file bih.1.e; output a total of 47 states so far
Solving for state 48, time 4.600000e-09
  NL step  0, |residual|_2 = 3.753053e+00
  NL step  1, |residual|_2 = 1.751063e-06
  NL step  2, |residual|_2 = 2.560313e-15
Writing state 48 at time 4.700000e-09 to file bih.1.e; output a total of 48 states so far
Solving for state 49, time 4.700000e-09
  NL step  0, |residual|_2 = 3.751014e+00
  NL step  1, |residual|_2 = 1.693285e-06
  NL step  2, |residual|_2 = 2.477943e-15
Writing state 49 at time 4.800000e-09 to file bih.1.e; output a total of 49 states so far
Solving for state 50, time 4.800000e-09
  NL step  0, |residual|_2 = 3.748985e+00
  NL step  1, |residual|_2 = 1.706077e-06
  NL step  2, |residual|_2 = 2.455758e-15
Writing state 50 at time 4.900000e-09 to file bih.1.e; output a total of 50 states so far
Solving for state 51, time 4.900000e-09
  NL step  0, |residual|_2 = 3.746966e+00
  NL step  1, |residual|_2 = 1.658773e-06
  NL step  2, |residual|_2 = 2.465064e-15
Writing state 51 at time 5.000000e-09 to file bih.1.e; output a total of 51 states so far
Solving for state 52, time 5.000000e-09
  NL step  0, |residual|_2 = 3.744956e+00
  NL step  1, |residual|_2 = 1.650880e-06
  NL step  2, |residual|_2 = 2.398226e-15
Writing state 52 at time 5.100000e-09 to file bih.1.e; output a total of 52 states so far
Solving for state 53, time 5.100000e-09
  NL step  0, |residual|_2 = 3.742955e+00
  NL step  1, |residual|_2 = 1.615116e-06
  NL step  2, |residual|_2 = 2.594251e-15
Writing state 53 at time 5.200000e-09 to file bih.1.e; output a total of 53 states so far
Solving for state 54, time 5.200000e-09
  NL step  0, |residual|_2 = 3.740963e+00
  NL step  1, |residual|_2 = 1.587613e-06
  NL step  2, |residual|_2 = 2.435817e-15
Writing state 54 at time 5.300000e-09 to file bih.1.e; output a total of 54 states so far
Solving for state 55, time 5.300000e-09
  NL step  0, |residual|_2 = 3.738979e+00
  NL step  1, |residual|_2 = 1.564448e-06
  NL step  2, |residual|_2 = 2.446763e-15
Writing state 55 at time 5.400000e-09 to file bih.1.e; output a total of 55 states so far
Solving for state 56, time 5.400000e-09
  NL step  0, |residual|_2 = 3.737005e+00
  NL step  1, |residual|_2 = 1.518104e-06
  NL step  2, |residual|_2 = 2.566135e-15
Writing state 56 at time 5.500000e-09 to file bih.1.e; output a total of 56 states so far
Solving for state 57, time 5.500000e-09
  NL step  0, |residual|_2 = 3.735039e+00
  NL step  1, |residual|_2 = 1.508598e-06
  NL step  2, |residual|_2 = 2.470940e-15
Writing state 57 at time 5.600000e-09 to file bih.1.e; output a total of 57 states so far
Solving for state 58, time 5.600000e-09
  NL step  0, |residual|_2 = 3.733082e+00
  NL step  1, |residual|_2 = 1.443902e-06
  NL step  2, |residual|_2 = 2.521951e-15
Writing state 58 at time 5.700000e-09 to file bih.1.e; output a total of 58 states so far
Solving for state 59, time 5.700000e-09
  NL step  0, |residual|_2 = 3.731132e+00
  NL step  1, |residual|_2 = 1.449135e-06
  NL step  2, |residual|_2 = 2.381519e-15
Writing state 59 at time 5.800000e-09 to file bih.1.e; output a total of 59 states so far
Solving for state 60, time 5.800000e-09
  NL step  0, |residual|_2 = 3.729191e+00
  NL step  1, |residual|_2 = 1.366309e-06
  NL step  2, |residual|_2 = 2.411224e-15
Writing state 60 at time 5.900000e-09 to file bih.1.e; output a total of 60 states so far
Solving for state 61, time 5.900000e-09
  NL step  0, |residual|_2 = 3.727257e+00
  NL step  1, |residual|_2 = 1.387402e-06
  NL step  2, |residual|_2 = 2.695918e-15
Writing state 61 at time 6.000000e-09 to file bih.1.e; output a total of 61 states so far
Solving for state 62, time 6.000000e-09
  NL step  0, |residual|_2 = 3.725332e+00
  NL step  1, |residual|_2 = 1.286417e-06
  NL step  2, |residual|_2 = 2.556915e-15
Writing state 62 at time 6.100000e-09 to file bih.1.e; output a total of 62 states so far
Solving for state 63, time 6.100000e-09
  NL step  0, |residual|_2 = 3.723414e+00
  NL step  1, |residual|_2 = 1.324557e-06
  NL step  2, |residual|_2 = 2.478048e-15
Writing state 63 at time 6.200000e-09 to file bih.1.e; output a total of 63 states so far
Solving for state 64, time 6.200000e-09
  NL step  0, |residual|_2 = 3.721504e+00
  NL step  1, |residual|_2 = 1.205140e-06
  NL step  2, |residual|_2 = 2.471400e-15
Writing state 64 at time 6.300000e-09 to file bih.1.e; output a total of 64 states so far
Solving for state 65, time 6.300000e-09
  NL step  0, |residual|_2 = 3.719601e+00
  NL step  1, |residual|_2 = 1.261600e-06
  NL step  2, |residual|_2 = 2.480741e-15
Writing state 65 at time 6.400000e-09 to file bih.1.e; output a total of 65 states so far
Solving for state 66, time 6.400000e-09
  NL step  0, |residual|_2 = 3.717706e+00
  NL step  1, |residual|_2 = 1.123245e-06
  NL step  2, |residual|_2 = 2.413220e-15
Writing state 66 at time 6.500000e-09 to file bih.1.e; output a total of 66 states so far
Solving for state 67, time 6.500000e-09
  NL step  0, |residual|_2 = 3.715818e+00
  NL step  1, |residual|_2 = 1.199408e-06
  NL step  2, |residual|_2 = 2.465159e-15
Writing state 67 at time 6.600000e-09 to file bih.1.e; output a total of 67 states so far
Solving for state 68, time 6.600000e-09
  NL step  0, |residual|_2 = 3.713937e+00
  NL step  1, |residual|_2 = 1.041374e-06
  NL step  2, |residual|_2 = 2.547517e-15
Writing state 68 at time 6.700000e-09 to file bih.1.e; output a total of 68 states so far
Solving for state 69, time 6.700000e-09
  NL step  0, |residual|_2 = 3.712063e+00
  NL step  1, |residual|_2 = 1.138751e-06
  NL step  2, |residual|_2 = 2.356805e-15
Writing state 69 at time 6.800000e-09 to file bih.1.e; output a total of 69 states so far
Solving for state 70, time 6.800000e-09
  NL step  0, |residual|_2 = 3.710197e+00
  NL step  1, |residual|_2 = 9.600750e-07
  NL step  2, |residual|_2 = 2.643616e-15
Writing state 70 at time 6.900000e-09 to file bih.1.e; output a total of 70 states so far
Solving for state 71, time 6.900000e-09
  NL step  0, |residual|_2 = 3.708336e+00
  NL step  1, |residual|_2 = 1.080324e-06
  NL step  2, |residual|_2 = 2.533206e-15
Writing state 71 at time 7.000000e-09 to file bih.1.e; output a total of 71 states so far
Solving for state 72, time 7.000000e-09
  NL step  0, |residual|_2 = 3.706484e+00
  NL step  1, |residual|_2 = 8.798224e-07
  NL step  2, |residual|_2 = 2.386217e-15
Writing state 72 at time 7.100000e-09 to file bih.1.e; output a total of 72 states so far
Solving for state 73, time 7.100000e-09
  NL step  0, |residual|_2 = 3.704637e+00
  NL step  1, |residual|_2 = 1.024755e-06
  NL step  2, |residual|_2 = 2.584296e-15
Writing state 73 at time 7.200000e-09 to file bih.1.e; output a total of 73 states so far
Solving for state 74, time 7.200000e-09
  NL step  0, |residual|_2 = 3.702797e+00
  NL step  1, |residual|_2 = 8.010389e-07
  NL step  2, |residual|_2 = 2.538293e-15
Writing state 74 at time 7.300000e-09 to file bih.1.e; output a total of 74 states so far
Solving for state 75, time 7.300000e-09
  NL step  0, |residual|_2 = 3.700963e+00
  NL step  1, |residual|_2 = 9.726188e-07
  NL step  2, |residual|_2 = 2.462608e-15
Writing state 75 at time 7.400000e-09 to file bih.1.e; output a total of 75 states so far
Solving for state 76, time 7.400000e-09
  NL step  0, |residual|_2 = 3.699137e+00
  NL step  1, |residual|_2 = 7.241212e-07
  NL step  2, |residual|_2 = 2.587927e-15
Writing state 76 at time 7.500000e-09 to file bih.1.e; output a total of 76 states so far
Solving for state 77, time 7.500000e-09
  NL step  0, |residual|_2 = 3.697316e+00
  NL step  1, |residual|_2 = 9.244483e-07
  NL step  2, |residual|_2 = 2.348597e-15
Writing state 77 at time 7.600000e-09 to file bih.1.e; output a total of 77 states so far
Solving for state 78, time 7.600000e-09
  NL step  0, |residual|_2 = 3.695502e+00
  NL step  1, |residual|_2 = 6.494678e-07
  NL step  2, |residual|_2 = 2.486855e-15
Writing state 78 at time 7.700000e-09 to file bih.1.e; output a total of 78 states so far
Solving for state 79, time 7.700000e-09
  NL step  0, |residual|_2 = 3.693694e+00
  NL step  1, |residual|_2 = 8.807293e-07
  NL step  2, |residual|_2 = 2.498792e-15
Writing state 79 at time 7.800000e-09 to file bih.1.e; output a total of 79 states so far
Solving for state 80, time 7.800000e-09
  NL step  0, |residual|_2 = 3.691892e+00
  NL step  1, |residual|_2 = 5.775163e-07
  NL step  2, |residual|_2 = 2.668040e-15
Writing state 80 at time 7.900000e-09 to file bih.1.e; output a total of 80 states so far
Solving for state 81, time 7.900000e-09
  NL step  0, |residual|_2 = 3.690096e+00
  NL step  1, |residual|_2 = 8.418965e-07
  NL step  2, |residual|_2 = 2.555825e-15
Writing state 81 at time 8.000000e-09 to file bih.1.e; output a total of 81 states so far
Solving for state 82, time 8.000000e-09
  NL step  0, |residual|_2 = 3.688306e+00
  NL step  1, |residual|_2 = 5.087959e-07
  NL step  2, |residual|_2 = 2.585262e-15
Writing state 82 at time 8.100000e-09 to file bih.1.e; output a total of 82 states so far
Solving for state 83, time 8.100000e-09
  NL step  0, |residual|_2 = 3.686522e+00
  NL step  1, |residual|_2 = 8.083197e-07
  NL step  2, |residual|_2 = 2.383634e-15
Writing state 83 at time 8.200000e-09 to file bih.1.e; output a total of 83 states so far
Solving for state 84, time 8.200000e-09
  NL step  0, |residual|_2 = 3.684744e+00
  NL step  1, |residual|_2 = 4.440076e-07
  NL step  2, |residual|_2 = 2.496927e-15
Writing state 84 at time 8.300000e-09 to file bih.1.e; output a total of 84 states so far
Solving for state 85, time 8.300000e-09
  NL step  0, |residual|_2 = 3.682971e+00
  NL step  1, |residual|_2 = 7.802849e-07
  NL step  2, |residual|_2 = 2.492994e-15
Writing state 85 at time 8.400000e-09 to file bih.1.e; output a total of 85 states so far
Solving for state 86, time 8.400000e-09
  NL step  0, |residual|_2 = 3.681205e+00
  NL step  1, |residual|_2 = 3.841496e-07
  NL step  2, |residual|_2 = 2.448370e-15
Writing state 86 at time 8.500000e-09 to file bih.1.e; output a total of 86 states so far
Solving for state 87, time 8.500000e-09
  NL step  0, |residual|_2 = 3.679443e+00
  NL step  1, |residual|_2 = 7.579725e-07
  NL step  2, |residual|_2 = 2.356195e-15
Writing state 87 at time 8.600000e-09 to file bih.1.e; output a total of 87 states so far
Solving for state 88, time 8.600000e-09
  NL step  0, |residual|_2 = 3.677688e+00
  NL step  1, |residual|_2 = 3.307043e-07
  NL step  2, |residual|_2 = 2.527272e-15
Writing state 88 at time 8.700000e-09 to file bih.1.e; output a total of 88 states so far
Solving for state 89, time 8.700000e-09
  NL step  0, |residual|_2 = 3.675938e+00
  NL step  1, |residual|_2 = 7.414374e-07
  NL step  2, |residual|_2 = 2.411819e-15
Writing state 89 at time 8.800000e-09 to file bih.1.e; output a total of 89 states so far
Solving for state 90, time 8.800000e-09
  NL step  0, |residual|_2 = 3.674193e+00
  NL step  1, |residual|_2 = 2.858712e-07
  NL step  2, |residual|_2 = 2.505766e-15
Writing state 90 at time 8.900000e-09 to file bih.1.e; output a total of 90 states so far
Solving for state 91, time 8.900000e-09
  NL step  0, |residual|_2 = 3.672454e+00
  NL step  1, |residual|_2 = 7.305966e-07
  NL step  2, |residual|_2 = 2.588068e-15
Writing state 91 at time 9.000000e-09 to file bih.1.e; output a total of 91 states so far
Solving for state 92, time 9.000000e-09
  NL step  0, |residual|_2 = 3.670720e+00
  NL step  1, |residual|_2 = 2.526672e-07
  NL step  2, |residual|_2 = 2.370847e-15
Writing state 92 at time 9.100000e-09 to file bih.1.e; output a total of 92 states so far
Solving for state 93, time 9.100000e-09
  NL step  0, |residual|_2 = 3.668992e+00
  NL step  1, |residual|_2 = 7.252278e-07
  NL step  2, |residual|_2 = 2.463545e-15
Writing state 93 at time 9.200000e-09 to file bih.1.e; output a total of 93 states so far
Solving for state 94, time 9.200000e-09
  NL step  0, |residual|_2 = 3.667269e+00
  NL step  1, |residual|_2 = 2.343756e-07
  NL step  2, |residual|_2 = 2.661797e-15
Writing state 94 at time 9.300000e-09 to file bih.1.e; output a total of 94 states so far
Solving for state 95, time 9.300000e-09
  NL step  0, |residual|_2 = 3.665550e+00
  NL step  1, |residual|_2 = 7.249825e-07
  NL step  2, |residual|_2 = 2.552010e-15
Writing state 95 at time 9.400000e-09 to file bih.1.e; output a total of 95 states so far
Solving for state 96, time 9.400000e-09
  NL step  0, |residual|_2 = 3.663837e+00
  NL step  1, |residual|_2 = 2.328701e-07
  NL step  2, |residual|_2 = 2.491362e-15
Writing state 96 at time 9.500000e-09 to file bih.1.e; output a total of 96 states so far
Solving for state 97, time 9.500000e-09
  NL step  0, |residual|_2 = 3.662129e+00
  NL step  1, |residual|_2 = 7.294108e-07
  NL step  2, |residual|_2 = 2.517984e-15
Writing state 97 at time 9.600000e-09 to file bih.1.e; output a total of 97 states so far
Solving for state 98, time 9.600000e-09
  NL step  0, |residual|_2 = 3.660426e+00
  NL step  1, |residual|_2 = 2.469601e-07
  NL step  2, |residual|_2 = 2.485030e-15
Writing state 98 at time 9.700000e-09 to file bih.1.e; output a total of 98 states so far
Solving for state 99, time 9.700000e-09
  NL step  0, |residual|_2 = 3.658728e+00
  NL step  1, |residual|_2 = 7.379944e-07
  NL step  2, |residual|_2 = 2.536646e-15
Writing state 99 at time 9.800000e-09 to file bih.1.e; output a total of 99 states so far
Solving for state 100, time 9.800000e-09
  NL step  0, |residual|_2 = 3.657035e+00
  NL step  1, |residual|_2 = 2.729540e-07
  NL step  2, |residual|_2 = 2.508287e-15
Writing state 100 at time 9.900000e-09 to file bih.1.e; output a total of 100 states so far
Solving for state 101, time 9.900000e-09
  NL step  0, |residual|_2 = 3.655347e+00
  NL step  1, |residual|_2 = 7.501811e-07
  NL step  2, |residual|_2 = 2.443282e-15
Writing state 101 at time 1.000000e-08 to file bih.1.e; output a total of 101 states so far
Solving for state 102, time 1.000000e-08
  NL step  0, |residual|_2 = 3.653664e+00
  NL step  1, |residual|_2 = 3.067491e-07
  NL step  2, |residual|_2 = 2.447045e-15
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./miscellaneous_ex7-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:22:03 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.332e+00      1.00012   2.332e+00
Objects:              4.716e+03      1.00000   4.716e+03
Flops:                1.092e+07      1.01236   1.084e+07  6.503e+07
Flops/sec:            4.681e+06      1.01241   4.647e+06  2.788e+07
MPI Messages:         1.163e+04      1.03134   1.151e+04  6.908e+04
MPI Message Lengths:  1.034e+07      1.00477   8.955e+02  6.187e+07
MPI Reductions:       1.585e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.3321e+00 100.0%  6.5026e+07 100.0%  6.908e+04 100.0%  8.955e+02      100.0%  1.493e+04  94.2% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecDot               198 1.0 1.1677e-02 2.0 1.36e+05 1.0 0.0e+00 0.0e+00 2.0e+02  0  1  0  0  1   0  1  0  0  1    69
VecMDot              578 1.0 8.5711e-02 2.9 1.04e+06 1.0 0.0e+00 0.0e+00 5.8e+02  2 10  0  0  4   2 10  0  0  4    72
VecNorm             1172 1.0 7.1512e-02 1.5 8.06e+05 1.0 0.0e+00 0.0e+00 1.2e+03  2  7  0  0  7   2  7  0  0  8    67
VecScale             776 1.0 5.1022e-04 1.4 2.67e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0  3118
VecCopy             1889 1.0 1.2975e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              1478 1.0 5.4383e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              198 1.0 3.0231e-04 1.5 1.36e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  2685
VecWAXPY             198 1.0 1.8263e-04 1.4 6.81e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  2223
VecMAXPY             776 1.0 8.0895e-04 1.3 1.44e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0 10613
VecAssemblyBegin    3290 1.0 9.0213e-01 1.8 0.00e+00 0.0 3.6e+03 1.2e+01 8.4e+03 31  0  5  0 53  31  0  5  0 56     0
VecAssemblyEnd      3290 1.0 2.3174e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1975 1.0 6.8219e-03 1.4 0.00e+00 0.0 3.2e+04 1.3e+03 0.0e+00  0  0 46 67  0   0  0 46 67  0     0
VecScatterEnd       1975 1.0 5.6910e-0112.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 10  0  0  0  0  10  0  0  0  0     0
VecReduceArith       202 1.0 6.9292e-0313.1 1.39e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   119
VecReduceComm        101 1.0 1.0018e-02 6.4 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  1   0  0  0  0  1     0
VecNormalize         776 1.0 4.3432e-02 2.8 6.65e+05 1.0 0.0e+00 0.0e+00 5.8e+02  1  6  0  0  4   1  6  0  0  4    91
MatMult              776 1.0 7.4355e-02 8.7 2.94e+06 1.0 9.3e+03 1.7e+01 0.0e+00  1 27 13  0  0   1 27 13  0  0   236
MatSolve             776 1.0 5.4607e-03 1.4 2.92e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0 27  0  0  0   0 27  0  0  0  3191
MatLUFactorNum       198 1.0 1.1976e-02 1.3 9.48e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  9  0  0  0   0  9  0  0  0   471
MatILUFactorSym      101 1.0 7.3094e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin     396 1.0 2.2760e-01 1.7 0.00e+00 0.0 3.6e+03 4.8e+01 7.9e+02  8  0  5  0  5   8  0  5  0  5     0
MatAssemblyEnd       396 1.0 2.7126e-02 1.3 0.00e+00 0.0 2.4e+01 6.3e+00 4.0e+02  1  0  0  0  3   1  0  0  0  3     0
MatGetRowIJ          101 1.0 3.0994e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       101 1.0 2.2457e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 4.0e+02  0  0  0  0  3   0  0  0  0  3     0
MatZeroEntries       200 1.0 3.4428e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
SNESSolve            101 1.0 2.0841e+00 1.0 1.09e+07 1.0 6.7e+04 9.3e+02 1.4e+04 89100 97100 90  89100 97100 95    31
SNESLineSearch       198 1.0 7.4729e-01 1.0 1.30e+06 1.0 2.5e+04 9.9e+02 5.3e+03 32 12 36 40 34  32 12 36 40 36    10
SNESFunctionEval     299 1.0 1.0620e+00 1.0 0.00e+00 0.0 3.4e+04 1.1e+03 6.9e+03 45  0 49 60 43  45  0 49 60 46     0
SNESJacobianEval     198 1.0 8.2387e-01 1.0 0.00e+00 0.0 2.4e+04 1.0e+03 4.6e+03 34  0 34 40 29  34  0 34 40 31     0
KSPGMRESOrthog       578 1.0 8.6764e-02 2.9 2.08e+06 1.0 0.0e+00 0.0e+00 5.8e+02  2 19  0  0  4   2 19  0  0  4   143
KSPSetup             396 1.0 1.9369e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             198 1.0 1.6107e-01 1.2 9.35e+06 1.0 6.9e+03 1.7e+01 1.7e+03  7 86 10  0 10   7 86 10  0 11   346
PCSetUp              396 1.0 2.9628e-02 1.3 9.48e+05 1.0 0.0e+00 0.0e+00 5.0e+02  1  9  0  0  3   1  9  0  0  3   190
PCSetUpOnBlocks      198 1.0 2.4673e-02 1.4 9.48e+05 1.0 0.0e+00 0.0e+00 5.0e+02  1  9  0  0  3   1  9  0  0  3   229
PCApply              776 1.0 1.0792e-02 1.4 2.92e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0 27  0  0  0   0 27  0  0  0  1614
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  2591           2591     16764672     0
         Vec Scatter   502            502       435736     0
           Index Set  1011           1011      1078432     0
   IS L to G Mapping     3              3         5364     0
              Matrix   104            104      2745740     0
                SNES   101            101       104232     0
       Krylov Solver   202            202      1906880     0
      Preconditioner   202            202       142208     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.282e-05
Average time for zero size MPI_Send(): 4.21604e-05
#PETSc Option Table entries:
--verbose dim=1
-ksp_right_pc
-log_summary
-pc_type bjacobi
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 102 states so far
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 103 states so far

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:22:03 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.52706, Active time=2.23806                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000075    0.0001      0.000078    0.00     0.00     |
|   build_constraint_matrix()        84490     0.0172      0.000000    0.0172      0.000000    0.77     0.77     |
|   build_sparsity()                 1         0.0009      0.000887    0.0012      0.001228    0.04     0.05     |
|   constrain_elem_matrix()          33660     0.0062      0.000000    0.0062      0.000000    0.28     0.28     |
|   constrain_elem_vector()          50830     0.0154      0.000000    0.0154      0.000000    0.69     0.69     |
|   create_dof_constraints()         1         0.0011      0.001080    0.0029      0.002922    0.05     0.13     |
|   distribute_dofs()                1         0.0003      0.000277    0.0024      0.002439    0.01     0.11     |
|   dof_indices()                    103200    0.0359      0.000000    0.0359      0.000000    1.61     1.61     |
|   enforce_constraints_exactly()    498       0.3038      0.000610    0.3038      0.000610    13.57    13.57    |
|   prepare_send_list()              1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                         1         0.0006      0.000554    0.0006      0.000554    0.02     0.02     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          103       0.0349      0.000339    0.1467      0.001425    1.56     6.56     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               103       0.0049      0.000047    0.0049      0.000047    0.22     0.22     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        501       0.0014      0.000003    0.0014      0.000003    0.06     0.06     |
|   init_shape_functions()           501       0.0045      0.000009    0.0045      0.000009    0.20     0.20     |
|   inverse_map()                    9         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             84494     0.0312      0.000000    0.0312      0.000000    1.39     1.39     |
|   compute_face_map()               2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_face_shape_functions()      2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_reference_to_physical_map() 501       0.0009      0.000002    0.0009      0.000002    0.04     0.04     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0006      0.000576    0.0008      0.000757    0.03     0.03     |
|   renumber_nodes_and_elem()        2         0.0000      0.000023    0.0000      0.000023    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0118      0.005903    0.0118      0.005903    0.53     0.53     |
|   find_global_indices()            2         0.0004      0.000204    0.0166      0.008299    0.02     0.74     |
|   parallel_sort()                  2         0.0015      0.000758    0.0035      0.001728    0.07     0.15     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         103       0.0007      0.000007    0.1524      0.001479    0.03     6.81     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0003      0.000274    0.0003      0.000274    0.01     0.01     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0013      0.001268    0.0056      0.005633    0.06     0.25     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      7         0.0013      0.000183    0.0013      0.000183    0.06     0.06     |
|   max(scalar)                      2         0.0027      0.001360    0.0027      0.001360    0.12     0.12     |
|   max(vector)                      4         0.0001      0.000027    0.0001      0.000027    0.00     0.00     |
|   min(vector)                      4         0.0005      0.000133    0.0005      0.000133    0.02     0.02     |
|   probe()                          50        0.0026      0.000052    0.0026      0.000052    0.12     0.12     |
|   receive()                        50        0.0001      0.000002    0.0027      0.000053    0.00     0.12     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()                   54        0.0001      0.000002    0.0029      0.000053    0.00     0.13     |
|   sum()                            813       0.1340      0.000165    0.1340      0.000165    5.99     5.99     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000108    0.0017      0.001721    0.00     0.08     |
|   set_parent_processor_ids()       1         0.0001      0.000062    0.0001      0.000062    0.00     0.00     |
|                                                                                                                |
| PetscNonlinearSolver                                                                                           |
|   jacobian()                       198       0.5995      0.003028    0.7906      0.003993    26.79    35.32    |
|   residual()                       299       0.7834      0.002620    1.0526      0.003521    35.00    47.03    |
|   solve()                          101       0.2347      0.002324    2.0780      0.020574    10.49    92.85    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  1         0.0014      0.001375    0.0017      0.001678    0.06     0.07     |
|   operator()                       4         0.0000      0.000008    0.0000      0.000010    0.00     0.00     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 1         0.0010      0.001000    0.0021      0.002065    0.04     0.09     |
|   solve()                          101       0.0007      0.000007    2.0787      0.020581    0.03     92.88    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            360805    2.2381                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./miscellaneous_ex7-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
</pre>
</div>
<?php make_footer() ?>
</body>
</html>
<?php if (0) { ?>
\#Local Variables:
\#mode: html
\#End:
<?php } ?>
