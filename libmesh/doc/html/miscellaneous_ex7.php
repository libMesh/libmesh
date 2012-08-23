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
Compiling C++ (in optimized mode) biharmonic.C...
Compiling C++ (in optimized mode) biharmonic_jr.C...
Compiling C++ (in optimized mode) miscellaneous_ex7.C...
Linking miscellaneous_ex7-opt...
***************************************************************
* Running Example  ./miscellaneous_ex7-opt
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
  NL step  1, |residual|_2 = 3.700105e-14
Writing state 2 at time 1.000000e-10 to file bih.1.e; output a total of 2 states so far
Solving for state 3, time 1.000000e-10
  NL step  0, |residual|_2 = 3.861870e+00
  NL step  1, |residual|_2 = 9.419949e-15
Writing state 3 at time 2.000000e-10 to file bih.1.e; output a total of 3 states so far
Solving for state 4, time 2.000000e-10
  NL step  0, |residual|_2 = 3.858634e+00
  NL step  1, |residual|_2 = 2.484184e-14
Writing state 4 at time 3.000000e-10 to file bih.1.e; output a total of 4 states so far
Solving for state 5, time 3.000000e-10
  NL step  0, |residual|_2 = 3.855483e+00
  NL step  1, |residual|_2 = 2.581990e-15
Writing state 5 at time 4.000000e-10 to file bih.1.e; output a total of 5 states so far
Solving for state 6, time 4.000000e-10
  NL step  0, |residual|_2 = 3.852464e+00
  NL step  1, |residual|_2 = 1.013857e-14
Writing state 6 at time 5.000000e-10 to file bih.1.e; output a total of 6 states so far
Solving for state 7, time 5.000000e-10
  NL step  0, |residual|_2 = 3.849501e+00
  NL step  1, |residual|_2 = 1.458791e-14
Writing state 7 at time 6.000000e-10 to file bih.1.e; output a total of 7 states so far
Solving for state 8, time 6.000000e-10
  NL step  0, |residual|_2 = 3.846621e+00
  NL step  1, |residual|_2 = 1.080643e-14
Writing state 8 at time 7.000000e-10 to file bih.1.e; output a total of 8 states so far
Solving for state 9, time 7.000000e-10
  NL step  0, |residual|_2 = 3.843784e+00
  NL step  1, |residual|_2 = 1.551695e-14
Writing state 9 at time 8.000000e-10 to file bih.1.e; output a total of 9 states so far
Solving for state 10, time 8.000000e-10
  NL step  0, |residual|_2 = 3.841009e+00
  NL step  1, |residual|_2 = 1.690656e-14
Writing state 10 at time 9.000000e-10 to file bih.1.e; output a total of 10 states so far
Solving for state 11, time 9.000000e-10
  NL step  0, |residual|_2 = 3.838268e+00
  NL step  1, |residual|_2 = 2.450236e-15
Writing state 11 at time 1.000000e-09 to file bih.1.e; output a total of 11 states so far
Solving for state 12, time 1.000000e-09
  NL step  0, |residual|_2 = 3.835577e+00
  NL step  1, |residual|_2 = 6.586079e-15
Writing state 12 at time 1.100000e-09 to file bih.1.e; output a total of 12 states so far
Solving for state 13, time 1.100000e-09
  NL step  0, |residual|_2 = 3.832915e+00
  NL step  1, |residual|_2 = 5.896559e-15
Writing state 13 at time 1.200000e-09 to file bih.1.e; output a total of 13 states so far
Solving for state 14, time 1.200000e-09
  NL step  0, |residual|_2 = 3.830295e+00
  NL step  1, |residual|_2 = 3.431924e-15
Writing state 14 at time 1.300000e-09 to file bih.1.e; output a total of 14 states so far
Solving for state 15, time 1.300000e-09
  NL step  0, |residual|_2 = 3.827700e+00
  NL step  1, |residual|_2 = 5.412929e-15
Writing state 15 at time 1.400000e-09 to file bih.1.e; output a total of 15 states so far
Solving for state 16, time 1.400000e-09
  NL step  0, |residual|_2 = 3.825140e+00
  NL step  1, |residual|_2 = 4.992858e-15
Writing state 16 at time 1.500000e-09 to file bih.1.e; output a total of 16 states so far
Solving for state 17, time 1.500000e-09
  NL step  0, |residual|_2 = 3.822604e+00
  NL step  1, |residual|_2 = 9.898635e-15
Writing state 17 at time 1.600000e-09 to file bih.1.e; output a total of 17 states so far
Solving for state 18, time 1.600000e-09
  NL step  0, |residual|_2 = 3.820098e+00
  NL step  1, |residual|_2 = 4.071119e-15
Writing state 18 at time 1.700000e-09 to file bih.1.e; output a total of 18 states so far
Solving for state 19, time 1.700000e-09
  NL step  0, |residual|_2 = 3.817613e+00
  NL step  1, |residual|_2 = 2.620636e-15
Writing state 19 at time 1.800000e-09 to file bih.1.e; output a total of 19 states so far
Solving for state 20, time 1.800000e-09
  NL step  0, |residual|_2 = 3.815156e+00
  NL step  1, |residual|_2 = 5.960056e-15
Writing state 20 at time 1.900000e-09 to file bih.1.e; output a total of 20 states so far
Solving for state 21, time 1.900000e-09
  NL step  0, |residual|_2 = 3.812718e+00
  NL step  1, |residual|_2 = 1.107967e-14
Writing state 21 at time 2.000000e-09 to file bih.1.e; output a total of 21 states so far
Solving for state 22, time 2.000000e-09
  NL step  0, |residual|_2 = 3.810305e+00
  NL step  1, |residual|_2 = 1.535662e-14
Writing state 22 at time 2.100000e-09 to file bih.1.e; output a total of 22 states so far
Solving for state 23, time 2.100000e-09
  NL step  0, |residual|_2 = 3.807909e+00
  NL step  1, |residual|_2 = 1.176115e-14
Writing state 23 at time 2.200000e-09 to file bih.1.e; output a total of 23 states so far
Solving for state 24, time 2.200000e-09
  NL step  0, |residual|_2 = 3.805536e+00
  NL step  1, |residual|_2 = 1.691962e-14
Writing state 24 at time 2.300000e-09 to file bih.1.e; output a total of 24 states so far
Solving for state 25, time 2.300000e-09
  NL step  0, |residual|_2 = 3.803179e+00
  NL step  1, |residual|_2 = 3.363071e-15
Writing state 25 at time 2.400000e-09 to file bih.1.e; output a total of 25 states so far
Solving for state 26, time 2.400000e-09
  NL step  0, |residual|_2 = 3.800842e+00
  NL step  1, |residual|_2 = 1.200336e-14
Writing state 26 at time 2.500000e-09 to file bih.1.e; output a total of 26 states so far
Solving for state 27, time 2.500000e-09
  NL step  0, |residual|_2 = 3.798522e+00
  NL step  1, |residual|_2 = 8.651098e-15
Writing state 27 at time 2.600000e-09 to file bih.1.e; output a total of 27 states so far
Solving for state 28, time 2.600000e-09
  NL step  0, |residual|_2 = 3.796220e+00
  NL step  1, |residual|_2 = 6.598501e-15
Writing state 28 at time 2.700000e-09 to file bih.1.e; output a total of 28 states so far
Solving for state 29, time 2.700000e-09
  NL step  0, |residual|_2 = 3.793932e+00
  NL step  1, |residual|_2 = 6.938762e-15
Writing state 29 at time 2.800000e-09 to file bih.1.e; output a total of 29 states so far
Solving for state 30, time 2.800000e-09
  NL step  0, |residual|_2 = 3.791663e+00
  NL step  1, |residual|_2 = 7.534192e-15
Writing state 30 at time 2.900000e-09 to file bih.1.e; output a total of 30 states so far
Solving for state 31, time 2.900000e-09
  NL step  0, |residual|_2 = 3.789407e+00
  NL step  1, |residual|_2 = 6.826261e-15
Writing state 31 at time 3.000000e-09 to file bih.1.e; output a total of 31 states so far
Solving for state 32, time 3.000000e-09
  NL step  0, |residual|_2 = 3.787167e+00
  NL step  1, |residual|_2 = 1.366032e-14
Writing state 32 at time 3.100000e-09 to file bih.1.e; output a total of 32 states so far
Solving for state 33, time 3.100000e-09
  NL step  0, |residual|_2 = 3.784941e+00
  NL step  1, |residual|_2 = 7.775207e-15
Writing state 33 at time 3.200000e-09 to file bih.1.e; output a total of 33 states so far
Solving for state 34, time 3.200000e-09
  NL step  0, |residual|_2 = 3.782730e+00
  NL step  1, |residual|_2 = 5.915989e-15
Writing state 34 at time 3.300000e-09 to file bih.1.e; output a total of 34 states so far
Solving for state 35, time 3.300000e-09
  NL step  0, |residual|_2 = 3.780531e+00
  NL step  1, |residual|_2 = 1.583271e-14
Writing state 35 at time 3.400000e-09 to file bih.1.e; output a total of 35 states so far
Solving for state 36, time 3.400000e-09
  NL step  0, |residual|_2 = 3.778347e+00
  NL step  1, |residual|_2 = 2.867504e-15
Writing state 36 at time 3.500000e-09 to file bih.1.e; output a total of 36 states so far
Solving for state 37, time 3.500000e-09
  NL step  0, |residual|_2 = 3.776174e+00
  NL step  1, |residual|_2 = 6.916565e-15
Writing state 37 at time 3.600000e-09 to file bih.1.e; output a total of 37 states so far
Solving for state 38, time 3.600000e-09
  NL step  0, |residual|_2 = 3.774015e+00
  NL step  1, |residual|_2 = 8.571053e-15
Writing state 38 at time 3.700000e-09 to file bih.1.e; output a total of 38 states so far
Solving for state 39, time 3.700000e-09
  NL step  0, |residual|_2 = 3.771868e+00
  NL step  1, |residual|_2 = 2.199033e-14
Writing state 39 at time 3.800000e-09 to file bih.1.e; output a total of 39 states so far
Solving for state 40, time 3.800000e-09
  NL step  0, |residual|_2 = 3.769733e+00
  NL step  1, |residual|_2 = 6.211343e-14
Writing state 40 at time 3.900000e-09 to file bih.1.e; output a total of 40 states so far
Solving for state 41, time 3.900000e-09
  NL step  0, |residual|_2 = 3.767609e+00
  NL step  1, |residual|_2 = 1.775819e-13
Writing state 41 at time 4.000000e-09 to file bih.1.e; output a total of 41 states so far
Solving for state 42, time 4.000000e-09
  NL step  0, |residual|_2 = 3.765498e+00
  NL step  1, |residual|_2 = 5.111619e-13
Writing state 42 at time 4.100000e-09 to file bih.1.e; output a total of 42 states so far
Solving for state 43, time 4.100000e-09
  NL step  0, |residual|_2 = 3.763397e+00
  NL step  1, |residual|_2 = 1.471075e-12
Writing state 43 at time 4.200000e-09 to file bih.1.e; output a total of 43 states so far
Solving for state 44, time 4.200000e-09
  NL step  0, |residual|_2 = 3.761307e+00
  NL step  1, |residual|_2 = 4.234083e-12
Writing state 44 at time 4.300000e-09 to file bih.1.e; output a total of 44 states so far
Solving for state 45, time 4.300000e-09
  NL step  0, |residual|_2 = 3.759228e+00
  NL step  1, |residual|_2 = 1.218675e-11
Writing state 45 at time 4.400000e-09 to file bih.1.e; output a total of 45 states so far
Solving for state 46, time 4.400000e-09
  NL step  0, |residual|_2 = 3.757160e+00
  NL step  1, |residual|_2 = 3.507717e-11
Writing state 46 at time 4.500000e-09 to file bih.1.e; output a total of 46 states so far
Solving for state 47, time 4.500000e-09
  NL step  0, |residual|_2 = 3.755101e+00
  NL step  1, |residual|_2 = 1.009625e-10
Writing state 47 at time 4.600000e-09 to file bih.1.e; output a total of 47 states so far
Solving for state 48, time 4.600000e-09
  NL step  0, |residual|_2 = 3.753053e+00
  NL step  1, |residual|_2 = 2.906006e-10
Writing state 48 at time 4.700000e-09 to file bih.1.e; output a total of 48 states so far
Solving for state 49, time 4.700000e-09
  NL step  0, |residual|_2 = 3.751014e+00
  NL step  1, |residual|_2 = 8.364360e-10
Writing state 49 at time 4.800000e-09 to file bih.1.e; output a total of 49 states so far
Solving for state 50, time 4.800000e-09
  NL step  0, |residual|_2 = 3.748985e+00
  NL step  1, |residual|_2 = 2.407514e-09
Writing state 50 at time 4.900000e-09 to file bih.1.e; output a total of 50 states so far
Solving for state 51, time 4.900000e-09
  NL step  0, |residual|_2 = 3.746966e+00
  NL step  1, |residual|_2 = 6.929549e-09
Writing state 51 at time 5.000000e-09 to file bih.1.e; output a total of 51 states so far
Solving for state 52, time 5.000000e-09
  NL step  0, |residual|_2 = 3.744956e+00
  NL step  1, |residual|_2 = 1.994532e-08
Writing state 52 at time 5.100000e-09 to file bih.1.e; output a total of 52 states so far
Solving for state 53, time 5.100000e-09
  NL step  0, |residual|_2 = 3.742955e+00
  NL step  1, |residual|_2 = 5.740863e-08
  NL step  2, |residual|_2 = 2.384681e-15
Writing state 53 at time 5.200000e-09 to file bih.1.e; output a total of 53 states so far
Solving for state 54, time 5.200000e-09
  NL step  0, |residual|_2 = 3.740963e+00
  NL step  1, |residual|_2 = 5.660993e-08
  NL step  2, |residual|_2 = 2.419184e-15
Writing state 54 at time 5.300000e-09 to file bih.1.e; output a total of 54 states so far
Solving for state 55, time 5.300000e-09
  NL step  0, |residual|_2 = 3.738979e+00
  NL step  1, |residual|_2 = 5.606989e-08
  NL step  2, |residual|_2 = 2.469188e-15
Writing state 55 at time 5.400000e-09 to file bih.1.e; output a total of 55 states so far
Solving for state 56, time 5.400000e-09
  NL step  0, |residual|_2 = 3.737005e+00
  NL step  1, |residual|_2 = 5.559189e-08
  NL step  2, |residual|_2 = 2.442416e-15
Writing state 56 at time 5.500000e-09 to file bih.1.e; output a total of 56 states so far
Solving for state 57, time 5.500000e-09
  NL step  0, |residual|_2 = 3.735039e+00
  NL step  1, |residual|_2 = 5.518270e-08
  NL step  2, |residual|_2 = 2.431436e-15
Writing state 57 at time 5.600000e-09 to file bih.1.e; output a total of 57 states so far
Solving for state 58, time 5.600000e-09
  NL step  0, |residual|_2 = 3.733082e+00
  NL step  1, |residual|_2 = 5.480558e-08
  NL step  2, |residual|_2 = 2.400013e-15
Writing state 58 at time 5.700000e-09 to file bih.1.e; output a total of 58 states so far
Solving for state 59, time 5.700000e-09
  NL step  0, |residual|_2 = 3.731132e+00
  NL step  1, |residual|_2 = 5.446405e-08
  NL step  2, |residual|_2 = 2.488918e-15
Writing state 59 at time 5.800000e-09 to file bih.1.e; output a total of 59 states so far
Solving for state 60, time 5.800000e-09
  NL step  0, |residual|_2 = 3.729191e+00
  NL step  1, |residual|_2 = 5.414304e-08
  NL step  2, |residual|_2 = 2.415734e-15
Writing state 60 at time 5.900000e-09 to file bih.1.e; output a total of 60 states so far
Solving for state 61, time 5.900000e-09
  NL step  0, |residual|_2 = 3.727257e+00
  NL step  1, |residual|_2 = 5.384478e-08
  NL step  2, |residual|_2 = 2.465045e-15
Writing state 61 at time 6.000000e-09 to file bih.1.e; output a total of 61 states so far
Solving for state 62, time 6.000000e-09
  NL step  0, |residual|_2 = 3.725332e+00
  NL step  1, |residual|_2 = 5.356114e-08
  NL step  2, |residual|_2 = 2.404997e-15
Writing state 62 at time 6.100000e-09 to file bih.1.e; output a total of 62 states so far
Solving for state 63, time 6.100000e-09
  NL step  0, |residual|_2 = 3.723414e+00
  NL step  1, |residual|_2 = 5.329367e-08
  NL step  2, |residual|_2 = 2.347143e-15
Writing state 63 at time 6.200000e-09 to file bih.1.e; output a total of 63 states so far
Solving for state 64, time 6.200000e-09
  NL step  0, |residual|_2 = 3.721504e+00
  NL step  1, |residual|_2 = 5.303732e-08
  NL step  2, |residual|_2 = 2.466288e-15
Writing state 64 at time 6.300000e-09 to file bih.1.e; output a total of 64 states so far
Solving for state 65, time 6.300000e-09
  NL step  0, |residual|_2 = 3.719601e+00
  NL step  1, |residual|_2 = 5.279325e-08
  NL step  2, |residual|_2 = 2.396288e-15
Writing state 65 at time 6.400000e-09 to file bih.1.e; output a total of 65 states so far
Solving for state 66, time 6.400000e-09
  NL step  0, |residual|_2 = 3.717706e+00
  NL step  1, |residual|_2 = 5.255803e-08
  NL step  2, |residual|_2 = 2.367857e-15
Writing state 66 at time 6.500000e-09 to file bih.1.e; output a total of 66 states so far
Solving for state 67, time 6.500000e-09
  NL step  0, |residual|_2 = 3.715818e+00
  NL step  1, |residual|_2 = 5.233254e-08
  NL step  2, |residual|_2 = 2.396083e-15
Writing state 67 at time 6.600000e-09 to file bih.1.e; output a total of 67 states so far
Solving for state 68, time 6.600000e-09
  NL step  0, |residual|_2 = 3.713937e+00
  NL step  1, |residual|_2 = 5.211431e-08
  NL step  2, |residual|_2 = 2.331777e-15
Writing state 68 at time 6.700000e-09 to file bih.1.e; output a total of 68 states so far
Solving for state 69, time 6.700000e-09
  NL step  0, |residual|_2 = 3.712063e+00
  NL step  1, |residual|_2 = 5.190405e-08
  NL step  2, |residual|_2 = 2.171709e-15
Writing state 69 at time 6.800000e-09 to file bih.1.e; output a total of 69 states so far
Solving for state 70, time 6.800000e-09
  NL step  0, |residual|_2 = 3.710197e+00
  NL step  1, |residual|_2 = 5.169989e-08
  NL step  2, |residual|_2 = 2.386625e-15
Writing state 70 at time 6.900000e-09 to file bih.1.e; output a total of 70 states so far
Solving for state 71, time 6.900000e-09
  NL step  0, |residual|_2 = 3.708336e+00
  NL step  1, |residual|_2 = 5.150241e-08
  NL step  2, |residual|_2 = 2.390530e-15
Writing state 71 at time 7.000000e-09 to file bih.1.e; output a total of 71 states so far
Solving for state 72, time 7.000000e-09
  NL step  0, |residual|_2 = 3.706484e+00
  NL step  1, |residual|_2 = 5.131016e-08
  NL step  2, |residual|_2 = 2.452716e-15
Writing state 72 at time 7.100000e-09 to file bih.1.e; output a total of 72 states so far
Solving for state 73, time 7.100000e-09
  NL step  0, |residual|_2 = 3.704637e+00
  NL step  1, |residual|_2 = 5.112361e-08
  NL step  2, |residual|_2 = 2.443674e-15
Writing state 73 at time 7.200000e-09 to file bih.1.e; output a total of 73 states so far
Solving for state 74, time 7.200000e-09
  NL step  0, |residual|_2 = 3.702797e+00
  NL step  1, |residual|_2 = 5.094161e-08
  NL step  2, |residual|_2 = 2.504119e-15
Writing state 74 at time 7.300000e-09 to file bih.1.e; output a total of 74 states so far
Solving for state 75, time 7.300000e-09
  NL step  0, |residual|_2 = 3.700963e+00
  NL step  1, |residual|_2 = 5.076456e-08
  NL step  2, |residual|_2 = 2.483832e-15
Writing state 75 at time 7.400000e-09 to file bih.1.e; output a total of 75 states so far
Solving for state 76, time 7.400000e-09
  NL step  0, |residual|_2 = 3.699137e+00
  NL step  1, |residual|_2 = 5.059150e-08
  NL step  2, |residual|_2 = 2.485467e-15
Writing state 76 at time 7.500000e-09 to file bih.1.e; output a total of 76 states so far
Solving for state 77, time 7.500000e-09
  NL step  0, |residual|_2 = 3.697316e+00
  NL step  1, |residual|_2 = 5.042279e-08
  NL step  2, |residual|_2 = 2.485405e-15
Writing state 77 at time 7.600000e-09 to file bih.1.e; output a total of 77 states so far
Solving for state 78, time 7.600000e-09
  NL step  0, |residual|_2 = 3.695502e+00
  NL step  1, |residual|_2 = 5.025763e-08
  NL step  2, |residual|_2 = 2.462417e-15
Writing state 78 at time 7.700000e-09 to file bih.1.e; output a total of 78 states so far
Solving for state 79, time 7.700000e-09
  NL step  0, |residual|_2 = 3.693694e+00
  NL step  1, |residual|_2 = 5.009634e-08
  NL step  2, |residual|_2 = 2.330280e-15
Writing state 79 at time 7.800000e-09 to file bih.1.e; output a total of 79 states so far
Solving for state 80, time 7.800000e-09
  NL step  0, |residual|_2 = 3.691892e+00
  NL step  1, |residual|_2 = 4.993823e-08
  NL step  2, |residual|_2 = 2.557443e-15
Writing state 80 at time 7.900000e-09 to file bih.1.e; output a total of 80 states so far
Solving for state 81, time 7.900000e-09
  NL step  0, |residual|_2 = 3.690096e+00
  NL step  1, |residual|_2 = 4.978357e-08
  NL step  2, |residual|_2 = 2.510954e-15
Writing state 81 at time 8.000000e-09 to file bih.1.e; output a total of 81 states so far
Solving for state 82, time 8.000000e-09
  NL step  0, |residual|_2 = 3.688306e+00
  NL step  1, |residual|_2 = 4.963179e-08
  NL step  2, |residual|_2 = 2.532870e-15
Writing state 82 at time 8.100000e-09 to file bih.1.e; output a total of 82 states so far
Solving for state 83, time 8.100000e-09
  NL step  0, |residual|_2 = 3.686522e+00
  NL step  1, |residual|_2 = 4.948314e-08
  NL step  2, |residual|_2 = 2.409232e-15
Writing state 83 at time 8.200000e-09 to file bih.1.e; output a total of 83 states so far
Solving for state 84, time 8.200000e-09
  NL step  0, |residual|_2 = 3.684744e+00
  NL step  1, |residual|_2 = 4.933709e-08
  NL step  2, |residual|_2 = 2.373077e-15
Writing state 84 at time 8.300000e-09 to file bih.1.e; output a total of 84 states so far
Solving for state 85, time 8.300000e-09
  NL step  0, |residual|_2 = 3.682971e+00
  NL step  1, |residual|_2 = 4.919388e-08
  NL step  2, |residual|_2 = 2.567430e-15
Writing state 85 at time 8.400000e-09 to file bih.1.e; output a total of 85 states so far
Solving for state 86, time 8.400000e-09
  NL step  0, |residual|_2 = 3.681205e+00
  NL step  1, |residual|_2 = 4.905306e-08
  NL step  2, |residual|_2 = 2.455994e-15
Writing state 86 at time 8.500000e-09 to file bih.1.e; output a total of 86 states so far
Solving for state 87, time 8.500000e-09
  NL step  0, |residual|_2 = 3.679443e+00
  NL step  1, |residual|_2 = 4.891484e-08
  NL step  2, |residual|_2 = 2.374471e-15
Writing state 87 at time 8.600000e-09 to file bih.1.e; output a total of 87 states so far
Solving for state 88, time 8.600000e-09
  NL step  0, |residual|_2 = 3.677688e+00
  NL step  1, |residual|_2 = 4.877880e-08
  NL step  2, |residual|_2 = 2.483941e-15
Writing state 88 at time 8.700000e-09 to file bih.1.e; output a total of 88 states so far
Solving for state 89, time 8.700000e-09
  NL step  0, |residual|_2 = 3.675938e+00
  NL step  1, |residual|_2 = 4.864515e-08
  NL step  2, |residual|_2 = 2.398725e-15
Writing state 89 at time 8.800000e-09 to file bih.1.e; output a total of 89 states so far
Solving for state 90, time 8.800000e-09
  NL step  0, |residual|_2 = 3.674193e+00
  NL step  1, |residual|_2 = 4.851352e-08
  NL step  2, |residual|_2 = 2.530767e-15
Writing state 90 at time 8.900000e-09 to file bih.1.e; output a total of 90 states so far
Solving for state 91, time 8.900000e-09
  NL step  0, |residual|_2 = 3.672454e+00
  NL step  1, |residual|_2 = 4.838409e-08
  NL step  2, |residual|_2 = 2.615689e-15
Writing state 91 at time 9.000000e-09 to file bih.1.e; output a total of 91 states so far
Solving for state 92, time 9.000000e-09
  NL step  0, |residual|_2 = 3.670720e+00
  NL step  1, |residual|_2 = 4.825654e-08
  NL step  2, |residual|_2 = 2.568743e-15
Writing state 92 at time 9.100000e-09 to file bih.1.e; output a total of 92 states so far
Solving for state 93, time 9.100000e-09
  NL step  0, |residual|_2 = 3.668992e+00
  NL step  1, |residual|_2 = 4.813102e-08
  NL step  2, |residual|_2 = 2.470177e-15
Writing state 93 at time 9.200000e-09 to file bih.1.e; output a total of 93 states so far
Solving for state 94, time 9.200000e-09
  NL step  0, |residual|_2 = 3.667269e+00
  NL step  1, |residual|_2 = 4.800725e-08
  NL step  2, |residual|_2 = 2.370272e-15
Writing state 94 at time 9.300000e-09 to file bih.1.e; output a total of 94 states so far
Solving for state 95, time 9.300000e-09
  NL step  0, |residual|_2 = 3.665550e+00
  NL step  1, |residual|_2 = 4.788536e-08
  NL step  2, |residual|_2 = 2.495545e-15
Writing state 95 at time 9.400000e-09 to file bih.1.e; output a total of 95 states so far
Solving for state 96, time 9.400000e-09
  NL step  0, |residual|_2 = 3.663837e+00
  NL step  1, |residual|_2 = 4.776511e-08
  NL step  2, |residual|_2 = 2.493295e-15
Writing state 96 at time 9.500000e-09 to file bih.1.e; output a total of 96 states so far
Solving for state 97, time 9.500000e-09
  NL step  0, |residual|_2 = 3.662129e+00
  NL step  1, |residual|_2 = 4.764662e-08
  NL step  2, |residual|_2 = 2.383296e-15
Writing state 97 at time 9.600000e-09 to file bih.1.e; output a total of 97 states so far
Solving for state 98, time 9.600000e-09
  NL step  0, |residual|_2 = 3.660426e+00
  NL step  1, |residual|_2 = 4.752966e-08
  NL step  2, |residual|_2 = 2.499442e-15
Writing state 98 at time 9.700000e-09 to file bih.1.e; output a total of 98 states so far
Solving for state 99, time 9.700000e-09
  NL step  0, |residual|_2 = 3.658728e+00
  NL step  1, |residual|_2 = 4.741435e-08
  NL step  2, |residual|_2 = 2.524939e-15
Writing state 99 at time 9.800000e-09 to file bih.1.e; output a total of 99 states so far
Solving for state 100, time 9.800000e-09
  NL step  0, |residual|_2 = 3.657035e+00
  NL step  1, |residual|_2 = 4.730047e-08
  NL step  2, |residual|_2 = 2.548990e-15
Writing state 100 at time 9.900000e-09 to file bih.1.e; output a total of 100 states so far
Solving for state 101, time 9.900000e-09
  NL step  0, |residual|_2 = 3.655347e+00
  NL step  1, |residual|_2 = 4.718814e-08
  NL step  2, |residual|_2 = 2.498620e-15
Writing state 101 at time 1.000000e-08 to file bih.1.e; output a total of 101 states so far
Solving for state 102, time 1.000000e-08
  NL step  0, |residual|_2 = 3.653664e+00
  NL step  1, |residual|_2 = 4.707716e-08
  NL step  2, |residual|_2 = 2.539593e-15
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 102 states so far
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 103 states so far

-------------------------------------------------------------------------------------------------------------------
| Time:           Thu Aug 23 08:05:33 2012                                                                         |
| OS:             Darwin                                                                                           |
| HostName:       www.example.com                                                                                  |
| OS Release:     10.8.0                                                                                           |
| OS Version:     Darwin Kernel Version 10.8.0: Tue Jun  7 16:32:41 PDT 2011; root:xnu-1504.15.3~1/RELEASE_X86_64  |
| Machine:        x86_64                                                                                           |
| Username:       jwpeterson                                                                                       |
| Configuration:  ./configure run on Wed Aug 22 13:54:56 MDT 2012                                                  |
-------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.22066, Active time=3.04215                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0002      0.000170    0.0002      0.000170    0.01     0.01     |
|   build_constraint_matrix()        412672    0.1473      0.000000    0.1473      0.000000    4.84     4.84     |
|   build_sparsity()                 1         0.0020      0.001968    0.0026      0.002556    0.06     0.08     |
|   constrain_elem_matrix()          154624    0.0471      0.000000    0.0471      0.000000    1.55     1.55     |
|   constrain_elem_vector()          258048    0.0815      0.000000    0.0815      0.000000    2.68     2.68     |
|   create_dof_constraints()         1         0.0019      0.001926    0.0052      0.005171    0.06     0.17     |
|   distribute_dofs()                1         0.0003      0.000258    0.0013      0.001324    0.01     0.04     |
|   dof_indices()                    520196    0.2793      0.000001    0.2793      0.000001    9.18     9.18     |
|   enforce_constraints_exactly()    404       0.0156      0.000039    0.0156      0.000039    0.51     0.51     |
|   prepare_send_list()              1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                         1         0.0011      0.001065    0.0011      0.001065    0.04     0.04     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          103       0.2330      0.002262    0.2921      0.002836    7.66     9.60     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               103       0.0062      0.000060    0.0062      0.000060    0.20     0.20     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        407       0.0011      0.000003    0.0011      0.000003    0.04     0.04     |
|   init_shape_functions()           407       0.0046      0.000011    0.0046      0.000011    0.15     0.15     |
|   inverse_map()                    8         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             412676    0.2214      0.000001    0.2214      0.000001    7.28     7.28     |
|   compute_face_map()               2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_face_shape_functions()      2         0.0000      0.000013    0.0000      0.000013    0.00     0.00     |
|   init_reference_to_physical_map() 407       0.0011      0.000003    0.0011      0.000003    0.04     0.04     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0010      0.001026    0.0010      0.001026    0.03     0.03     |
|   renumber_nodes_and_elem()        2         0.0001      0.000044    0.0001      0.000044    0.00     0.00     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         103       0.0007      0.000007    0.2990      0.002903    0.02     9.83     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0005      0.000543    0.0005      0.000543    0.02     0.02     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   single_partition()               1         0.0001      0.000071    0.0001      0.000071    0.00     0.00     |
|                                                                                                                |
| PetscNonlinearSolver                                                                                           |
|   jacobian()                       151       0.7402      0.004902    1.0756      0.007123    24.33    35.36    |
|   residual()                       252       1.1155      0.004427    1.6836      0.006681    36.67    55.34    |
|   solve()                          101       0.1345      0.001331    2.8937      0.028650    4.42     95.12    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  1         0.0030      0.002970    0.0030      0.002970    0.10     0.10     |
|   operator()                       2         0.0000      0.000021    0.0001      0.000030    0.00     0.00     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 1         0.0022      0.002185    0.0028      0.002839    0.07     0.09     |
|   solve()                          101       0.0007      0.000007    2.8944      0.028657    0.02     95.14    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1760783   3.0422                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./miscellaneous_ex7-opt
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
