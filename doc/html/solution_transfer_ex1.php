<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("solution_transfer_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file main.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/explicit_system.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/libmesh_config.h"
        
        #ifdef LIBMESH_HAVE_DTK
        #include "libmesh/dtk_solution_transfer.h"
        #endif
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        Number initial_value(const Point& p,
                             const Parameters& /* parameters */,
                             const std::string&,
                             const std::string&)
        {
          return p(0)*p(0) + 1; // x^2 + 1
        }
        
        void initialize(EquationSystems& es,
                        const std::string& system_name)
        {
          ExplicitSystem & system = es.get_system&lt;ExplicitSystem&gt;(system_name);
          es.parameters.set&lt;Real&gt; ("time") = system.time = 0;
          system.project_solution(initial_value, NULL, es.parameters);
        }
        
        int main(int argc, char* argv[])
        {
          LibMeshInit init (argc, argv);
        
        #ifdef LIBMESH_HAVE_DTK
        
          Mesh from_mesh(init.comm());
          MeshTools::Generation::build_cube(from_mesh, 4, 4, 4, 0, 1, 0, 1, 0, 1, HEX8);
          from_mesh.print_info();
          EquationSystems from_es(from_mesh);
          System & from_sys = from_es.add_system&lt;ExplicitSystem&gt;("From");
          unsigned int from_var = from_sys.add_variable("from");
          from_sys.attach_init_function(initialize);
          from_es.init();
        
          ExodusII_IO(from_mesh).write_equation_systems("from.e", from_es);
        
          Mesh to_mesh;
          MeshTools::Generation::build_cube(to_mesh, 5, 5, 5, 0, 1, 0, 1, 0, 1, TET4);
          to_mesh.print_info();
          EquationSystems to_es(to_mesh);
          System & to_sys = to_es.add_system&lt;ExplicitSystem&gt;("To");
          unsigned int to_var = to_sys.add_variable("to");
          to_es.init();
        
          DTKSolutionTransfer dtk_transfer;
        
          dtk_transfer.transfer(from_sys.variable(from_var), to_sys.variable(to_var));
        
          to_es.update();
          ExodusII_IO(to_mesh).write_equation_systems("to.e", to_es);
        
        #endif
        
          return 0;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file main.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/explicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_config.h&quot;</FONT></B>
  
  #ifdef LIBMESH_HAVE_DTK
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dtk_solution_transfer.h&quot;</FONT></B>
  #endif
  
  using namespace libMesh;
  
  
  Number initial_value(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; <I><FONT COLOR="#B22222">/* parameters */</FONT></I>,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> p(0)*p(0) + 1; <I><FONT COLOR="#B22222">// x^2 + 1
</FONT></I>  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> initialize(EquationSystems&amp; es,
                  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    ExplicitSystem &amp; system = es.get_system&lt;ExplicitSystem&gt;(system_name);
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time = 0;
    system.project_solution(initial_value, NULL, es.parameters);
  }
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>* argv[])
  {
    LibMeshInit init (argc, argv);
  
  #ifdef LIBMESH_HAVE_DTK
  
    Mesh from_mesh(init.comm());
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube(from_mesh, 4, 4, 4, 0, 1, 0, 1, 0, 1, HEX8);
    from_mesh.print_info();
    EquationSystems from_es(from_mesh);
    System &amp; from_sys = from_es.add_system&lt;ExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;From&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> from_var = from_sys.add_variable(<B><FONT COLOR="#BC8F8F">&quot;from&quot;</FONT></B>);
    from_sys.attach_init_function(initialize);
    from_es.init();
  
    ExodusII_IO(from_mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;from.e&quot;</FONT></B>, from_es);
  
    Mesh to_mesh;
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube(to_mesh, 5, 5, 5, 0, 1, 0, 1, 0, 1, TET4);
    to_mesh.print_info();
    EquationSystems to_es(to_mesh);
    System &amp; to_sys = to_es.add_system&lt;ExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;To&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> to_var = to_sys.add_variable(<B><FONT COLOR="#BC8F8F">&quot;to&quot;</FONT></B>);
    to_es.init();
  
    DTKSolutionTransfer dtk_transfer;
  
    dtk_transfer.transfer(from_sys.variable(from_var), to_sys.variable(to_var));
  
    to_es.update();
    ExodusII_IO(to_mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;to.e&quot;</FONT></B>, to_es);
  
  #endif
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/solution_transfer/solution_transfer_ex1'
***************************************************************
* Running Example solution_transfer_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:50:29 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.053403, Active time=0.002006                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Parallel                                                                                                  |
|   max(scalar)                 8         0.0004      0.000049    0.0004      0.000049    19.54    19.54    |
|   max(vector)                 2         0.0001      0.000045    0.0004      0.000188    4.54     18.79    |
|   min(bool)                   10        0.0004      0.000044    0.0004      0.000044    21.88    21.88    |
|   min(scalar)                 8         0.0009      0.000112    0.0009      0.000112    44.77    44.77    |
|   min(vector)                 2         0.0002      0.000093    0.0005      0.000228    9.27     22.73    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       30        0.0020                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example solution_transfer_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/solution_transfer/solution_transfer_ex1'
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
