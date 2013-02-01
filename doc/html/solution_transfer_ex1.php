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
                             const Parameters& parameters,
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
        
          Mesh from_mesh;
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
                       <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
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
  
    Mesh from_mesh;
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
***************************************************************
* Running Example solution_transfer_ex1:
*  mpirun -np 2 example-dbg  
***************************************************************
 

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh4ElemE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh9DofObjectE reference count information:
|  Creations:    1
|  Destructions: 1
 ---------------------------------------------------------------------------- 

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:31:25 2013                    |
| OS:             Linux                                       |
| HostName:       lkirk-home                                  |
| OS Release:     3.2.0-35-generic                            |
| OS Version:     #55-Ubuntu SMP Wed Dec 5 17:42:16 UTC 2012  |
| Machine:        x86_64                                      |
| Username:       benkirk                                     |
| Configuration:  ./configure  '--prefix=/home/benkirk/codes/install'|
|  '--disable-glibcxx-debugging'                              |
|  '--enable-everything'                                      |
 -------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.05201, Active time=0.00031                                              |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Parallel                                                                                                  |
|   max(scalar)                 8         0.0001      0.000007    0.0001      0.000007    17.10    17.10    |
|   max(vector)                 2         0.0000      0.000015    0.0001      0.000035    9.68     22.58    |
|   min(bool)                   10        0.0001      0.000007    0.0001      0.000007    22.90    22.90    |
|   min(scalar)                 8         0.0001      0.000015    0.0001      0.000015    38.06    38.06    |
|   min(vector)                 2         0.0000      0.000019    0.0001      0.000042    12.26    27.10    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       30        0.0003                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example solution_transfer_ex1:
*  mpirun -np 2 example-dbg  
***************************************************************
***************************************************************
* Running Example solution_transfer_ex1:
*  mpirun -np 2 example-devel  
***************************************************************
 

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:31:26 2013                    |
| OS:             Linux                                       |
| HostName:       lkirk-home                                  |
| OS Release:     3.2.0-35-generic                            |
| OS Version:     #55-Ubuntu SMP Wed Dec 5 17:42:16 UTC 2012  |
| Machine:        x86_64                                      |
| Username:       benkirk                                     |
| Configuration:  ./configure  '--prefix=/home/benkirk/codes/install'|
|  '--disable-glibcxx-debugging'                              |
|  '--enable-everything'                                      |
 -------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.044522, Active time=0.000148                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
|                                                                                                           |
| Parallel                                                                                                  |
|   max(scalar)                 8         0.0000      0.000003    0.0000      0.000003    16.22    16.22    |
|   max(vector)                 2         0.0000      0.000005    0.0000      0.000014    6.76     18.92    |
|   min(bool)                   10        0.0000      0.000003    0.0000      0.000003    22.30    22.30    |
|   min(scalar)                 8         0.0001      0.000009    0.0001      0.000009    46.62    46.62    |
|   min(vector)                 2         0.0000      0.000006    0.0000      0.000015    8.11     20.27    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       30        0.0001                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example solution_transfer_ex1:
*  mpirun -np 2 example-devel  
***************************************************************
***************************************************************
* Running Example solution_transfer_ex1:
*  mpirun -np 2 example-opt  
***************************************************************
 
 
***************************************************************
* Done Running Example solution_transfer_ex1:
*  mpirun -np 2 example-opt  
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
