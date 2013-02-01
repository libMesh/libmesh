<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex8",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file meshless_interpolation_function.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
        #define LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
        
</pre>
</div>
<div class = "comment">
Local Includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/function_base.h"
        #include "libmesh/meshfree_interpolation.h"
        #include "libmesh/threads.h"
        
        
</pre>
</div>
<div class = "comment">
C++ includes
</div>

<div class ="fragment">
<pre>
        #include &lt;cstddef&gt;
        
        namespace libMesh
        {
        
        
        
</pre>
</div>
<div class = "comment">
Forward Declarations
</div>

<div class ="fragment">
<pre>
        template &lt;typename T&gt;
        class DenseVector;
        
        
</pre>
</div>
<div class = "comment">
------------------------------------------------------------
MeshlessInterpolationFunction class definition
</div>

<div class ="fragment">
<pre>
        class MeshlessInterpolationFunction : public FunctionBase&lt;Number&gt;
        {
        private:
          const MeshfreeInterpolation &_mfi;
          mutable std::vector&lt;Point&gt; _pts;
          mutable std::vector&lt;Number&gt; _vals;
          Threads::spin_mutex &_mutex;
          
        public:
        
          /**
           * Constructor.  Requires a \p \pMeshlessInterpolation object.
           */
          MeshlessInterpolationFunction (const MeshfreeInterpolation &mfi,
        				 Threads::spin_mutex &mutex) :
          _mfi (mfi),
          _mutex(mutex)
          {}
          
        
          /**
           * The actual initialization process.
           */
          void init ();
        
          /**
           * Clears the function.
           */
          void clear ();
        
          /**
           * Returns a new deep copy of the function.
           */
          virtual AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone () const;
        
          /**
           * @returns the value at point \p p and time
           * \p time, which defaults to zero.
           */
          Number operator() (const Point& p,
        		     const Real time=0.);
        
          /**
           * Like before, but returns the values in a
           * writable reference.
           */
          void operator() (const Point& p,
        		   const Real time,
        		   DenseVector&lt;Number&gt;& output);
        
        };
        
        
        
</pre>
</div>
<div class = "comment">
------------------------------------------------------------
MeshlessInterpolationFunction inline methods
</div>

<div class ="fragment">
<pre>
        inline
        Number MeshlessInterpolationFunction::operator() (const Point& p,
        						  const Real /* time */)
        {
          _pts.clear();
          _pts.push_back(p);
          _vals.resize(1);
        
          Threads::spin_mutex::scoped_lock lock(_mutex);
        
          _mfi.interpolate_field_data (_mfi.field_variables(),
        			       _pts, _vals);
        
          return _vals.front();
        }
        
        
        
        inline
        void MeshlessInterpolationFunction::operator() (const Point& p,
        						const Real time,
        						DenseVector&lt;Number&gt;& output)
        {
          output.resize(1);
          output(0) = (*this)(p,time);
          return;
        }
        
        
        
        inline
        void MeshlessInterpolationFunction::init ()
        {
        }
        
        
        
        inline
        void MeshlessInterpolationFunction::clear ()
        {
        }
        
        
        
        inline
        AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
        MeshlessInterpolationFunction::clone () const
        { 
          return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (new MeshlessInterpolationFunction (_mfi, _mutex) );
        }
        
        
        } // namespace libMesh
        
        
        #endif // LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
        
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex8.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/meshfree_interpolation.h"
        #include "libmesh/mesh.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/tecplot_io.h"
        #include "libmesh/threads.h"
        #include "meshless_interpolation_function.h"
        
</pre>
</div>
<div class = "comment">
C++ includes
</div>

<div class ="fragment">
<pre>
        #include &lt;cstdlib&gt;
        
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        void create_random_point_cloud (const unsigned int Npts,
        				std::vector&lt;Point&gt; &pts,
        				const Real max_range = 10)
        {
          std::cout &lt;&lt; "Generating "&lt;&lt; Npts &lt;&lt; " point cloud...";
          pts.resize(Npts);
        
          for (size_t i=0;i&lt;Npts;i++)
            {
              pts[i](0) = max_range * (std::rand() % 1000) / Real(1000);
              pts[i](1) = max_range * (std::rand() % 1000) / Real(1000);
              pts[i](2) = max_range * (std::rand() % 1000) / Real(1000);
            }
          std::cout &lt;&lt; "done\n";
        }
        
        
        
        Real exact_solution_u (const Point &p)
        {
          const Real
            x = p(0),
            y = p(1),
            z = p(2);
          
          return (x*x*x   +
        	  y*y*y*y +
        	  z*z*z*z*z);
        }
        
        
        
        Real exact_solution_v (const Point &p)
        {
          const Real
            x = p(0),
            y = p(1),
            z = p(2);
          
          return (x*x   +
        	  y*y +
        	  z*z*z);
        }
        
        Number exact_value (const Point& p,
                            const Parameters&,
                            const std::string&,
                            const std::string&)
        {
          return exact_solution_v(p);
        }
        
</pre>
</div>
<div class = "comment">
We now define the function which provides the
initialization routines for the "Convection-Diffusion"
system.  This handles things like setting initial
conditions and boundary conditions.
</div>

<div class ="fragment">
<pre>
        void init_sys(EquationSystems& es,
                      const std::string& system_name)
        {
</pre>
</div>
<div class = "comment">
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          System & system =
            es.get_system&lt;System&gt;(system_name);
          
          system.project_solution(exact_value, NULL, es.parameters);
        }
        
        
        
        
        int main(int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Skip this example if we do not meet certain requirements
</div>

<div class ="fragment">
<pre>
        #if LIBMESH_DIM != 3
          std::cout &lt;&lt; "This example requires 3D support - skipping.\n";
          return 77;
        #endif
        #ifndef LIBMESH_HAVE_ZLIB_H
          std::cout &lt;&lt; "This example requires zlib support - skipping.\n";
          return 77;
        #endif
        
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
          {
</pre>
</div>
<div class = "comment">
Demonstration case 1
</div>

<div class ="fragment">
<pre>
            {
              std::vector&lt;Point&gt;       tgt_pts;
              std::vector&lt;Number&gt;      tgt_data;
              std::vector&lt;std::string&gt; field_vars;
              
              field_vars.push_back("u");
              field_vars.push_back("v");
              
              InverseDistanceInterpolation&lt;3&gt; idi (/* n_interp_pts = */ 8,
        					   /* power =        */ 2);
              
              idi.set_field_variables (field_vars);
              
              create_random_point_cloud (1e5,
        				 idi.get_source_points());
              
</pre>
</div>
<div class = "comment">
Explicitly set the data values we will interpolate from
</div>

<div class ="fragment">
<pre>
              {
        	const std::vector&lt;Point&gt; &src_pts  (idi.get_source_points());
        	std::vector&lt;Number&gt;      &src_vals (idi.get_source_vals());
        	
        	src_vals.clear(); src_vals.reserve(src_pts.size());
        	
        	for (std::vector&lt;Point&gt;::const_iterator pt_it=src_pts.begin();
        	     pt_it != src_pts.end(); ++pt_it)
        	  {
        	    src_vals.push_back (exact_solution_u (*pt_it));
        	    src_vals.push_back (exact_solution_v (*pt_it));
        	  }	  
              }
        
              std::cout &lt;&lt; idi;
              
</pre>
</div>
<div class = "comment">
Interpolate to some other random points, and evaluate the result
</div>

<div class ="fragment">
<pre>
              {
        	create_random_point_cloud (10,
        				   tgt_pts);
        	
        	idi.interpolate_field_data (field_vars,
        				    tgt_pts,
        				    tgt_data);
        	
              
        	std::vector&lt;Number&gt;::const_iterator v_it=tgt_data.begin();
        	
        	for (std::vector&lt;Point&gt;::const_iterator  p_it=tgt_pts.begin();
        	     p_it!=tgt_pts.end(); ++p_it)
        	  {
        	    std::cout &lt;&lt; "\nAt target point " &lt;&lt; *p_it
        		      &lt;&lt; "\n u_interp=" &lt;&lt; *v_it
        		      &lt;&lt; ", u_exact="  &lt;&lt; exact_solution_u(*p_it);
        	    ++v_it;
        	    std::cout &lt;&lt; "\n v_interp=" &lt;&lt; *v_it
        		      &lt;&lt; ", v_exact="  &lt;&lt; exact_solution_v(*p_it)
        		      &lt;&lt; std::endl;
        	    ++v_it;
        	  }
              }
            }
        
        
</pre>
</div>
<div class = "comment">
Demonstration case 2
</div>

<div class ="fragment">
<pre>
            {
              Mesh mesh_a, mesh_b;
        
              mesh_a.read("struct.ucd.gz"); mesh_b.read("unstruct.ucd.gz");
        
</pre>
</div>
<div class = "comment">
Create equation systems objects.
</div>

<div class ="fragment">
<pre>
              EquationSystems
        	es_a(mesh_a), es_b(mesh_b);
        
              System
        	&sys_a = es_a.add_system&lt;System&gt;("src_system"),
        	&sys_b = es_b.add_system&lt;System&gt;("dest_system");
        
              sys_a.add_variable ("Cp", FIRST);
              sys_b.add_variable ("Cp", FIRST);
        
              sys_a.attach_init_function (init_sys);
              es_a.init();
              
</pre>
</div>
<div class = "comment">
Write out the initial conditions.
</div>

<div class ="fragment">
<pre>
              TecplotIO(mesh_a).write_equation_systems ("src.dat",
        						es_a);
        
              InverseDistanceInterpolation&lt;3&gt; idi (/* n_interp_pts = */ 4,
        					   /* power =        */ 2);
        
              std::vector&lt;Point&gt;  &src_pts  (idi.get_source_points());
              std::vector&lt;Number&gt; &src_vals (idi.get_source_vals());
              std::vector&lt;std::string&gt; field_vars;      
              field_vars.push_back("Cp");
              idi.set_field_variables(field_vars);
        
</pre>
</div>
<div class = "comment">
We now will loop over every node in the source mesh
and add it to a source point list, along with the solution
</div>

<div class ="fragment">
<pre>
              {
        	MeshBase::const_node_iterator nd  = mesh_a.local_nodes_begin();
        	MeshBase::const_node_iterator end = mesh_a.local_nodes_end();
        
        	for (; nd!=end; ++nd)
        	  {
        	    const Node *node(*nd);
        	    src_pts.push_back(*node);
        	    src_vals.push_back(sys_a.current_solution(node-&gt;dof_number(0,0,0)));
        	  }	  	
              }
        
</pre>
</div>
<div class = "comment">
We have only set local values - prepare for use by gathering remote gata
</div>

<div class ="fragment">
<pre>
              idi.prepare_for_use();
        
</pre>
</div>
<div class = "comment">
Create a MeshlessInterpolationFunction that uses our InverseDistanceInterpolation
object.  Since each MeshlessInterpolationFunction shares the same InverseDistanceInterpolation
object in a threaded environment we must also provide a locking mechanism.
</div>

<div class ="fragment">
<pre>
              Threads::spin_mutex mutex;
              MeshlessInterpolationFunction mif(idi, mutex);
        
</pre>
</div>
<div class = "comment">
project the solution onto system b
</div>

<div class ="fragment">
<pre>
              es_b.init();
              sys_b.project_solution (&mif);
              
</pre>
</div>
<div class = "comment">
Write the result
</div>

<div class ="fragment">
<pre>
              TecplotIO(mesh_b).write_equation_systems ("dest.dat",
        						es_b);
            }
        
        
            
          }
          return 0;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file meshless_interpolation_function.h without comments: </h1> 
<pre> 
  #ifndef LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
  #define LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/function_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/meshfree_interpolation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/threads.h&quot;</FONT></B>
  
  
  #include &lt;cstddef&gt;
  
  namespace libMesh
  {
  
  
  
  <B><FONT COLOR="#228B22">template</FONT></B> &lt;typename T&gt;
  <B><FONT COLOR="#228B22">class</FONT></B> DenseVector;
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> MeshlessInterpolationFunction : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">private</FONT></B>:
    <B><FONT COLOR="#228B22">const</FONT></B> MeshfreeInterpolation &amp;_mfi;
    mutable std::vector&lt;Point&gt; _pts;
    mutable std::vector&lt;Number&gt; _vals;
    <B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex &amp;_mutex;
    
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    <I><FONT COLOR="#B22222">/**
     * Constructor.  Requires a \p \pMeshlessInterpolation object.
     */</FONT></I>
    MeshlessInterpolationFunction (<B><FONT COLOR="#228B22">const</FONT></B> MeshfreeInterpolation &amp;mfi,
  				 <B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex &amp;mutex) :
    _mfi (mfi),
    _mutex(mutex)
    {}
    
  
    <I><FONT COLOR="#B22222">/**
     * The actual initialization process.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> init ();
  
    <I><FONT COLOR="#B22222">/**
     * Clears the function.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> clear ();
  
    <I><FONT COLOR="#B22222">/**
     * Returns a new deep copy of the function.
     */</FONT></I>
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone () <B><FONT COLOR="#228B22">const</FONT></B>;
  
    <I><FONT COLOR="#B22222">/**
     * @returns the value at point \p p and time
     * \p time, which defaults to zero.
     */</FONT></I>
    Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real time=0.);
  
    <I><FONT COLOR="#B22222">/**
     * Like before, but returns the values in a
     * writable reference.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  		   <B><FONT COLOR="#228B22">const</FONT></B> Real time,
  		   DenseVector&lt;Number&gt;&amp; output);
  
  };
  
  
  
  <B><FONT COLOR="#228B22">inline</FONT></B>
  Number MeshlessInterpolationFunction::<B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  						  <B><FONT COLOR="#228B22">const</FONT></B> Real <I><FONT COLOR="#B22222">/* time */</FONT></I>)
  {
    _pts.clear();
    _pts.push_back(p);
    _vals.resize(1);
  
    <B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex::scoped_lock lock(_mutex);
  
    _mfi.interpolate_field_data (_mfi.field_variables(),
  			       _pts, _vals);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> _vals.front();
  }
  
  
  
  <B><FONT COLOR="#228B22">inline</FONT></B>
  <B><FONT COLOR="#228B22">void</FONT></B> MeshlessInterpolationFunction::<B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  						<B><FONT COLOR="#228B22">const</FONT></B> Real time,
  						DenseVector&lt;Number&gt;&amp; output)
  {
    output.resize(1);
    output(0) = (*<B><FONT COLOR="#A020F0">this</FONT></B>)(p,time);
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  
  
  <B><FONT COLOR="#228B22">inline</FONT></B>
  <B><FONT COLOR="#228B22">void</FONT></B> MeshlessInterpolationFunction::init ()
  {
  }
  
  
  
  <B><FONT COLOR="#228B22">inline</FONT></B>
  <B><FONT COLOR="#228B22">void</FONT></B> MeshlessInterpolationFunction::clear ()
  {
  }
  
  
  
  <B><FONT COLOR="#228B22">inline</FONT></B>
  AutoPtr&lt;FunctionBase&lt;Number&gt; &gt;
  <B><FONT COLOR="#5F9EA0">MeshlessInterpolationFunction</FONT></B>::clone () <B><FONT COLOR="#228B22">const</FONT></B>
  { 
    <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (<B><FONT COLOR="#A020F0">new</FONT></B> MeshlessInterpolationFunction (_mfi, _mutex) );
  }
  
  
  } <I><FONT COLOR="#B22222">// namespace libMesh
</FONT></I>  
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_MESHLESS_INTERPOLATION_FUNCTION_H
</FONT></I>  
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex8.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/meshfree_interpolation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/threads.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;meshless_interpolation_function.h&quot;</FONT></B>
  
  #include &lt;cstdlib&gt;
  
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> create_random_point_cloud (<B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> Npts,
  				<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt; &amp;pts,
  				<B><FONT COLOR="#228B22">const</FONT></B> Real max_range = 10)
  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Generating &quot;</FONT></B>&lt;&lt; Npts &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; point cloud...&quot;</FONT></B>;
    pts.resize(Npts);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (size_t i=0;i&lt;Npts;i++)
      {
        pts[i](0) = max_range * (std::rand() % 1000) / Real(1000);
        pts[i](1) = max_range * (std::rand() % 1000) / Real(1000);
        pts[i](2) = max_range * (std::rand() % 1000) / Real(1000);
      }
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;done\n&quot;</FONT></B>;
  }
  
  
  
  Real exact_solution_u (<B><FONT COLOR="#228B22">const</FONT></B> Point &amp;p)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real
      x = p(0),
      y = p(1),
      z = p(2);
    
    <B><FONT COLOR="#A020F0">return</FONT></B> (x*x*x   +
  	  y*y*y*y +
  	  z*z*z*z*z);
  }
  
  
  
  Real exact_solution_v (<B><FONT COLOR="#228B22">const</FONT></B> Point &amp;p)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real
      x = p(0),
      y = p(1),
      z = p(2);
    
    <B><FONT COLOR="#A020F0">return</FONT></B> (x*x   +
  	  y*y +
  	  z*z*z);
  }
  
  Number exact_value (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> exact_solution_v(p);
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_sys(EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    System &amp; system =
      es.get_system&lt;System&gt;(system_name);
    
    system.project_solution(exact_value, NULL, es.parameters);
  }
  
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM != 3
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;This example requires 3D support - skipping.\n&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">return</FONT></B> 77;
  #endif
  #ifndef LIBMESH_HAVE_ZLIB_H
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;This example requires zlib support - skipping.\n&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">return</FONT></B> 77;
  #endif
  
    LibMeshInit init (argc, argv);
    {
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt;       tgt_pts;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;      tgt_data;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; field_vars;
        
        field_vars.push_back(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
        field_vars.push_back(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
        
        InverseDistanceInterpolation&lt;3&gt; idi (<I><FONT COLOR="#B22222">/* n_interp_pts = */</FONT></I> 8,
  					   <I><FONT COLOR="#B22222">/* power =        */</FONT></I> 2);
        
        idi.set_field_variables (field_vars);
        
        create_random_point_cloud (1e5,
  				 idi.get_source_points());
        
        {
  	<B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;src_pts  (idi.get_source_points());
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;      &amp;src_vals (idi.get_source_vals());
  	
  	src_vals.clear(); src_vals.reserve(src_pts.size());
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;Point&gt;::const_iterator pt_it=src_pts.begin();
  	     pt_it != src_pts.end(); ++pt_it)
  	  {
  	    src_vals.push_back (exact_solution_u (*pt_it));
  	    src_vals.push_back (exact_solution_v (*pt_it));
  	  }	  
        }
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; idi;
        
        {
  	create_random_point_cloud (10,
  				   tgt_pts);
  	
  	idi.interpolate_field_data (field_vars,
  				    tgt_pts,
  				    tgt_data);
  	
        
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;::const_iterator v_it=tgt_data.begin();
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;Point&gt;::const_iterator  p_it=tgt_pts.begin();
  	     p_it!=tgt_pts.end(); ++p_it)
  	  {
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\nAt target point &quot;</FONT></B> &lt;&lt; *p_it
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n u_interp=&quot;</FONT></B> &lt;&lt; *v_it
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, u_exact=&quot;</FONT></B>  &lt;&lt; exact_solution_u(*p_it);
  	    ++v_it;
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n v_interp=&quot;</FONT></B> &lt;&lt; *v_it
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, v_exact=&quot;</FONT></B>  &lt;&lt; exact_solution_v(*p_it)
  		      &lt;&lt; std::endl;
  	    ++v_it;
  	  }
        }
      }
  
  
      {
        Mesh mesh_a, mesh_b;
  
        mesh_a.read(<B><FONT COLOR="#BC8F8F">&quot;struct.ucd.gz&quot;</FONT></B>); mesh_b.read(<B><FONT COLOR="#BC8F8F">&quot;unstruct.ucd.gz&quot;</FONT></B>);
  
        EquationSystems
  	es_a(mesh_a), es_b(mesh_b);
  
        System
  	&amp;sys_a = es_a.add_system&lt;System&gt;(<B><FONT COLOR="#BC8F8F">&quot;src_system&quot;</FONT></B>),
  	&amp;sys_b = es_b.add_system&lt;System&gt;(<B><FONT COLOR="#BC8F8F">&quot;dest_system&quot;</FONT></B>);
  
        sys_a.add_variable (<B><FONT COLOR="#BC8F8F">&quot;Cp&quot;</FONT></B>, FIRST);
        sys_b.add_variable (<B><FONT COLOR="#BC8F8F">&quot;Cp&quot;</FONT></B>, FIRST);
  
        sys_a.attach_init_function (init_sys);
        es_a.init();
        
        TecplotIO(mesh_a).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;src.dat&quot;</FONT></B>,
  						es_a);
  
        InverseDistanceInterpolation&lt;3&gt; idi (<I><FONT COLOR="#B22222">/* n_interp_pts = */</FONT></I> 4,
  					   <I><FONT COLOR="#B22222">/* power =        */</FONT></I> 2);
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt;  &amp;src_pts  (idi.get_source_points());
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; &amp;src_vals (idi.get_source_vals());
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; field_vars;      
        field_vars.push_back(<B><FONT COLOR="#BC8F8F">&quot;Cp&quot;</FONT></B>);
        idi.set_field_variables(field_vars);
  
        {
  	<B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator nd  = mesh_a.local_nodes_begin();
  	<B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator end = mesh_a.local_nodes_end();
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (; nd!=end; ++nd)
  	  {
  	    <B><FONT COLOR="#228B22">const</FONT></B> Node *node(*nd);
  	    src_pts.push_back(*node);
  	    src_vals.push_back(sys_a.current_solution(node-&gt;dof_number(0,0,0)));
  	  }	  	
        }
  
        idi.prepare_for_use();
  
        <B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex mutex;
        MeshlessInterpolationFunction mif(idi, mutex);
  
        es_b.init();
        sys_b.project_solution (&amp;mif);
        
        TecplotIO(mesh_b).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;dest.dat&quot;</FONT></B>,
  						es_b);
      }
  
  
      
    }
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example miscellaneous_ex8:
*  mpirun -np 2 example-devel  
***************************************************************
 
Generating 100000 point cloud...done
MeshfreeInterpolation
 n_source_points()=100000
 n_field_variables()=2
  variables = u v 
Generating 10 point cloud...done
*** Warning, This code is untested, experimental, or likely to see future API changes: src/solution_transfer/meshfree_interpolation.C, line 283, compiled Feb  1 2013 at 09:03:51 ***

At target point (x,y,z)=(    9.61,     7.45,     5.24)
 u_interp=8241.22, u_exact=7918.57
 v_interp=298.061, v_exact=291.732

At target point (x,y,z)=(     3.5,     9.71,     8.11)
 u_interp=44209.4, u_exact=44016
 v_interp=641.42, v_exact=639.946

At target point (x,y,z)=(     5.3,      8.4,     4.36)
 u_interp=6706.14, u_exact=6703.14
 v_interp=183.748, v_exact=181.532

At target point (x,y,z)=(    3.57,     5.66,     5.31)
 u_interp=5324.74, u_exact=5293.34
 v_interp=194.632, v_exact=194.502

At target point (x,y,z)=(    7.66,     3.88,     8.31)
 u_interp=39033.6, u_exact=40304.4
 v_interp=635.944, v_exact=647.586

At target point (x,y,z)=(    9.54,     7.83,     1.32)
 u_interp=4706.06, u_exact=4631.04
 v_interp=155.44, v_exact=154.62

At target point (x,y,z)=(    8.32,     1.21,     6.09)
 u_interp=8806.23, u_exact=8955.03
 v_interp=293.686, v_exact=296.553

At target point (x,y,z)=(    0.32,     2.44,     1.59)
 u_interp=45.7344, u_exact=45.6403
 v_interp=10.1975, v_exact=10.0757

At target point (x,y,z)=(    2.14,     4.04,     8.72)
 u_interp=51856, u_exact=50693.8
 v_interp=692.691, v_exact=683.956

At target point (x,y,z)=(     9.4,     0.51,     2.22)
 u_interp=883.149, u_exact=884.574
 v_interp=99.3335, v_exact=99.5611

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:35:07 2013                    |
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
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.660645, Active time=0.567406                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 2         0.0097      0.004854    0.0106      0.005279    1.71     1.86     |
|   create_dof_constraints()     2         0.0031      0.001526    0.0031      0.001526    0.54     0.54     |
|   distribute_dofs()            2         0.0137      0.006849    0.0364      0.018178    2.41     6.41     |
|   dof_indices()                9157      0.0250      0.000003    0.0250      0.000003    4.40     4.40     |
|   prepare_send_list()          2         0.0000      0.000017    0.0000      0.000017    0.01     0.01     |
|   reinit()                     2         0.0220      0.011000    0.0220      0.011000    3.88     3.88     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      2         0.0071      0.003529    0.0192      0.009587    1.24     3.38     |
|                                                                                                            |
| InverseDistanceInterpolation<>                                                                             |
|   construct_kd_tree()          2         0.1806      0.090307    0.1806      0.090307    31.83    31.83    |
|   interpolate_field_data()     8437      0.0198      0.000002    0.0198      0.000002    3.49     3.49     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0200      0.009977    0.0206      0.010304    3.52     3.63     |
|   read()                       2         0.0431      0.021550    0.0431      0.021550    7.60     7.60     |
|   renumber_nodes_and_elem()    4         0.0021      0.000535    0.0021      0.000535    0.38     0.38     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   broadcast()                  2         0.0102      0.005122    0.0180      0.008993    1.81     3.17     |
|   compute_hilbert_indices()    4         0.0736      0.018388    0.0736      0.018388    12.96    12.96    |
|   find_global_indices()        4         0.0093      0.002322    0.0867      0.021683    1.64     15.29    |
|   parallel_sort()              4         0.0017      0.000435    0.0033      0.000821    0.31     0.58     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     2         0.0001      0.000046    0.0648      0.032406    0.02     11.42    |
|                                                                                                            |
| MeshfreeInterpolation                                                                                      |
|   gather_remote_data()         1         0.0001      0.000132    0.0003      0.000274    0.02     0.05     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0416      0.020810    0.0831      0.041538    7.34     14.64    |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  18        0.0002      0.000010    0.0002      0.000012    0.03     0.04     |
|   broadcast()                  12        0.0020      0.000163    0.0020      0.000163    0.35     0.35     |
|   max(scalar)                  204       0.0006      0.000003    0.0006      0.000003    0.10     0.10     |
|   max(vector)                  47        0.0002      0.000005    0.0006      0.000013    0.04     0.11     |
|   min(bool)                    241       0.0006      0.000003    0.0006      0.000003    0.11     0.11     |
|   min(scalar)                  198       0.0075      0.000038    0.0075      0.000038    1.32     1.32     |
|   min(vector)                  47        0.0003      0.000007    0.0007      0.000016    0.06     0.13     |
|   probe()                      21        0.0002      0.000009    0.0002      0.000009    0.03     0.03     |
|   receive()                    21        0.0006      0.000028    0.0008      0.000037    0.10     0.14     |
|   send()                       21        0.0001      0.000005    0.0001      0.000005    0.02     0.02     |
|   send_receive()               28        0.0001      0.000005    0.0009      0.000034    0.02     0.17     |
|   sum()                        20        0.0019      0.000093    0.0020      0.000099    0.33     0.35     |
|                                                                                                            |
| Parallel::Request                                                                                          |
|   wait()                       21        0.0001      0.000002    0.0001      0.000002    0.01     0.01     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0052      0.002580    0.0055      0.002729    0.91     0.96     |
|   set_parent_processor_ids()   2         0.0020      0.001011    0.0020      0.001011    0.36     0.36     |
|                                                                                                            |
| System                                                                                                     |
|   project_vector()             2         0.0176      0.008803    0.0555      0.027763    3.10     9.79     |
|                                                                                                            |
| TecplotIO                                                                                                  |
|   write_nodal_data()           2         0.0454      0.022722    0.0454      0.022722    8.01     8.01     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        18542     0.5674                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex8:
*  mpirun -np 2 example-devel  
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
