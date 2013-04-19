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
        #include "libmesh/radial_basis_interpolation.h"
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
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "3D support");
        #ifndef LIBMESH_HAVE_EIGEN
          libmesh_example_assert(false, "--enable-eigen");
        #endif
        #ifndef LIBMESH_HAVE_ZLIB_H
          libmesh_example_assert(false, "--enable-zlib");
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
              std::vector&lt;Number&gt;      tgt_data_idi, tgt_data_rbi;
              std::vector&lt;std::string&gt; field_vars;
        
              field_vars.push_back("u");
              field_vars.push_back("v");
        
              InverseDistanceInterpolation&lt;3&gt; idi (init.comm(),
        					   /* n_interp_pts = */ 8,
        					   /* power =        */ 2);
        
              RadialBasisInterpolation&lt;3&gt; rbi (init.comm());
        
              idi.set_field_variables (field_vars);
              rbi.set_field_variables (field_vars);
        
              create_random_point_cloud (100,
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
        
        	src_vals.clear(); src_vals.reserve(2*src_pts.size());
        
        	for (std::vector&lt;Point&gt;::const_iterator pt_it=src_pts.begin();
        	     pt_it != src_pts.end(); ++pt_it)
        	  {
        	    src_vals.push_back (exact_solution_u (*pt_it));
        	    src_vals.push_back (exact_solution_v (*pt_it));
        	  }
              }
        
</pre>
</div>
<div class = "comment">
give rbi the same info as idi
</div>

<div class ="fragment">
<pre>
              rbi.get_source_points() = idi.get_source_points();
              rbi.get_source_vals()   = idi.get_source_vals();
        
              idi.prepare_for_use();
              rbi.prepare_for_use();
        
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
        
</pre>
</div>
<div class = "comment">
tgt_pts = rbi.get_source_points();


<br><br></div>

<div class ="fragment">
<pre>
                idi.interpolate_field_data (field_vars,
        				    tgt_pts,
        				    tgt_data_idi);
        
        	rbi.interpolate_field_data (field_vars,
        				    tgt_pts,
        				    tgt_data_rbi);
        
              	std::vector&lt;Number&gt;::const_iterator
        	  v_idi=tgt_data_idi.begin(),
        	  v_rbi=tgt_data_rbi.begin();
        
        	for (std::vector&lt;Point&gt;::const_iterator  p_it=tgt_pts.begin();
        	     p_it!=tgt_pts.end(); ++p_it)
        	  {
        	    std::cout &lt;&lt; "\nAt target point " &lt;&lt; *p_it
        		      &lt;&lt; "\n u_interp_idi="   &lt;&lt; *v_idi
        		      &lt;&lt; ", u_interp_rbi="    &lt;&lt; *v_rbi
        		      &lt;&lt; ", u_exact="         &lt;&lt; exact_solution_u(*p_it);
        	    ++v_idi;
        	    ++v_rbi;
        	    std::cout &lt;&lt; "\n v_interp_idi=" &lt;&lt; *v_idi
        		      &lt;&lt; ", v_interp_rbi="  &lt;&lt; *v_rbi
        		      &lt;&lt; ", v_exact="       &lt;&lt; exact_solution_v(*p_it)
        		      &lt;&lt; std::endl;
        	    ++v_idi;
        	    ++v_rbi;
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
              Mesh mesh_a(init.comm()), mesh_b(init.comm());
        
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
        
              InverseDistanceInterpolation&lt;3&gt; idi (init.comm(),
        					   /* n_interp_pts = */ 4,
        					   /* power =        */ 2);
              RadialBasisInterpolation&lt;3&gt; rbi (init.comm());
        
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
        
        	rbi.set_field_variables(field_vars);
        	rbi.get_source_points() = idi.get_source_points();
        	rbi.get_source_vals()   = idi.get_source_vals();
              }
        
</pre>
</div>
<div class = "comment">
We have only set local values - prepare for use by gathering remote gata
</div>

<div class ="fragment">
<pre>
              idi.prepare_for_use();
              rbi.prepare_for_use();
        
</pre>
</div>
<div class = "comment">
Create a MeshlessInterpolationFunction that uses our InverseDistanceInterpolation
object.  Since each MeshlessInterpolationFunction shares the same InverseDistanceInterpolation
object in a threaded environment we must also provide a locking mechanism.
</div>

<div class ="fragment">
<pre>
              {
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
                TecplotIO(mesh_b).write_equation_systems ("dest_idi.dat",
        						  es_b);
              }
        
</pre>
</div>
<div class = "comment">
Create a MeshlessInterpolationFunction that uses our RadialBasisInterpolation
object.  Since each MeshlessInterpolationFunction shares the same RadialBasisInterpolation
object in a threaded environment we must also provide a locking mechanism.
</div>

<div class ="fragment">
<pre>
              {
        	Threads::spin_mutex mutex;
        	MeshlessInterpolationFunction mif(rbi, mutex);
        
</pre>
</div>
<div class = "comment">
project the solution onto system b
</div>

<div class ="fragment">
<pre>
                sys_b.project_solution (&mif);
        
</pre>
</div>
<div class = "comment">
Write the result
</div>

<div class ="fragment">
<pre>
                TecplotIO(mesh_b).write_equation_systems ("dest_rbi.dat",
        						  es_b);
              }
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
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex8.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/meshfree_interpolation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/radial_basis_interpolation.h&quot;</FONT></B>
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
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  #ifndef LIBMESH_HAVE_EIGEN
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-eigen&quot;</FONT></B>);
  #endif
  #ifndef LIBMESH_HAVE_ZLIB_H
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-zlib&quot;</FONT></B>);
  #endif
  
    LibMeshInit init (argc, argv);
    {
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt;       tgt_pts;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;      tgt_data_idi, tgt_data_rbi;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;std::string&gt; field_vars;
  
        field_vars.push_back(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
        field_vars.push_back(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
  
        InverseDistanceInterpolation&lt;3&gt; idi (init.comm(),
  					   <I><FONT COLOR="#B22222">/* n_interp_pts = */</FONT></I> 8,
  					   <I><FONT COLOR="#B22222">/* power =        */</FONT></I> 2);
  
        RadialBasisInterpolation&lt;3&gt; rbi (init.comm());
  
        idi.set_field_variables (field_vars);
        rbi.set_field_variables (field_vars);
  
        create_random_point_cloud (100,
  				 idi.get_source_points());
  
  
        {
  	<B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt; &amp;src_pts  (idi.get_source_points());
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;      &amp;src_vals (idi.get_source_vals());
  
  	src_vals.clear(); src_vals.reserve(2*src_pts.size());
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;Point&gt;::const_iterator pt_it=src_pts.begin();
  	     pt_it != src_pts.end(); ++pt_it)
  	  {
  	    src_vals.push_back (exact_solution_u (*pt_it));
  	    src_vals.push_back (exact_solution_v (*pt_it));
  	  }
        }
  
        rbi.get_source_points() = idi.get_source_points();
        rbi.get_source_vals()   = idi.get_source_vals();
  
        idi.prepare_for_use();
        rbi.prepare_for_use();
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; idi;
  
        {
  	create_random_point_cloud (10,
  				   tgt_pts);
  
  
  	idi.interpolate_field_data (field_vars,
  				    tgt_pts,
  				    tgt_data_idi);
  
  	rbi.interpolate_field_data (field_vars,
  				    tgt_pts,
  				    tgt_data_rbi);
  
        	<B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt;::const_iterator
  	  v_idi=tgt_data_idi.begin(),
  	  v_rbi=tgt_data_rbi.begin();
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (std::vector&lt;Point&gt;::const_iterator  p_it=tgt_pts.begin();
  	     p_it!=tgt_pts.end(); ++p_it)
  	  {
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\nAt target point &quot;</FONT></B> &lt;&lt; *p_it
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n u_interp_idi=&quot;</FONT></B>   &lt;&lt; *v_idi
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, u_interp_rbi=&quot;</FONT></B>    &lt;&lt; *v_rbi
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, u_exact=&quot;</FONT></B>         &lt;&lt; exact_solution_u(*p_it);
  	    ++v_idi;
  	    ++v_rbi;
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n v_interp_idi=&quot;</FONT></B> &lt;&lt; *v_idi
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, v_interp_rbi=&quot;</FONT></B>  &lt;&lt; *v_rbi
  		      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, v_exact=&quot;</FONT></B>       &lt;&lt; exact_solution_v(*p_it)
  		      &lt;&lt; std::endl;
  	    ++v_idi;
  	    ++v_rbi;
  	  }
        }
      }
  
  
      {
        Mesh mesh_a(init.comm()), mesh_b(init.comm());
  
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
  
        InverseDistanceInterpolation&lt;3&gt; idi (init.comm(),
  					   <I><FONT COLOR="#B22222">/* n_interp_pts = */</FONT></I> 4,
  					   <I><FONT COLOR="#B22222">/* power =        */</FONT></I> 2);
        RadialBasisInterpolation&lt;3&gt; rbi (init.comm());
  
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
  
  	rbi.set_field_variables(field_vars);
  	rbi.get_source_points() = idi.get_source_points();
  	rbi.get_source_vals()   = idi.get_source_vals();
        }
  
        idi.prepare_for_use();
        rbi.prepare_for_use();
  
        {
  	<B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex mutex;
  	MeshlessInterpolationFunction mif(idi, mutex);
  
  	es_b.init();
  	sys_b.project_solution (&amp;mif);
  
  	TecplotIO(mesh_b).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;dest_idi.dat&quot;</FONT></B>,
  						  es_b);
        }
  
        {
  	<B><FONT COLOR="#5F9EA0">Threads</FONT></B>::spin_mutex mutex;
  	MeshlessInterpolationFunction mif(rbi, mutex);
  
  	sys_b.project_solution (&amp;mif);
  
  	TecplotIO(mesh_b).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;dest_rbi.dat&quot;</FONT></B>,
  						  es_b);
        }
      }
  
  
  
    }
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex8'
***************************************************************
* Running Example miscellaneous_ex8:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/radial_basis_interpolation.h, line 79, compiled Apr 19 2013 at 11:51:51 ***
Generating 100 point cloud...done
bounding box is 
(x,y,z)=(    0.34,     0.19,     0.11)
(x,y,z)=(    9.65,     9.96,     9.96)
*** Warning, This code is untested, experimental, or likely to see future API changes: ./include/libmesh/radial_basis_functions.h, line 83, compiled Apr 19 2013 at 11:42:53 ***
bounding box is 
(x,y,z)=(    0.34,     0.19,     0.11)
(x,y,z)=(    9.65,     9.96,     9.96)
r_bbox = 16.7078
rbf(r_bbox/2) = 0.1875
MeshfreeInterpolation
 n_source_points()=400
 n_field_variables()=2
  variables = u v 
Generating 10 point cloud...done
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/solution_transfer/meshfree_interpolation.C, line 283, compiled Apr 19 2013 at 11:42:51 ***
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/solution_transfer/radial_basis_interpolation.C, line 158, compiled Apr 19 2013 at 11:42:53 ***

At target point (x,y,z)=(     0.9,     6.84,     3.76)
 u_interp_idi=2931.24, u_interp_rbi=3006.88, u_exact=2941.14
 v_interp_idi=100.789, v_interp_rbi=101.434, v_exact=100.753

At target point (x,y,z)=(    5.42,     9.36,     1.07)
 u_interp_idi=6246.44, u_interp_rbi=6877.55, u_exact=7836.06
 v_interp_idi=115.044, v_interp_rbi=113.128, v_exact=118.211

At target point (x,y,z)=(    4.45,     7.56,     1.79)
 u_interp_idi=2958.51, u_interp_rbi=3427.22, u_exact=3373.03
 v_interp_idi=70.7865, v_interp_rbi=81.8618, v_exact=82.6914

At target point (x,y,z)=(    4.18,     8.87,     4.12)
 u_interp_idi=4632.05, u_interp_rbi=6651.05, u_exact=7450.19
 v_interp_idi=156.888, v_interp_rbi=162.735, v_exact=166.084

At target point (x,y,z)=(    3.48,     1.72,     6.59)
 u_interp_idi=19625.3, u_interp_rbi=10933.5, u_exact=12479.6
 v_interp_idi=356.174, v_interp_rbi=293.021, v_exact=301.26

At target point (x,y,z)=(    0.09,     3.36,      2.1)
 u_interp_idi=498.469, u_interp_rbi=806.722, u_exact=168.297
 v_interp_idi=39.6305, v_interp_rbi=21.7629, v_exact=20.5587

At target point (x,y,z)=(    3.42,     5.87,     2.06)
 u_interp_idi=971.312, u_interp_rbi=1072.91, u_exact=1264.38
 v_interp_idi=48.0353, v_interp_rbi=53.0049, v_exact=54.8951

At target point (x,y,z)=(    3.01,     7.13,     3.72)
 u_interp_idi=4955.74, u_interp_rbi=3311.65, u_exact=3324.05
 v_interp_idi=162.394, v_interp_rbi=110.872, v_exact=111.376

At target point (x,y,z)=(    3.21,     2.55,     8.19)
 u_interp_idi=26700.4, u_interp_rbi=37098, u_exact=36923.8
 v_interp_idi=466.577, v_interp_rbi=565.731, v_exact=566.16

At target point (x,y,z)=(    5.99,     7.21,     9.04)
 u_interp_idi=53370.1, u_interp_rbi=62518.5, u_exact=63290.2
 v_interp_idi=755.935, v_interp_rbi=822.984, v_exact=826.627
bounding box is 
(x,y,z)=(      -1,        0,        0)
(x,y,z)=(       0,        1,        1)
bounding box is 
(x,y,z)=(      -1,        0,        0)
(x,y,z)=(       0,        1,        1)
r_bbox = 1.73205
rbf(r_bbox/2) = 0.1875

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:52:17 2013                                                                          |
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
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=17.3246, Active time=17.2494                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 2         0.0080      0.004017    0.0091      0.004570    0.05     0.05     |
|   create_dof_constraints()     2         0.0033      0.001668    0.0033      0.001668    0.02     0.02     |
|   distribute_dofs()            2         0.0157      0.007827    0.0465      0.023253    0.09     0.27     |
|   dof_indices()                7528      0.0272      0.000004    0.0272      0.000004    0.16     0.16     |
|   prepare_send_list()          2         0.0000      0.000013    0.0000      0.000013    0.00     0.00     |
|   reinit()                     2         0.0269      0.013447    0.0269      0.013447    0.16     0.16     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      3         0.0086      0.002855    0.0250      0.008318    0.05     0.14     |
|                                                                                                            |
| InverseDistanceInterpolation<>                                                                             |
|   construct_kd_tree()          4         0.0064      0.001610    0.0064      0.001610    0.04     0.04     |
|   interpolate_field_data()     4219      0.0128      0.000003    0.0128      0.000003    0.07     0.07     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0209      0.010456    0.0345      0.017260    0.12     0.20     |
|   read()                       2         0.0407      0.020368    0.0407      0.020368    0.24     0.24     |
|   renumber_nodes_and_elem()    4         0.0020      0.000507    0.0020      0.000507    0.01     0.01     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   broadcast()                  2         0.0106      0.005299    0.0252      0.012591    0.06     0.15     |
|   compute_hilbert_indices()    4         0.0723      0.018084    0.0723      0.018084    0.42     0.42     |
|   find_global_indices()        4         0.0092      0.002311    0.1095      0.027365    0.05     0.63     |
|   parallel_sort()              4         0.0017      0.000415    0.0267      0.006686    0.01     0.16     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     3         0.0001      0.000040    0.0997      0.033224    0.00     0.58     |
|                                                                                                            |
| MeshfreeInterpolation                                                                                      |
|   gather_remote_data()         4         0.0003      0.000085    0.0009      0.000237    0.00     0.01     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0937      0.046840    0.1419      0.070953    0.54     0.82     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  20        0.0030      0.000150    0.0034      0.000170    0.02     0.02     |
|   broadcast()                  48        0.0044      0.000091    0.0043      0.000090    0.03     0.03     |
|   max(scalar)                  243       0.0045      0.000018    0.0045      0.000018    0.03     0.03     |
|   max(vector)                  56        0.0009      0.000016    0.0032      0.000058    0.01     0.02     |
|   min(bool)                    289       0.0037      0.000013    0.0037      0.000013    0.02     0.02     |
|   min(scalar)                  237       0.9160      0.003865    0.9160      0.003865    5.31     5.31     |
|   min(vector)                  56        0.0010      0.000018    0.0058      0.000104    0.01     0.03     |
|   probe()                      72        0.0024      0.000033    0.0024      0.000033    0.01     0.01     |
|   receive()                    72        0.0004      0.000005    0.0028      0.000038    0.00     0.02     |
|   send()                       72        0.0003      0.000004    0.0003      0.000004    0.00     0.00     |
|   send_receive()               68        0.0003      0.000004    0.0028      0.000041    0.00     0.02     |
|   sum()                        24        0.0254      0.001059    0.0284      0.001181    0.15     0.16     |
|                                                                                                            |
| Parallel::Request                                                                                          |
|   wait()                       72        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0054      0.002697    0.0215      0.010748    0.03     0.12     |
|   set_parent_processor_ids()   2         0.0022      0.001107    0.0022      0.001107    0.01     0.01     |
|                                                                                                            |
| RadialBasisInterpolation<>                                                                                 |
|   interpolate_field_data()     4219      0.2936      0.000070    0.2936      0.000070    1.70     1.70     |
|   prepare_for_use()            2         15.3048     7.652407    15.3048     7.652407    88.73    88.73    |
|                                                                                                            |
| System                                                                                                     |
|   project_vector()             3         0.2461      0.082041    0.5710      0.190338    1.43     3.31     |
|                                                                                                            |
| TecplotIO                                                                                                  |
|   write_nodal_data()           3         0.0744      0.024791    0.0744      0.024791    0.43     0.43     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        17355     17.2494                                         100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex8:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex8'
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
