<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("vector_fe_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file laplace_exact_solution.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        
        using namespace libMesh;
        
        #ifndef __laplace_exact_solution_h__
        #define __laplace_exact_solution_h__
        
        class LaplaceExactSolution
        {
        public:
          LaplaceExactSolution(){}
        
          ~LaplaceExactSolution(){}
        
          Real operator()( unsigned int component,
        		   Real x, Real y, Real z = 0.0)
          {
            const Real hp = 0.5*pi;
        
            switch(component)
            {
            case 0:
              return cos(hp*x)*sin(hp*y)*cos(hp*z);
        
            case 1:
              return sin(hp*x)*cos(hp*y)*cos(hp*z);
        
            case 2:
              return sin(hp*x)*cos(hp*y)*sin(hp*z);
        
            default:
              libmesh_error();
            }
          }
        };
        
        
        class LaplaceExactGradient
        {
        public:
          LaplaceExactGradient(){}
        
          ~LaplaceExactGradient(){}
        
          RealGradient operator()( unsigned int component,
        			   Real x, Real y, Real z = 0.0)
          {
            const Real hp = 0.5*pi;
        
            switch(component)
            {
            case 0:
              return RealGradient( -hp*sin(hp*x)*sin(hp*y)*cos(hp*z),
        			   cos(hp*x)*(hp)*cos(hp*y)*cos(hp*z),
        			   cos(hp*x)*sin(hp*y)*(-hp)*sin(hp*z) );
        
            case 1:
              return RealGradient( hp*cos(hp*x)*cos(hp*y)*cos(hp*z),
        			   sin(hp*x)*(-hp)*sin(hp*y)*cos(hp*z),
        			   sin(hp*x)*cos(hp*y)*(-hp)*sin(hp*z) );
        
            case 2:
              return RealGradient( hp*cos(hp*x)*cos(hp*y)*sin(hp*z),
        			   sin(hp*x)*(-hp)*sin(hp*y)*sin(hp*z),
        			   sin(hp*x)*cos(hp*y)*(hp)*cos(hp*z) );
        
            default:
              libmesh_error();
            }
          }
        };
        
        #endif // __laplace_exact_solution_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file laplace_system.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fem_system.h"
        #include "libmesh/vector_value.h"
        #include "libmesh/tensor_value.h"
        #include "libmesh/dirichlet_boundaries.h"
        
        #include "solution_function.h"
        
        using namespace libMesh;
        
        #ifndef __laplace_system_h__
        #define __laplace_system_h__
        
</pre>
</div>
<div class = "comment">
FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
but we must specify element residuals
</div>

<div class ="fragment">
<pre>
        class LaplaceSystem : public FEMSystem
        {
        
        public:
        
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          LaplaceSystem( EquationSystems& es,
        		 const std::string& name,
        		 const unsigned int number);
        
</pre>
</div>
<div class = "comment">
System initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_data ();
        
</pre>
</div>
<div class = "comment">
Context initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_context(DiffContext &context);
        
</pre>
</div>
<div class = "comment">
Element residual and jacobian calculations
Time dependent parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_time_derivative (bool request_jacobian,
                                                DiffContext& context);
        
</pre>
</div>
<div class = "comment">
Constraint part
</div>

<div class ="fragment">
<pre>
          /*
          virtual bool side_constraint (bool request_jacobian,
        				DiffContext& context);
          */
        
        protected:
</pre>
</div>
<div class = "comment">
Indices for each variable;
</div>

<div class ="fragment">
<pre>
          unsigned int u_var;
        
          void init_dirichlet_bcs();
        
</pre>
</div>
<div class = "comment">
Returns the value of a forcing function at point p.
</div>

<div class ="fragment">
<pre>
          RealGradient forcing(const Point& p);
        
          LaplaceExactSolution exact_solution;
        };
        
        #endif //__laplace_system_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file solution_function.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/function_base.h"
        #include "laplace_exact_solution.h"
        
        using namespace libMesh;
        
        #ifndef __solution_function_h__
        #define __solution_function_h__
        
        class SolutionFunction : public FunctionBase&lt;Number&gt;
        {
        public:
        
          SolutionFunction( const unsigned int u_var )
          : _u_var(u_var) {}
          ~SolutionFunction( ){}
        
          virtual Number operator() (const Point&, const Real = 0)
            { libmesh_not_implemented(); }
        
          virtual void operator() (const Point& p,
                                   const Real,
                                   DenseVector&lt;Number&gt;& output)
          {
            output.zero();
            const Real x=p(0), y=p(1), z=p(2);
</pre>
</div>
<div class = "comment">
libMesh assumes each component of the vector-valued variable is stored
contiguously.
</div>

<div class ="fragment">
<pre>
            output(_u_var)   = soln( 0, x, y, z );
            output(_u_var+1) = soln( 1, x, y, z );
            output(_u_var+2) = soln( 2, x, y, z );
          }
        
          virtual Number component( unsigned int component_in, const Point& p,
        			    const Real )
          {
            const Real x=p(0), y=p(1), z=p(2);
            return soln( component_in, x, y, z );
          }
        
          virtual AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() const
          { return AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (new SolutionFunction(_u_var)); }
        
        private:
        
          const unsigned int _u_var;
          LaplaceExactSolution soln;
        };
        
</pre>
</div>
<div class = "comment">
FIXME: PB: We ought to be able to merge the above class with this one
through templating, but I'm being lazy.
</div>

<div class ="fragment">
<pre>
        class SolutionGradient : public FunctionBase&lt;Gradient&gt;
        {
        public:
        
          SolutionGradient( const unsigned int u_var )
          : _u_var(u_var) {}
          ~SolutionGradient( ){}
        
          virtual Gradient operator() (const Point&, const Real = 0)
            { libmesh_not_implemented(); }
        
          virtual void operator() (const Point& p,
                                   const Real,
                                   DenseVector&lt;Gradient&gt;& output)
          {
            output.zero();
            const Real x=p(0), y=p(1), z=p(2);
            output(_u_var)   = soln( 0, x, y, z );
            output(_u_var+1) = soln( 1, x, y, z );
            output(_u_var+2) = soln( 2, x, y, z );
          }
        
          virtual Gradient component( unsigned int component_in, const Point& p,
        			    const Real )
          {
            const Real x=p(0), y=p(1), z=p(2);
            return soln( component_in, x, y, z );
          }
        
          virtual AutoPtr&lt;FunctionBase&lt;Gradient&gt; &gt; clone() const
          { return AutoPtr&lt;FunctionBase&lt;Gradient&gt; &gt; (new SolutionGradient(_u_var)); }
        
        private:
        
          const unsigned int _u_var;
          LaplaceExactGradient soln;
        };
        
        #endif // __solution_function_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file laplace_system.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/getpot.h"
        
        #include "laplace_system.h"
        
        #include "libmesh/boundary_info.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/mesh.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/string_to_enum.h"
        
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        LaplaceSystem::LaplaceSystem( EquationSystems& es,
        			      const std::string& name_in,
        			      const unsigned int number_in)
          : FEMSystem(es, name_in, number_in)
        {
          return;
        }
        
        void LaplaceSystem::init_data ()
        {
</pre>
</div>
<div class = "comment">
Check the input file for Reynolds number, application type,
approximation type
</div>

<div class ="fragment">
<pre>
          GetPot infile("laplace.in");
        
</pre>
</div>
<div class = "comment">
Add the solution variable
</div>

<div class ="fragment">
<pre>
          u_var = this-&gt;add_variable ("u", FIRST, LAGRANGE_VEC);
        
          this-&gt;time_evolving(u_var);
        
</pre>
</div>
<div class = "comment">
Useful debugging options
Set verify_analytic_jacobians to 1e-6 to use
</div>

<div class ="fragment">
<pre>
          this-&gt;verify_analytic_jacobians = infile("verify_analytic_jacobians", 0.);
          this-&gt;print_jacobians = infile("print_jacobians", false);
          this-&gt;print_element_jacobians = infile("print_element_jacobians", false);
        
          this-&gt;init_dirichlet_bcs();
        
</pre>
</div>
<div class = "comment">
Do the parent's initialization after variables and boundary constraints are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        }
        
        void LaplaceSystem::init_dirichlet_bcs()
        {
          const boundary_id_type all_ids[6] = {0,1,2,3,4,5};
          std::set&lt;boundary_id_type&gt; boundary_ids( all_ids, all_ids+6 );
        
          std::vector&lt;unsigned int&gt; vars;
          vars.push_back( u_var );
        
</pre>
</div>
<div class = "comment">
Note that for vector-valued variables, it is assumed each component is stored contiguously.
For 2-D elements in 3-D space, only two components should be returned.
</div>

<div class ="fragment">
<pre>
          SolutionFunction func( u_var );
        
          this-&gt;get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( boundary_ids, vars, &func ) );
          return;
        }
        
        void LaplaceSystem::init_context(DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Get finite element object
</div>

<div class ="fragment">
<pre>
          FEGenericBase&lt;RealGradient&gt;* fe;
          c.get_element_fe&lt;RealGradient&gt;( u_var, fe );
        
</pre>
</div>
<div class = "comment">
We should prerequest all the data
we will need to build the linear system.
</div>

<div class ="fragment">
<pre>
          fe-&gt;get_JxW();
          fe-&gt;get_phi();
          fe-&gt;get_dphi();
          fe-&gt;get_xyz();
        
          FEGenericBase&lt;RealGradient&gt;* side_fe;
          c.get_side_fe&lt;RealGradient&gt;( u_var, side_fe );
        
          side_fe-&gt;get_JxW();
          side_fe-&gt;get_phi();
        }
        
        
        bool LaplaceSystem::element_time_derivative (bool request_jacobian,
                                                    DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Get finite element object
</div>

<div class ="fragment">
<pre>
          FEGenericBase&lt;RealGradient&gt;* fe = NULL;
          c.get_element_fe&lt;RealGradient&gt;( u_var, fe );
        
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& phi = fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealTensor&gt; &gt;& grad_phi = fe-&gt;get_dphi();
        
          const std::vector&lt;Point&gt;& qpoint = fe-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
          DenseSubMatrix&lt;Number&gt; &Kuu = *c.elem_subjacobians[u_var][u_var];
        
          DenseSubVector&lt;Number&gt; &Fu = *c.elem_subresiduals[u_var];
        
</pre>
</div>
<div class = "comment">
Now we will build the element Jacobian and residual.
Constructing the residual requires the solution and its
gradient from the previous timestep.  This must be
calculated at each quadrature point by summing the
solution degree-of-freedom values by the appropriate
weight functions.
</div>

<div class ="fragment">
<pre>
          const unsigned int n_qpoints = (c.get_element_qrule())-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              Tensor grad_u;
        
              c.interior_gradient( u_var, qp, grad_u );
        
        
</pre>
</div>
<div class = "comment">
Value of the forcing function at this quadrature point
</div>

<div class ="fragment">
<pre>
              RealGradient f = this-&gt;forcing(qpoint[qp]);
        
</pre>
</div>
<div class = "comment">
First, an i-loop over the velocity degrees of freedom.
We know that n_u_dofs == n_v_dofs so we can compute contributions
for both at the same time.
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_u_dofs; i++)
                {
                  Fu(i) += ( grad_u.contract(grad_phi[i][qp]) - f*phi[i][qp] )*JxW[qp];
        
                  if (request_jacobian)
                    {
</pre>
</div>
<div class = "comment">
Matrix contributions for the uu and vv couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j != n_u_dofs; j++)
                        {
                          Kuu(i,j) += grad_phi[j][qp].contract(grad_phi[i][qp])*JxW[qp];
        
        		}
        	    }
                }
            } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
        /*
        bool LaplaceSystem::side_constraint (bool request_jacobian,
        				     DiffContext &context)
        {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Get finite element object
</div>

<div class ="fragment">
<pre>
          FEGenericBase&lt;RealGradient&gt;* side_fe = NULL;
          c.get_side_fe&lt;RealGradient&gt;( u_var, side_fe );
        
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = side_fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& phi = side_fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
        
          const std::vector&lt;Point&gt;& qpoint = side_fe-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The penalty value.  \frac{1}{\epsilon}
in the discussion above.
</div>

<div class ="fragment">
<pre>
          const Real penalty = 1.e10;
        
          DenseSubMatrix&lt;Number&gt; &Kuu = *c.elem_subjacobians[u_var][u_var];
          DenseSubVector&lt;Number&gt; &Fu = *c.elem_subresiduals[u_var];
        
          const unsigned int n_qpoints = (c.get_side_qrule())-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
              Gradient u;
              c.side_value( u_var, qp, u );
        
              Gradient u_exact( this-&gt;exact_solution( 0, qpoint[qp](0), qpoint[qp](1) ),
        			this-&gt;exact_solution( 1, qpoint[qp](0), qpoint[qp](1) ));
        
              for (unsigned int i=0; i != n_u_dofs; i++)
        	{
        	  Fu(i) += penalty*(u - u_exact)*phi[i][qp]*JxW[qp];
        
        	  if (request_jacobian)
        	    {
        	      for (unsigned int j=0; j != n_u_dofs; j++)
        		Kuu(i,j) += penalty*phi[j][qp]*phi[i][qp]*JxW[qp];
        	    }
        	}
            }
        
          return request_jacobian;
        }
        */
        
        RealGradient LaplaceSystem::forcing( const Point& p )
        {
          Real x = p(0); Real y = p(1); Real z = p(2);
        
          const Real eps = 1.e-3;
        
          const Real fx = -(exact_solution(0,x,y,z-eps) +
        		    exact_solution(0,x,y,z+eps) +
        		    exact_solution(0,x,y-eps,z) +
        		    exact_solution(0,x,y+eps,z) +
        		    exact_solution(0,x-eps,y,z) +
        		    exact_solution(0,x+eps,y,z) -
        		    6.*exact_solution(0,x,y,z))/eps/eps;
        
          const Real fy = -(exact_solution(1,x,y,z-eps) +
        		    exact_solution(1,x,y,z+eps) +
        		    exact_solution(1,x,y-eps,z) +
        		    exact_solution(1,x,y+eps,z) +
        		    exact_solution(1,x-eps,y,z) +
        		    exact_solution(1,x+eps,y,z) -
        		    6.*exact_solution(1,x,y,z))/eps/eps;
        
          const Real fz = -(exact_solution(2,x,y,z-eps) +
        		    exact_solution(2,x,y,z+eps) +
        		    exact_solution(2,x,y-eps,z) +
        		    exact_solution(2,x,y+eps,z) +
        		    exact_solution(2,x-eps,y,z) +
        		    exact_solution(2,x+eps,y,z) -
        		    6.*exact_solution(2,x,y,z))/eps/eps;
        
          return RealGradient( fx, fy, fz );
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file vector_fe_ex2.C with comments: </h1> 
<div class = "comment">
<h1>FEMSystem Example 1 - Unsteady Navier-Stokes Equations with
FEMSystem</h1>

<br><br>This example shows how the transient nonlinear problem from
example 13 can be solved using the
DifferentiableSystem class framework


<br><br>Basic include files
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        #include "libmesh/getpot.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exact_solution.h"
        #include "libmesh/ucd_io.h"
        
</pre>
</div>
<div class = "comment">
The systems and solvers we may use
</div>

<div class ="fragment">
<pre>
        #include "laplace_system.h"
        #include "libmesh/diff_solver.h"
        #include "libmesh/steady_solver.h"
        
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
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Parse the input file
</div>

<div class ="fragment">
<pre>
          GetPot infile("vector_fe_ex2.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
          const unsigned int grid_size = infile( "grid_size", 2 );
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "2D/3D support");
        
</pre>
</div>
<div class = "comment">
Create a mesh, with dimension to be overridden later, on the
default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_cube (mesh,
        				     grid_size,
        				     grid_size,
        				     grid_size,
        				     -1., 1.,
        				     -1., 1.,
        				     -1., 1.,
        				     HEX8);
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
          mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
</pre>
</div>
<div class = "comment">
Declare the system "Navier-Stokes" and its variables.
</div>

<div class ="fragment">
<pre>
          LaplaceSystem & system =
            equation_systems.add_system&lt;LaplaceSystem&gt; ("Laplace");
        
</pre>
</div>
<div class = "comment">
This example only implements the steady-state problem
</div>

<div class ="fragment">
<pre>
          system.time_solver =
            AutoPtr&lt;TimeSolver&gt;(new SteadySolver(system));
        
</pre>
</div>
<div class = "comment">
Initialize the system
</div>

<div class ="fragment">
<pre>
          equation_systems.init();
        
</pre>
</div>
<div class = "comment">
And the nonlinear solver options
</div>

<div class ="fragment">
<pre>
          DiffSolver &solver = *(system.time_solver-&gt;diff_solver().get());
          solver.quiet = infile("solver_quiet", true);
          solver.verbose = !solver.quiet;
          solver.max_nonlinear_iterations =
            infile("max_nonlinear_iterations", 15);
          solver.relative_step_tolerance =
            infile("relative_step_tolerance", 1.e-3);
          solver.relative_residual_tolerance =
            infile("relative_residual_tolerance", 0.0);
          solver.absolute_residual_tolerance =
            infile("absolute_residual_tolerance", 0.0);
        
</pre>
</div>
<div class = "comment">
And the linear solver options
</div>

<div class ="fragment">
<pre>
          solver.max_linear_iterations =
            infile("max_linear_iterations", 50000);
          solver.initial_linear_tolerance =
            infile("initial_linear_tolerance", 1.e-3);
        
</pre>
</div>
<div class = "comment">
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
          system.solve();
        
          ExactSolution exact_sol( equation_systems );
        
          std::vector&lt;FunctionBase&lt;Number&gt;* &gt; sols;
          std::vector&lt;FunctionBase&lt;Gradient&gt;* &gt; grads;
        
          sols.push_back( new SolutionFunction(system.variable_number("u")) );
          grads.push_back( new SolutionGradient(system.variable_number("u")) );
        
          exact_sol.attach_exact_values(sols);
          exact_sol.attach_exact_derivs(grads);
        
</pre>
</div>
<div class = "comment">
Use higher quadrature order for more accurate error results
</div>

<div class ="fragment">
<pre>
          int extra_error_quadrature = infile("extra_error_quadrature",2);
          exact_sol.extra_quadrature_order(extra_error_quadrature);
        
</pre>
</div>
<div class = "comment">
Compute the error.
</div>

<div class ="fragment">
<pre>
          exact_sol.compute_error("Laplace", "u");
        
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "L2-Error is: "
        	    &lt;&lt; exact_sol.l2_error("Laplace", "u")
        	    &lt;&lt; std::endl;
          std::cout &lt;&lt; "H1-Error is: "
        	    &lt;&lt; exact_sol.h1_error("Laplace", "u")
        	    &lt;&lt; std::endl;
        
        #ifdef LIBMESH_HAVE_EXODUS_API
        
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
          ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
        
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
          UCDIO(mesh).write_equation_systems("out.inp", equation_systems);
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return 0;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file laplace_exact_solution.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  #ifndef __laplace_exact_solution_h__
  #define __laplace_exact_solution_h__
  
  <B><FONT COLOR="#228B22">class</FONT></B> LaplaceExactSolution
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    LaplaceExactSolution(){}
  
    ~LaplaceExactSolution(){}
  
    Real <B><FONT COLOR="#A020F0">operator</FONT></B>()( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component,
  		   Real x, Real y, Real z = 0.0)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real hp = 0.5*pi;
  
      <B><FONT COLOR="#A020F0">switch</FONT></B>(component)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> cos(hp*x)*sin(hp*y)*cos(hp*z);
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> sin(hp*x)*cos(hp*y)*cos(hp*z);
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> sin(hp*x)*cos(hp*y)*sin(hp*z);
  
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        libmesh_error();
      }
    }
  };
  
  
  <B><FONT COLOR="#228B22">class</FONT></B> LaplaceExactGradient
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    LaplaceExactGradient(){}
  
    ~LaplaceExactGradient(){}
  
    RealGradient <B><FONT COLOR="#A020F0">operator</FONT></B>()( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component,
  			   Real x, Real y, Real z = 0.0)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real hp = 0.5*pi;
  
      <B><FONT COLOR="#A020F0">switch</FONT></B>(component)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">0</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> RealGradient( -hp*sin(hp*x)*sin(hp*y)*cos(hp*z),
  			   cos(hp*x)*(hp)*cos(hp*y)*cos(hp*z),
  			   cos(hp*x)*sin(hp*y)*(-hp)*sin(hp*z) );
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> RealGradient( hp*cos(hp*x)*cos(hp*y)*cos(hp*z),
  			   sin(hp*x)*(-hp)*sin(hp*y)*cos(hp*z),
  			   sin(hp*x)*cos(hp*y)*(-hp)*sin(hp*z) );
  
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
        <B><FONT COLOR="#A020F0">return</FONT></B> RealGradient( hp*cos(hp*x)*cos(hp*y)*sin(hp*z),
  			   sin(hp*x)*(-hp)*sin(hp*y)*sin(hp*z),
  			   sin(hp*x)*cos(hp*y)*(hp)*cos(hp*z) );
  
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        libmesh_error();
      }
    }
  };
  
  #endif <I><FONT COLOR="#B22222">// __laplace_exact_solution_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file laplace_system.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tensor_value.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;solution_function.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  #ifndef __laplace_system_h__
  #define __laplace_system_h__
  
  <B><FONT COLOR="#228B22">class</FONT></B> LaplaceSystem : <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem
  {
  
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    LaplaceSystem( EquationSystems&amp; es,
  		 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
  		 <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data ();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                          DiffContext&amp; context);
  
    <I><FONT COLOR="#B22222">/*
    virtual bool side_constraint (bool request_jacobian,
  				DiffContext&amp; context);
    */</FONT></I>
  
  <B><FONT COLOR="#228B22">protected</FONT></B>:
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var;
  
    <B><FONT COLOR="#228B22">void</FONT></B> init_dirichlet_bcs();
  
    RealGradient forcing(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p);
  
    LaplaceExactSolution exact_solution;
  };
  
  #endif <I><FONT COLOR="#B22222">//__laplace_system_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file solution_function.h without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/function_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;laplace_exact_solution.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  #ifndef __solution_function_h__
  #define __solution_function_h__
  
  <B><FONT COLOR="#228B22">class</FONT></B> SolutionFunction : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Number&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SolutionFunction( <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var )
    : _u_var(u_var) {}
    ~SolutionFunction( ){}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Real = 0)
      { libmesh_not_implemented(); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                             <B><FONT COLOR="#228B22">const</FONT></B> Real,
                             DenseVector&lt;Number&gt;&amp; output)
    {
      output.zero();
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      output(_u_var)   = soln( 0, x, y, z );
      output(_u_var+1) = soln( 1, x, y, z );
      output(_u_var+2) = soln( 2, x, y, z );
    }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number component( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component_in, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			    <B><FONT COLOR="#228B22">const</FONT></B> Real )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      <B><FONT COLOR="#A020F0">return</FONT></B> soln( component_in, x, y, z );
    }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; clone() <B><FONT COLOR="#228B22">const</FONT></B>
    { <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Number&gt; &gt; (<B><FONT COLOR="#A020F0">new</FONT></B> SolutionFunction(_u_var)); }
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _u_var;
    LaplaceExactSolution soln;
  };
  
  <B><FONT COLOR="#228B22">class</FONT></B> SolutionGradient : <B><FONT COLOR="#228B22">public</FONT></B> FunctionBase&lt;Gradient&gt;
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
  
    SolutionGradient( <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var )
    : _u_var(u_var) {}
    ~SolutionGradient( ){}
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Gradient <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Real = 0)
      { libmesh_not_implemented(); }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> <B><FONT COLOR="#A020F0">operator</FONT></B>() (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                             <B><FONT COLOR="#228B22">const</FONT></B> Real,
                             DenseVector&lt;Gradient&gt;&amp; output)
    {
      output.zero();
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      output(_u_var)   = soln( 0, x, y, z );
      output(_u_var+1) = soln( 1, x, y, z );
      output(_u_var+2) = soln( 2, x, y, z );
    }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Gradient component( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component_in, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			    <B><FONT COLOR="#228B22">const</FONT></B> Real )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      <B><FONT COLOR="#A020F0">return</FONT></B> soln( component_in, x, y, z );
    }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> AutoPtr&lt;FunctionBase&lt;Gradient&gt; &gt; clone() <B><FONT COLOR="#228B22">const</FONT></B>
    { <B><FONT COLOR="#A020F0">return</FONT></B> AutoPtr&lt;FunctionBase&lt;Gradient&gt; &gt; (<B><FONT COLOR="#A020F0">new</FONT></B> SolutionGradient(_u_var)); }
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> _u_var;
    LaplaceExactGradient soln;
  };
  
  #endif <I><FONT COLOR="#B22222">// __solution_function_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file laplace_system.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;laplace_system.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  
  
  using namespace libMesh;
  
  <B><FONT COLOR="#5F9EA0">LaplaceSystem</FONT></B>::LaplaceSystem( EquationSystems&amp; es,
  			      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
  			      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in)
    : FEMSystem(es, name_in, number_in)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::init_data ()
  {
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;laplace.in&quot;</FONT></B>);
  
    u_var = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST, LAGRANGE_VEC);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(u_var);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;verify_analytic_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;verify_analytic_jacobians&quot;</FONT></B>, 0.);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_jacobians&quot;</FONT></B>, false);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;print_element_jacobians = infile(<B><FONT COLOR="#BC8F8F">&quot;print_element_jacobians&quot;</FONT></B>, false);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;init_dirichlet_bcs();
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::init_dirichlet_bcs()
  {
    <B><FONT COLOR="#228B22">const</FONT></B> boundary_id_type all_ids[6] = {0,1,2,3,4,5};
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; boundary_ids( all_ids, all_ids+6 );
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; vars;
    vars.push_back( u_var );
  
    SolutionFunction func( u_var );
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_dof_map().add_dirichlet_boundary( libMesh::DirichletBoundary( boundary_ids, vars, &amp;func ) );
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> LaplaceSystem::init_context(DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEGenericBase&lt;RealGradient&gt;* fe;
    c.get_element_fe&lt;RealGradient&gt;( u_var, fe );
  
    fe-&gt;get_JxW();
    fe-&gt;get_phi();
    fe-&gt;get_dphi();
    fe-&gt;get_xyz();
  
    FEGenericBase&lt;RealGradient&gt;* side_fe;
    c.get_side_fe&lt;RealGradient&gt;( u_var, side_fe );
  
    side_fe-&gt;get_JxW();
    side_fe-&gt;get_phi();
  }
  
  
  <B><FONT COLOR="#228B22">bool</FONT></B> LaplaceSystem::element_time_derivative (<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
                                              DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEGenericBase&lt;RealGradient&gt;* fe = NULL;
    c.get_element_fe&lt;RealGradient&gt;( u_var, fe );
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealTensor&gt; &gt;&amp; grad_phi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qpoint = fe-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[u_var].size();
  
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = *c.elem_subjacobians[u_var][u_var];
  
    DenseSubVector&lt;Number&gt; &amp;Fu = *c.elem_subresiduals[u_var];
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = (c.get_element_qrule())-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Tensor grad_u;
  
        c.interior_gradient( u_var, qp, grad_u );
  
  
        RealGradient f = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;forcing(qpoint[qp]);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i != n_u_dofs; i++)
          {
            Fu(i) += ( grad_u.contract(grad_phi[i][qp]) - f*phi[i][qp] )*JxW[qp];
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian)
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j != n_u_dofs; j++)
                  {
                    Kuu(i,j) += grad_phi[j][qp].contract(grad_phi[i][qp])*JxW[qp];
  
  		}
  	    }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  <I><FONT COLOR="#B22222">/*
  bool LaplaceSystem::side_constraint (bool request_jacobian,
  				     DiffContext &amp;context)
  {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    FEGenericBase&lt;RealGradient&gt;* side_fe = NULL;
    c.get_side_fe&lt;RealGradient&gt;( u_var, side_fe );
  
  
    const std::vector&lt;Real&gt; &amp;JxW = side_fe-&gt;get_JxW();
  
    const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; phi = side_fe-&gt;get_phi();
  
    const unsigned int n_u_dofs = c.dof_indices_var[u_var].size();
  
    const std::vector&lt;Point&gt;&amp; qpoint = side_fe-&gt;get_xyz();
  
    const Real penalty = 1.e10;
  
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = *c.elem_subjacobians[u_var][u_var];
    DenseSubVector&lt;Number&gt; &amp;Fu = *c.elem_subresiduals[u_var];
  
    const unsigned int n_qpoints = (c.get_side_qrule())-&gt;n_points();
  
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Gradient u;
        c.side_value( u_var, qp, u );
  
        Gradient u_exact( this-&gt;exact_solution( 0, qpoint[qp](0), qpoint[qp](1) ),
  			this-&gt;exact_solution( 1, qpoint[qp](0), qpoint[qp](1) ));
  
        for (unsigned int i=0; i != n_u_dofs; i++)
  	{
  	  Fu(i) += penalty*(u - u_exact)*phi[i][qp]*JxW[qp];
  
  	  if (request_jacobian)
  	    {
  	      for (unsigned int j=0; j != n_u_dofs; j++)
  		Kuu(i,j) += penalty*phi[j][qp]*phi[i][qp]*JxW[qp];
  	    }
  	}
      }
  
    return request_jacobian;
  }
  */</FONT></I>
  
  RealGradient LaplaceSystem::forcing( <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p )
  {
    Real x = p(0); Real y = p(1); Real z = p(2);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real eps = 1.e-3;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real fx = -(exact_solution(0,x,y,z-eps) +
  		    exact_solution(0,x,y,z+eps) +
  		    exact_solution(0,x,y-eps,z) +
  		    exact_solution(0,x,y+eps,z) +
  		    exact_solution(0,x-eps,y,z) +
  		    exact_solution(0,x+eps,y,z) -
  		    6.*exact_solution(0,x,y,z))/eps/eps;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real fy = -(exact_solution(1,x,y,z-eps) +
  		    exact_solution(1,x,y,z+eps) +
  		    exact_solution(1,x,y-eps,z) +
  		    exact_solution(1,x,y+eps,z) +
  		    exact_solution(1,x-eps,y,z) +
  		    exact_solution(1,x+eps,y,z) -
  		    6.*exact_solution(1,x,y,z))/eps/eps;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real fz = -(exact_solution(2,x,y,z-eps) +
  		    exact_solution(2,x,y,z+eps) +
  		    exact_solution(2,x,y-eps,z) +
  		    exact_solution(2,x,y+eps,z) +
  		    exact_solution(2,x-eps,y,z) +
  		    exact_solution(2,x+eps,y,z) -
  		    6.*exact_solution(2,x,y,z))/eps/eps;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> RealGradient( fx, fy, fz );
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file vector_fe_ex2.C without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exact_solution.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/ucd_io.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;laplace_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    GetPot infile(<B><FONT COLOR="#BC8F8F">&quot;vector_fe_ex2.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> grid_size = infile( <B><FONT COLOR="#BC8F8F">&quot;grid_size&quot;</FONT></B>, 2 );
  
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
    Mesh mesh(init.comm());
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
  				     grid_size,
  				     grid_size,
  				     grid_size,
  				     -1., 1.,
  				     -1., 1.,
  				     -1., 1.,
  				     HEX8);
  
    mesh.print_info();
  
    EquationSystems equation_systems (mesh);
  
    LaplaceSystem &amp; system =
      equation_systems.add_system&lt;LaplaceSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
    system.time_solver =
      AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(system));
  
    equation_systems.init();
  
    DiffSolver &amp;solver = *(system.time_solver-&gt;diff_solver().get());
    solver.quiet = infile(<B><FONT COLOR="#BC8F8F">&quot;solver_quiet&quot;</FONT></B>, true);
    solver.verbose = !solver.quiet;
    solver.max_nonlinear_iterations =
      infile(<B><FONT COLOR="#BC8F8F">&quot;max_nonlinear_iterations&quot;</FONT></B>, 15);
    solver.relative_step_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;relative_step_tolerance&quot;</FONT></B>, 1.e-3);
    solver.relative_residual_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;relative_residual_tolerance&quot;</FONT></B>, 0.0);
    solver.absolute_residual_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;absolute_residual_tolerance&quot;</FONT></B>, 0.0);
  
    solver.max_linear_iterations =
      infile(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 50000);
    solver.initial_linear_tolerance =
      infile(<B><FONT COLOR="#BC8F8F">&quot;initial_linear_tolerance&quot;</FONT></B>, 1.e-3);
  
    equation_systems.print_info();
  
    system.solve();
  
    ExactSolution exact_sol( equation_systems );
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FunctionBase&lt;Number&gt;* &gt; sols;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;FunctionBase&lt;Gradient&gt;* &gt; grads;
  
    sols.push_back( <B><FONT COLOR="#A020F0">new</FONT></B> SolutionFunction(system.variable_number(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)) );
    grads.push_back( <B><FONT COLOR="#A020F0">new</FONT></B> SolutionGradient(system.variable_number(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)) );
  
    exact_sol.attach_exact_values(sols);
    exact_sol.attach_exact_derivs(grads);
  
    <B><FONT COLOR="#228B22">int</FONT></B> extra_error_quadrature = infile(<B><FONT COLOR="#BC8F8F">&quot;extra_error_quadrature&quot;</FONT></B>,2);
    exact_sol.extra_quadrature_order(extra_error_quadrature);
  
    exact_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Error is: &quot;</FONT></B>
  	    &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  	    &lt;&lt; std::endl;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1-Error is: &quot;</FONT></B>
  	    &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  	    &lt;&lt; std::endl;
  
  #ifdef LIBMESH_HAVE_EXODUS_API
  
    ExodusII_IO(mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out.e&quot;</FONT></B>, equation_systems);
  
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    UCDIO(mesh).write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out.inp&quot;</FONT></B>, equation_systems);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/vector_fe/vector_fe_ex2'
***************************************************************
* Running Example vector_fe_ex2:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=216
    n_local_nodes()=78
  n_elem()=125
    n_local_elem()=32
    n_active_elem()=125
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Laplace"
    Type "Implicit"
    Variables="u" 
    Finite Element Types="LAGRANGE_VEC", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=648
    n_local_dofs()=234
    n_constrained_dofs()=456
    n_local_constrained_dofs()=156
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 50.1111
      Average Off-Processor Bandwidth <= 16.1111
      Maximum  On-Processor Bandwidth <= 108
      Maximum Off-Processor Bandwidth <= 99
    DofMap Constraints
      Number of DoF Constraints = 456
      Number of Heterogenous Constraints= 456
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

Assembling the System
Nonlinear Residual: 4.00845
Linear solve starting, tolerance 1e-12
Linear solve finished, step 15, residual 1.93334e-12
Trying full Newton step
  Current Residual: 2.01512e-12
  Nonlinear solver converged, step 0, residual reduction 5.02718e-13 < 1e-12
  Nonlinear solver relative step size 0.723571 > 1e-06
L2-Error is: 0.142039
H1-Error is: 0.902143

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:58:48 2013                                                                          |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.358299, Active time=0.227149                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0007      0.000655    0.0021      0.002142    0.29     0.94     |
|   build_constraint_matrix()        64        0.0006      0.000010    0.0006      0.000010    0.27     0.27     |
|   build_sparsity()                 1         0.0011      0.001096    0.0024      0.002392    0.48     1.05     |
|   cnstrn_elem_mat_vec()            32        0.0058      0.000181    0.0058      0.000181    2.56     2.56     |
|   constrain_elem_vector()          32        0.0002      0.000005    0.0002      0.000005    0.08     0.08     |
|   create_dof_constraints()         1         0.0024      0.002364    0.0054      0.005391    1.04     2.37     |
|   distribute_dofs()                1         0.0006      0.000574    0.0018      0.001787    0.25     0.79     |
|   dof_indices()                    443       0.0107      0.000024    0.0107      0.000024    4.72     4.72     |
|   enforce_constraints_exactly()    3         0.0007      0.000242    0.0007      0.000242    0.32     0.32     |
|   prepare_send_list()              1         0.0000      0.000035    0.0000      0.000035    0.02     0.02     |
|   reinit()                         1         0.0008      0.000832    0.0008      0.000832    0.37     0.37     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0004      0.000192    0.0024      0.001222    0.17     1.08     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.1175      0.117547    0.1175      0.117547    51.75    51.75    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        176       0.0338      0.000192    0.0338      0.000192    14.89    14.89    |
|   init_shape_functions()           83        0.0019      0.000023    0.0019      0.000023    0.83     0.83     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       1         0.0080      0.008048    0.0197      0.019714    3.54     8.68     |
|   assembly(get_residual)           1         0.0023      0.002308    0.0078      0.007840    1.02     3.45     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             176       0.0017      0.000010    0.0017      0.000010    0.76     0.76     |
|   compute_face_map()               80        0.0005      0.000006    0.0005      0.000006    0.22     0.22     |
|   init_face_shape_functions()      2         0.0000      0.000017    0.0000      0.000017    0.01     0.01     |
|   init_reference_to_physical_map() 83        0.0014      0.000017    0.0014      0.000017    0.62     0.62     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0008      0.000801    0.0009      0.000905    0.35     0.40     |
|   renumber_nodes_and_elem()        2         0.0001      0.000039    0.0001      0.000039    0.03     0.03     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0011      0.000563    0.0011      0.000563    0.50     0.50     |
|   find_global_indices()            2         0.0003      0.000147    0.0022      0.001112    0.13     0.98     |
|   parallel_sort()                  2         0.0004      0.000188    0.0005      0.000246    0.17     0.22     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0031      0.001531    0.1232      0.061615    1.35     54.25    |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0003      0.000318    0.0003      0.000318    0.14     0.14     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0020      0.002023    0.0030      0.003043    0.89     1.34     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          1         0.0032      0.003225    0.0434      0.043382    1.42     19.10    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0002      0.000019    0.0002      0.000024    0.08     0.09     |
|   max(bool)                        1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                      150       0.0075      0.000050    0.0075      0.000050    3.30     3.30     |
|   max(vector)                      34        0.0004      0.000010    0.0011      0.000032    0.15     0.48     |
|   min(bool)                        175       0.0012      0.000007    0.0012      0.000007    0.52     0.52     |
|   min(scalar)                      143       0.0022      0.000016    0.0022      0.000016    0.99     0.99     |
|   min(vector)                      34        0.0005      0.000013    0.0012      0.000035    0.20     0.52     |
|   probe()                          36        0.0002      0.000007    0.0002      0.000007    0.11     0.11     |
|   receive()                        36        0.0001      0.000004    0.0004      0.000011    0.06     0.17     |
|   send()                           36        0.0001      0.000002    0.0001      0.000002    0.03     0.03     |
|   send_receive()                   40        0.0002      0.000006    0.0008      0.000019    0.11     0.34     |
|   sum()                            27        0.0003      0.000010    0.0004      0.000016    0.11     0.19     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           36        0.0001      0.000001    0.0001      0.000001    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0002      0.000200    0.0004      0.000438    0.09     0.19     |
|   set_parent_processor_ids()       1         0.0000      0.000046    0.0000      0.000046    0.02     0.02     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0114      0.011413    0.0114      0.011413    5.02     5.02     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1959      0.2271                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example vector_fe_ex2:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/vector_fe/vector_fe_ex2'
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
