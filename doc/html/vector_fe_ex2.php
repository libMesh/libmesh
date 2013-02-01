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
        
          virtual Number component( unsigned int component, const Point& p,
        			    const Real )
          {
            const Real x=p(0), y=p(1), z=p(2);
            return soln( component, x, y, z );
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
        
          virtual Gradient component( unsigned int component, const Point& p,
        			    const Real )
          {
            const Real x=p(0), y=p(1), z=p(2);
            return soln( component, x, y, z );
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
        			      const std::string& name,
        			      const unsigned int number)
          : FEMSystem(es, name, number)
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
Create a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
          
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
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Number component( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			    <B><FONT COLOR="#228B22">const</FONT></B> Real )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      <B><FONT COLOR="#A020F0">return</FONT></B> soln( component, x, y, z );
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
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> Gradient component( <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> component, <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			    <B><FONT COLOR="#228B22">const</FONT></B> Real )
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Real x=p(0), y=p(1), z=p(2);
      <B><FONT COLOR="#A020F0">return</FONT></B> soln( component, x, y, z );
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
  			      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
  			      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number)
    : FEMSystem(es, name, number)
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
  
    Mesh mesh;
    
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
***************************************************************
* Running Example vector_fe_ex2:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=216
    n_local_nodes()=43
  n_elem()=125
    n_local_elem()=10
    n_active_elem()=125
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=129
    n_constrained_dofs()=456
    n_local_constrained_dofs()=84
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 47.8056
      Average Off-Processor Bandwidth <= 41.5556
      Maximum  On-Processor Bandwidth <= 135
      Maximum Off-Processor Bandwidth <= 141
    DofMap Constraints
      Number of DoF Constraints = 456
      Number of Heterogenous Constraints= 456
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

Assembling the System
Nonlinear Residual: 4.00845
Linear solve starting, tolerance 1e-12
Linear solve finished, step 16, residual 4.38343e-13
Trying full Newton step
  Current Residual: 4.5843e-13
  Nonlinear solver converged, step 0, residual reduction 1.14366e-13 < 1e-12
  Nonlinear solver relative step size 0.723571 > 1e-06
L2-Error is: 0.142039
H1-Error is: 0.902143
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/vector_fe/vector_fe_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:22:49 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           6.424e-01      1.01752   6.394e-01
Objects:              6.700e+01      1.03077   6.683e+01
Flops:                1.008e+06      0.00000   2.478e+05  2.974e+06
Flops/sec:            1.569e+06      0.00000   3.870e+05  4.644e+06
MPI Messages:         3.860e+02      3.38596   3.092e+02  3.710e+03
MPI Message Lengths:  2.652e+05      3.56014   3.937e+02  1.461e+06
MPI Reductions:       1.390e+02      1.01460

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 6.3940e-01 100.0%  2.9737e+06 100.0%  3.710e+03 100.0%  3.937e+02      100.0%  1.378e+02  99.2% 

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
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

KSPGMRESOrthog        16 1.0 8.1372e-04 4.5 7.49e+04 0.0 0.0e+00 0.0e+00 1.6e+01  0 12  0  0 12   0 12  0  0 12   431
KSPSetUp               2 1.0 7.6056e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 4.6711e-03 1.0 1.01e+06 0.0 1.4e+03 1.8e+02 4.3e+01  1100 38 18 31   1100 38 18 31   635
PCSetUp                2 1.0 2.8551e-03 3.3 2.30e+05 0.0 0.0e+00 0.0e+00 8.8e+00  0 10  0  0  6   0 10  0  0  6   104
PCSetUpOnBlocks        1 1.0 2.3239e-03 6.7 2.30e+05 0.0 0.0e+00 0.0e+00 6.8e+00  0 10  0  0  5   0 10  0  0  5   127
PCApply               18 1.0 4.8661e-04 2.1 3.64e+05 0.0 0.0e+00 0.0e+00 0.0e+00  0 34  0  0  0   0 34  0  0  0  2096
VecMDot               16 1.0 7.4959e-04 6.5 3.74e+04 0.0 0.0e+00 0.0e+00 1.6e+01  0  6  0  0 12   0  6  0  0 12   233
VecNorm               22 1.0 2.3499e-03 8.1 6.07e+03 0.0 0.0e+00 0.0e+00 2.2e+01  0  1  0  0 16   0  1  0  0 16    12
VecScale              17 1.0 4.4107e-05 1.5 2.35e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   250
VecCopy                5 1.0 2.0981e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                30 1.0 1.8358e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                3 1.0 3.6240e-05 2.6 8.28e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   107
VecMAXPY              17 1.0 3.0279e-05 7.5 4.20e+04 0.0 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  6506
VecAssemblyBegin      18 1.0 4.6328e-0256.3 0.00e+00 0.0 2.1e+02 3.1e+02 4.5e+01  5  0  6  5 32   5  0  6  5 33     0
VecAssemblyEnd        18 1.0 6.5327e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       24 1.0 2.3103e-04 3.0 0.00e+00 0.0 2.2e+03 2.5e+02 0.0e+00  0  0 60 38  0   0  0 60 38  0     0
VecScatterEnd         24 1.0 2.0657e-0315.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          17 1.0 2.1441e-0312.0 7.04e+03 0.0 0.0e+00 0.0e+00 1.7e+01  0  1  0  0 12   0  1  0  0 12    15
MatMult               17 1.0 2.1527e-0375.9 3.25e+05 0.0 1.4e+03 1.8e+02 0.0e+00  0 42 38 18  0   0 42 38 18  0   577
MatSolve              18 0.0 2.3913e-04 0.0 3.64e+05 0.0 0.0e+00 0.0e+00 0.0e+00  0 34  0  0  0   0 34  0  0  0  4265
MatLUFactorNum         1 1.0 2.8515e-0413.4 2.30e+05 0.0 0.0e+00 0.0e+00 0.0e+00  0 10  0  0  0   0 10  0  0  0  1039
MatILUFactorSym        1 1.0 1.7409e-0324.2 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatAssemblyBegin       2 1.0 1.4100e-03 3.6 0.00e+00 0.0 1.6e+02 4.5e+03 4.0e+00  0  0  4 49  3   0  0  4 49  3     0
MatAssemblyEnd         2 1.0 1.9109e-03 1.9 0.00e+00 0.0 1.7e+02 4.7e+01 8.0e+00  0  0  5  1  6   0  0  5  1  6     0
MatGetRowIJ            1 0.0 3.8147e-06 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 0.0 1.5593e-04 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.7e+00  0  0  0  0  3   0  0  0  0  3     0
MatZeroEntries         3 0.0 8.7976e-05 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     3              3        20592     0
      Preconditioner     3              3         2608     0
              Vector    38             38       105792     0
      Vector Scatter     5              5         5180     0
           Index Set    12             12        11556     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4       263704     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.94139e-05
Average time for zero size MPI_Send(): 3.39945e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-pc_type bjacobi
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:22:49 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.18615, Active time=0.57907                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0025      0.002533    0.0219      0.021904    0.44     3.78     |
|   build_constraint_matrix()        20        0.0014      0.000071    0.0014      0.000071    0.25     0.25     |
|   build_sparsity()                 1         0.0037      0.003740    0.0116      0.011596    0.65     2.00     |
|   cnstrn_elem_mat_vec()            10        0.0093      0.000927    0.0093      0.000927    1.60     1.60     |
|   constrain_elem_vector()          10        0.0002      0.000022    0.0002      0.000022    0.04     0.04     |
|   create_dof_constraints()         1         0.0118      0.011836    0.0494      0.049390    2.04     8.53     |
|   distribute_dofs()                1         0.0041      0.004082    0.0130      0.013048    0.70     2.25     |
|   dof_indices()                    270       0.0809      0.000300    0.0809      0.000300    13.98    13.98    |
|   enforce_constraints_exactly()    3         0.0012      0.000399    0.0012      0.000399    0.21     0.21     |
|   prepare_send_list()              1         0.0003      0.000262    0.0003      0.000262    0.05     0.05     |
|   reinit()                         1         0.0071      0.007059    0.0071      0.007059    1.22     1.22     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0007      0.000375    0.0082      0.004124    0.13     1.42     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0042      0.004182    0.0042      0.004182    0.72     0.72     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        58        0.1789      0.003084    0.1789      0.003084    30.89    30.89    |
|   init_shape_functions()           31        0.0059      0.000190    0.0059      0.000190    1.02     1.02     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       1         0.0103      0.010284    0.0336      0.033574    1.78     5.80     |
|   assembly(get_residual)           1         0.0021      0.002141    0.0149      0.014921    0.37     2.58     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             58        0.0019      0.000033    0.0019      0.000033    0.33     0.33     |
|   compute_face_map()               28        0.0007      0.000024    0.0007      0.000024    0.12     0.12     |
|   init_face_shape_functions()      2         0.0001      0.000060    0.0001      0.000060    0.02     0.02     |
|   init_reference_to_physical_map() 31        0.0026      0.000085    0.0026      0.000085    0.45     0.45     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0060      0.006005    0.0216      0.021573    1.04     3.73     |
|   renumber_nodes_and_elem()        2         0.0002      0.000117    0.0002      0.000117    0.04     0.04     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0019      0.000947    0.0019      0.000947    0.33     0.33     |
|   find_global_indices()            2         0.0013      0.000635    0.0116      0.005787    0.22     2.00     |
|   parallel_sort()                  2         0.0026      0.001306    0.0068      0.003405    0.45     1.18     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0034      0.001695    0.0161      0.008064    0.59     2.79     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0010      0.000967    0.0010      0.000967    0.17     0.17     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0352      0.035157    0.0383      0.038331    6.07     6.62     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          1         0.0374      0.037432    0.0955      0.095525    6.46     16.50    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0008      0.000088    0.0010      0.000108    0.14     0.17     |
|   max(bool)                        1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      150       0.0776      0.000518    0.0776      0.000518    13.41    13.41    |
|   max(vector)                      34        0.0011      0.000031    0.0036      0.000107    0.18     0.63     |
|   min(bool)                        175       0.0043      0.000025    0.0043      0.000025    0.75     0.75     |
|   min(scalar)                      143       0.0575      0.000402    0.0575      0.000402    9.93     9.93     |
|   min(vector)                      34        0.0012      0.000034    0.0050      0.000147    0.20     0.86     |
|   probe()                          132       0.0020      0.000015    0.0020      0.000015    0.35     0.35     |
|   receive()                        132       0.0008      0.000006    0.0029      0.000022    0.14     0.50     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.07     0.07     |
|   send_receive()                   136       0.0014      0.000010    0.0050      0.000037    0.24     0.87     |
|   sum()                            27        0.0040      0.000147    0.0055      0.000203    0.68     0.94     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0008      0.000787    0.0016      0.001592    0.14     0.27     |
|   set_parent_processor_ids()       1         0.0003      0.000304    0.0003      0.000304    0.05     0.05     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0076      0.007648    0.0076      0.007648    1.32     1.32     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1786      0.5791                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example vector_fe_ex2:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
