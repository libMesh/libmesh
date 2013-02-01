<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("fem_system_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file nonlinear_neohooke_cc.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef NONLINEAR_NEOHOOKE_CC_H_
        #define NONLINEAR_NEOHOOKE_CC_H_
        
        #include "libmesh/dense_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/vector_value.h"
        #include "libmesh/tensor_value.h"
        #include "libmesh/getpot.h"
        
        using namespace libMesh;
        
        /**
         * This class implements a constitutive formulation for an Neo-Hookean elastic solid
         * in terms of the current configuration. This implementation is suitable for computing
         * large deformations. See e.g. "Nonlinear finite element methods" (P. Wriggers, 2008,
         * Springer) for details.
         */
        class NonlinearNeoHookeCurrentConfig {
        public:
          NonlinearNeoHookeCurrentConfig(
              const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi, GetPot& args) :
              dphi(dphi) {
            E = args("material/neohooke/e_modulus", 10000.0);
            nu = args("material/neohooke/nu", 0.3);
          }
        
          /**
           * Initialize the class for the given displacement gradient at the
           * specified quadrature point.
           */
          void init_for_qp(VectorValue&lt;Gradient&gt; & grad_u, unsigned int qp);
        
          /**
           * Return the residual vector for the current state.
           */
          void get_residual(DenseVector&lt;Real&gt; & residuum, unsigned int & i);
        
          /**
           * Return the stiffness matrix for the current state.
           */
          void get_linearized_stiffness(DenseMatrix&lt;Real&gt; & stiffness,
              unsigned int & i, unsigned int & j);
        
          /**
           * Flag to indicate if it is necessary to calculate values for stiffness
           * matrix during initialization.
           */
          bool calculate_linearized_stiffness;
        private:
          void build_b_0_mat(int i, DenseMatrix&lt;Real&gt;& b_l_mat);
          void calculate_stress();
          void calculate_tangent();
          static void tensor_to_voigt(const RealTensor &tensor, DenseVector&lt;Real&gt; &vec);
        
          unsigned int current_qp;
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi;
        
          DenseMatrix&lt;Real&gt; C_mat;
          Real E;
          Real nu;
          RealTensor F, S, tau, sigma;
          DenseMatrix&lt;Real&gt; B_L;
          DenseMatrix&lt;Real&gt; B_K;
        };
        
        #endif /* NONLINEAR_NEOHOOKE_CC_H_ */
        
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file solid_system.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef SOLID_SYSTEM_H_
        #define SOLID_SYSTEM_H_
        
        #include "libmesh/fem_system.h"
        #include "libmesh/auto_ptr.h"
        
        using namespace libMesh;
        
        class SolidSystem: public FEMSystem {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          SolidSystem(EquationSystems& es, const std::string& name,
              const unsigned int number);
        
</pre>
</div>
<div class = "comment">
System initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_data();
        
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
</div>

<div class ="fragment">
<pre>
          virtual bool element_time_derivative(bool request_jacobian,
              DiffContext& context);
        
</pre>
</div>
<div class = "comment">
Contributions for adding boundary conditions
</div>

<div class ="fragment">
<pre>
          virtual bool side_time_derivative(bool request_jacobian,
              DiffContext& context);
        
          virtual bool eulerian_residual(bool, DiffContext &) {
            return false;
          }
        
</pre>
</div>
<div class = "comment">
Simulation parameters
</div>

<div class ="fragment">
<pre>
          GetPot args;
        
</pre>
</div>
<div class = "comment">
Custom Identifier
</div>

<div class ="fragment">
<pre>
          virtual std::string system_type() const {
            return "Solid";
          }
        
</pre>
</div>
<div class = "comment">
override method to update mesh also
</div>

<div class ="fragment">
<pre>
          virtual void update();
        
</pre>
</div>
<div class = "comment">
save the undeformed mesh to an auxiliary system
</div>

<div class ="fragment">
<pre>
          void save_initial_mesh();
        
</pre>
</div>
<div class = "comment">
variable numbers of primary variables in the current system
</div>

<div class ="fragment">
<pre>
          unsigned int var[3];
        
</pre>
</div>
<div class = "comment">
variable numbers of primary variables in auxiliary system (for accessing the
undeformed configuration)
</div>

<div class ="fragment">
<pre>
          unsigned int undefo_var[3];
        };
        
        #endif /* SOLID_SYSTEM_H_ */
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file fem_system_ex2.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/equation_systems.h"
        #include "libmesh/getpot.h"
        #include "libmesh/libmesh.h"
        #include "libmesh/libmesh_logging.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/time_solver.h"
        #include "libmesh/transient_system.h"
        #include "libmesh/vtk_io.h"
        
        #include &lt;cstdio&gt;
        #include &lt;ctime&gt;
        #include &lt;fstream&gt;
        #include &lt;iomanip&gt;
        #include &lt;iostream&gt;
        #include &lt;sstream&gt;
        #include &lt;string&gt;
        #include &lt;unistd.h&gt;
        
        using namespace libMesh;
        
        #include "solid_system.h"
        
        void setup(EquationSystems& systems, Mesh& mesh, GetPot& args)
        {
          const unsigned int dim = mesh.mesh_dimension();
</pre>
</div>
<div class = "comment">
We currently invert tensors with the assumption that they're 3x3
</div>

<div class ="fragment">
<pre>
          libmesh_assert (dim == 3); 
        
</pre>
</div>
<div class = "comment">
Generating Mesh
</div>

<div class ="fragment">
<pre>
          ElemType eltype = Utility::string_to_enum&lt;ElemType&gt;(args("mesh/generation/element_type", "hex8"));
          int nx = args("mesh/generation/num_elem", 4, 0);
          int ny = args("mesh/generation/num_elem", 4, 1);
          int nz = dim &gt; 2 ? args("mesh/generation/num_elem", 4, 2) : 0;
          double origx = args("mesh/generation/origin", -1.0, 0);
          double origy = args("mesh/generation/origin", -1.0, 1);
          double origz = args("mesh/generation/origin", 0.0, 2);
          double sizex = args("mesh/generation/size", 2.0, 0);
          double sizey = args("mesh/generation/size", 2.0, 1);
          double sizez = args("mesh/generation/size", 2.0, 2);
          MeshTools::Generation::build_cube(mesh, nx, ny, nz,
              origx, origx+sizex, origy, origy+sizey, origz, origz+sizez, eltype);
        
</pre>
</div>
<div class = "comment">
Creating Systems
</div>

<div class ="fragment">
<pre>
          SolidSystem& imms = systems.add_system&lt;SolidSystem&gt; ("solid");
          imms.args = args;
        
</pre>
</div>
<div class = "comment">
Build up auxiliary system
</div>

<div class ="fragment">
<pre>
          ExplicitSystem& aux_sys = systems.add_system&lt;TransientExplicitSystem&gt;("auxiliary");
        
</pre>
</div>
<div class = "comment">
Initialize the system
</div>

<div class ="fragment">
<pre>
          systems.parameters.set&lt;unsigned int&gt;("phase") = 0;
          systems.init();
        
          imms.save_initial_mesh();
        
</pre>
</div>
<div class = "comment">
Fill global solution vector from local ones
</div>

<div class ="fragment">
<pre>
          aux_sys.reinit();
        }
        
        
        
        void run_timestepping(EquationSystems& systems, GetPot& args)
        {
          TransientExplicitSystem& aux_system = systems.get_system&lt;TransientExplicitSystem&gt;("auxiliary");
        
          SolidSystem& solid_system = systems.get_system&lt;SolidSystem&gt;("solid");
        
          AutoPtr&lt;VTKIO&gt; io = AutoPtr&lt;VTKIO&gt;(new VTKIO(systems.get_mesh()));
        
          Real duration = args("duration", 1.0);
        
          for (unsigned int t_step = 0; t_step &lt; duration/solid_system.deltat; t_step++) {
</pre>
</div>
<div class = "comment">
Progress in current phase [0..1]
</div>

<div class ="fragment">
<pre>
            Real progress = t_step * solid_system.deltat / duration;
            systems.parameters.set&lt;Real&gt;("progress") = progress;
            systems.parameters.set&lt;unsigned int&gt;("step") = t_step;
        
</pre>
</div>
<div class = "comment">
Update message


<br><br></div>

<div class ="fragment">
<pre>
            out &lt;&lt; "===== Time Step " &lt;&lt; std::setw(4) &lt;&lt; t_step;
            out &lt;&lt; " (" &lt;&lt; std::fixed &lt;&lt; std::setprecision(2) &lt;&lt; std::setw(6) &lt;&lt; (progress*100.) &lt;&lt; "%)";
            out &lt;&lt; ", time = " &lt;&lt; std::setw(7) &lt;&lt; solid_system.time;
            out &lt;&lt; " =====" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Advance timestep in auxiliary system
</div>

<div class ="fragment">
<pre>
            aux_system.current_local_solution-&gt;close();
            aux_system.old_local_solution-&gt;close();
            *aux_system.older_local_solution = *aux_system.old_local_solution;
            *aux_system.old_local_solution = *aux_system.current_local_solution;
        
            out &lt;&lt; "Solving Solid" &lt;&lt; std::endl;
            solid_system.solve();
            aux_system.reinit();
        
</pre>
</div>
<div class = "comment">
Carry out the adaptive mesh refinement/coarsening
</div>

<div class ="fragment">
<pre>
            out &lt;&lt; "Doing a reinit of the equation systems" &lt;&lt; std::endl;
            systems.reinit();
        
            if (t_step % args("output/frequency", 1) == 0) {
              std::string result;
              std::stringstream file_name;
              file_name &lt;&lt; args("results_directory", "./") &lt;&lt; "fem_";
              file_name &lt;&lt; std::setw(6) &lt;&lt; std::setfill('0') &lt;&lt; t_step;
              file_name &lt;&lt; ".pvtu";
        
        
              io-&gt;write_equation_systems(file_name.str(), systems);
            }
</pre>
</div>
<div class = "comment">
Advance to the next timestep in a transient problem
</div>

<div class ="fragment">
<pre>
            out &lt;&lt; "Advancing to next step" &lt;&lt; std::endl;
            solid_system.time_solver-&gt;advance_timestep();
          }
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
        #ifndef LIBMESH_HAVE_VTK
          libmesh_example_assert(false, "--enable-vtk");
        #endif
        
</pre>
</div>
<div class = "comment">
Initialize libMesh and any dependent libraries
</div>

<div class ="fragment">
<pre>
          LibMeshInit init(argc, argv);
        
</pre>
</div>
<div class = "comment">
Threaded assembly doesn't currently work with the moving mesh
code.
We'll skip this example for now.
</div>

<div class ="fragment">
<pre>
          if (libMesh::n_threads() &gt; 1)
            {
              std::cout &lt;&lt; "We skip fem_system_ex2 when using threads.\n"
                        &lt;&lt; std::endl;
              return 0;
            }
        
</pre>
</div>
<div class = "comment">
read simulation parameters from file
</div>

<div class ="fragment">
<pre>
          GetPot args = GetPot("solid.in");
        
</pre>
</div>
<div class = "comment">
Create System and Mesh
</div>

<div class ="fragment">
<pre>
          int dim = args("mesh/generation/dimension", 3);
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "3D support");
        
          Mesh mesh(dim);
          EquationSystems systems(mesh);
        
</pre>
</div>
<div class = "comment">
Create and set systems up
</div>

<div class ="fragment">
<pre>
          setup(systems, mesh, args);
        
</pre>
</div>
<div class = "comment">
run the systems
</div>

<div class ="fragment">
<pre>
          run_timestepping(systems, args);
        
          out &lt;&lt; "Finished calculations" &lt;&lt; std::endl;
          return 0;
        }
        
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file nonlinear_neohooke_cc.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "nonlinear_neohooke_cc.h"
        
        /**
         * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
         * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
         */
        template &lt;typename T&gt; TypeTensor&lt;T&gt; inv(const TypeTensor&lt;T&gt; &A ) {
          double Sub11, Sub12, Sub13;
          Sub11 = A._coords[4]*A._coords[8] - A._coords[5]*A._coords[7];
          Sub12 = A._coords[3]*A._coords[8] - A._coords[6]*A._coords[5];
          Sub13 = A._coords[3]*A._coords[7] - A._coords[6]*A._coords[4];
          double detA = A._coords[0]*Sub11 - A._coords[1]*Sub12 + A._coords[2]*Sub13;
          libmesh_assert( std::fabs(detA)&gt;1.e-15 );
        
          TypeTensor&lt;T&gt; Ainv(A);
        
          Ainv._coords[0] =  Sub11/detA;
          Ainv._coords[1] = (-A._coords[1]*A._coords[8]+A._coords[2]*A._coords[7])/detA;
          Ainv._coords[2] = ( A._coords[1]*A._coords[5]-A._coords[2]*A._coords[4])/detA;
          Ainv._coords[3] = -Sub12/detA;
          Ainv._coords[4] = ( A._coords[0]*A._coords[8]-A._coords[2]*A._coords[6])/detA;
          Ainv._coords[5] = (-A._coords[0]*A._coords[5]+A._coords[2]*A._coords[3])/detA;
          Ainv._coords[6] =  Sub13/detA;
          Ainv._coords[7] = (-A._coords[0]*A._coords[7]+A._coords[1]*A._coords[6])/detA;
          Ainv._coords[8] = ( A._coords[0]*A._coords[4]-A._coords[1]*A._coords[3])/detA;
        
          return Ainv;
        }
        
        void NonlinearNeoHookeCurrentConfig::init_for_qp(VectorValue&lt;Gradient&gt; & grad_u, unsigned int qp) {
        	this-&gt;current_qp = qp;
        	F.zero();
        	S.zero();
        
        	{
        	  RealTensor invF;
        	  invF.zero();
        	  for (unsigned int i = 0; i &lt; 3; ++i)
        	    for (unsigned int j = 0; j &lt; 3; ++j) {
        	      invF(i, j) += grad_u(i)(j);
        	    }
            F.add(inv(invF));
        	}
        
        	if (F.det() &lt; -TOLERANCE) {
        		std::cout &lt;&lt; "detF &lt; 0" &lt;&lt; std::endl;
        		libmesh_error();
        	}
        
        	if (this-&gt;calculate_linearized_stiffness) {
        		this-&gt;calculate_tangent();
        	}
        
        	this-&gt;calculate_stress();
        }
        
        
        
        void NonlinearNeoHookeCurrentConfig::calculate_tangent() {
        	Real mu = E / (2 * (1 + nu));
        	Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        
        	Real detF = F.det();
        
        	C_mat.resize(6, 6);
        	for (unsigned int i = 0; i &lt; 3; ++i) {
        		for (unsigned int j = 0; j &lt; 3; ++j) {
        			if (i == j) {
        				C_mat(i, j) = 2 * mu + lambda;
        				C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
        			} else {
        				C_mat(i, j) = lambda * detF * detF;
        			}
        		}
        	}
        }
        
        
        void NonlinearNeoHookeCurrentConfig::calculate_stress() {
        
          double mu = E / (2.0 * (1.0 + nu));
        	double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
        
        	Real detF = F.det();
        	RealTensor Ft = F.transpose();
        
        	RealTensor C = Ft * F;
        	RealTensor b = F * Ft;
        	RealTensor identity;
        	identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
        	RealTensor invC = inv(C);
        
        	S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);
        
        	tau = (F * S) * Ft;
        	sigma = 1/detF * tau;
        }
        
        void NonlinearNeoHookeCurrentConfig::get_residual(DenseVector&lt;Real&gt; & residuum, unsigned int & i) {
        	B_L.resize(3, 6);
        	DenseVector&lt;Real&gt; sigma_voigt(6);
        
        	this-&gt;build_b_0_mat(i, B_L);
        
        	tensor_to_voigt(sigma, sigma_voigt);
        
        	B_L.vector_mult(residuum, sigma_voigt);
        }
        
        void NonlinearNeoHookeCurrentConfig::tensor_to_voigt(const RealTensor &tensor, DenseVector&lt;Real&gt; &vec) {
          vec(0) = tensor(0, 0);
          vec(1) = tensor(1, 1);
          vec(2) = tensor(2, 2);
          vec(3) = tensor(0, 1);
          vec(4) = tensor(1, 2);
          vec(5) = tensor(0, 2);
        
        }
        
        void NonlinearNeoHookeCurrentConfig::get_linearized_stiffness(DenseMatrix&lt;Real&gt; & stiffness, unsigned int & i, unsigned int & j) {
        	stiffness.resize(3, 3);
        
        	double G_IK = (sigma * dphi[i][current_qp]) * dphi[j][current_qp];
        	stiffness(0, 0) += G_IK;
        	stiffness(1, 1) += G_IK;
        	stiffness(2, 2) += G_IK;
        
        	B_L.resize(3, 6);
        	this-&gt;build_b_0_mat(i, B_L);
        	B_K.resize(3, 6);
        	this-&gt;build_b_0_mat(j, B_K);
        
        	B_L.right_multiply(C_mat);
        	B_L.right_multiply_transpose(B_K);
        	B_L *= 1/F.det();
        
        	stiffness += B_L;
        }
        
        void NonlinearNeoHookeCurrentConfig::build_b_0_mat(int i, DenseMatrix&lt;Real&gt;& b_0_mat) {
        	for (unsigned int ii = 0; ii &lt; 3; ++ii) {
        		b_0_mat(ii, ii) = dphi[i][current_qp](ii);
        	}
        	b_0_mat(0, 3) = dphi[i][current_qp](1);
        	b_0_mat(1, 3) = dphi[i][current_qp](0);
        	b_0_mat(1, 4) = dphi[i][current_qp](2);
        	b_0_mat(2, 4) = dphi[i][current_qp](1);
        	b_0_mat(0, 5) = dphi[i][current_qp](2);
        	b_0_mat(2, 5) = dphi[i][current_qp](0);
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file solid_system.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/boundary_info.h"
        #include "libmesh/diff_solver.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe_base.h"
        #include "libmesh/fem_context.h"
        #include "libmesh/getpot.h"
        #include "libmesh/mesh.h"
        #include "libmesh/newton_solver.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/steady_solver.h"
        #include "libmesh/transient_system.h"
        
        #include "nonlinear_neohooke_cc.h"
        #include "solid_system.h"
        
</pre>
</div>
<div class = "comment">
Solaris Studio has no NAN
</div>

<div class ="fragment">
<pre>
        #ifdef __SUNPRO_CC
          #define NAN (1.0/0.0)
        #endif
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        SolidSystem::SolidSystem(EquationSystems& es, const std::string& name,
            const unsigned int number) :
            FEMSystem(es, name, number) {
        
</pre>
</div>
<div class = "comment">
Add a time solver. We are just looking at a steady state problem here.
</div>

<div class ="fragment">
<pre>
          this-&gt;time_solver = AutoPtr&lt;TimeSolver&gt;(new SteadySolver(*this));
        }
        
        void SolidSystem::save_initial_mesh() {
          System & aux_sys = this-&gt;get_equation_systems().get_system("auxiliary");
        
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Loop over all nodes and copy the location from the current system to
the auxiliary system.
</div>

<div class ="fragment">
<pre>
          const MeshBase::const_node_iterator nd_end =
              this-&gt;get_mesh().local_nodes_end();
          for (MeshBase::const_node_iterator nd = this-&gt;get_mesh().local_nodes_begin();
              nd != nd_end; ++nd) {
            const Node *node = *nd;
            for (unsigned int d = 0; d &lt; dim; ++d) {
              unsigned int source_dof = node-&gt;dof_number(this-&gt;number(), var[d], 0);
              unsigned int dest_dof = node-&gt;dof_number(aux_sys.number(), undefo_var[d],
                  0);
              Real value = this-&gt;current_local_solution-&gt;el(source_dof);
              aux_sys.current_local_solution-&gt;set(dest_dof, value);
            }
          }
        }
        
        void SolidSystem::init_data() {
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get the default order of the used elements. Assumption:
Just one type of elements in the mesh.
</div>

<div class ="fragment">
<pre>
          Order order = (*(this-&gt;get_mesh().elements_begin()))-&gt;default_order();
        
</pre>
</div>
<div class = "comment">
Add the node positions as primary variables.
</div>

<div class ="fragment">
<pre>
          var[0] = this-&gt;add_variable("x", order);
          var[1] = this-&gt;add_variable("y", order);
          if (dim == 3)
            var[2] = this-&gt;add_variable("z", order);
          else
            var[2] = var[1];
        
</pre>
</div>
<div class = "comment">
Add variables for storing the initial mesh to an auxiliary system.
</div>

<div class ="fragment">
<pre>
          System& aux_sys = this-&gt;get_equation_systems().get_system("auxiliary");
          undefo_var[0] = aux_sys.add_variable("undefo_x", order);
          undefo_var[1] = aux_sys.add_variable("undefo_y", order);
          undefo_var[2] = aux_sys.add_variable("undefo_z", order);
        
</pre>
</div>
<div class = "comment">
Set the time stepping options
</div>

<div class ="fragment">
<pre>
          this-&gt;deltat = args("schedule/dt", 0.2);
        
</pre>
</div>
<div class = "comment">
Do the parent's initialization after variables are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        
</pre>
</div>
<div class = "comment">
// Tell the system to march velocity forward in time, but
// leave p as a constraint only
this->time_evolving(var[0]);
this->time_evolving(var[1]);
if (dim == 3)
this->time_evolving(var[2]);


<br><br>Tell the system which variables are containing the node positions
</div>

<div class ="fragment">
<pre>
          set_mesh_system(this);
        
          this-&gt;set_mesh_x_var(var[0]);
          this-&gt;set_mesh_y_var(var[1]);
          if (dim == 3)
            this-&gt;set_mesh_z_var(var[2]);
        
</pre>
</div>
<div class = "comment">
Fill the variables with the position of the nodes
</div>

<div class ="fragment">
<pre>
          this-&gt;mesh_position_get();
        
          System::reinit();
        
</pre>
</div>
<div class = "comment">
Set some options for the DiffSolver
</div>

<div class ="fragment">
<pre>
          DiffSolver &solver = *(this-&gt;time_solver-&gt;diff_solver().get());
          solver.quiet = args("solver/quiet", false);
          solver.max_nonlinear_iterations = args(
              "solver/nonlinear/max_nonlinear_iterations", 100);
          solver.relative_step_tolerance = args(
              "solver/nonlinear/relative_step_tolerance", 1.e-3);
          solver.relative_residual_tolerance = args(
              "solver/nonlinear/relative_residual_tolerance", 1.e-8);
          solver.absolute_residual_tolerance = args(
              "solver/nonlinear/absolute_residual_tolerance", 1.e-8);
          solver.verbose = !args("solver/quiet", false);
        
          ((NewtonSolver&) solver).require_residual_reduction = args(
              "solver/nonlinear/require_reduction", false);
        
</pre>
</div>
<div class = "comment">
And the linear solver options
</div>

<div class ="fragment">
<pre>
          solver.max_linear_iterations = args("max_linear_iterations", 50000);
          solver.initial_linear_tolerance = args("initial_linear_tolerance", 1.e-3);
        }
        
        void SolidSystem::update() {
          System::update();
          this-&gt;mesh_position_set();
        }
        
        void SolidSystem::init_context(DiffContext &context) {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Pre-request all the data needed
</div>

<div class ="fragment">
<pre>
          c.element_fe_var[var[0]]-&gt;get_JxW();
          c.element_fe_var[var[0]]-&gt;get_phi();
          c.element_fe_var[var[0]]-&gt;get_dphi();
          c.element_fe_var[var[0]]-&gt;get_xyz();
        
          c.side_fe_var[var[0]]-&gt;get_JxW();
          c.side_fe_var[var[0]]-&gt;get_phi();
          c.side_fe_var[var[0]]-&gt;get_xyz();
        }
        
        /**
         * Compute contribution from internal forces in elem_residual and contribution from
         * linearized internal forces (stiffness matrix) in elem_jacobian.
         */
        bool SolidSystem::element_time_derivative(bool request_jacobian,
            DiffContext &context) {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = c.element_fe_var[var[0]]-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The gradients of the shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi =
              c.element_fe_var[var[0]]-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
Dimension of the mesh
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = this-&gt;get_mesh().mesh_dimension();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = c.dof_indices_var[var[0]].size();
          libmesh_assert(n_u_dofs == c.dof_indices_var[var[1]].size());
          if (dim == 3) {
            libmesh_assert(n_u_dofs == c.dof_indices_var[var[2]].size());
          }
        
          unsigned int n_qpoints = c.element_qrule-&gt;n_points();
        
</pre>
</div>
<div class = "comment">
Some matrices and vectors for storing the results of the constitutive
law
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Real&gt; stiff;
          DenseVector&lt;Real&gt; res;
          VectorValue&lt;Gradient&gt; grad_u;
        
</pre>
</div>
<div class = "comment">
Instantiate the constitutive law
</div>

<div class ="fragment">
<pre>
          NonlinearNeoHookeCurrentConfig material(dphi, args);
        
</pre>
</div>
<div class = "comment">
Just calculate jacobian contribution when we need to
</div>

<div class ="fragment">
<pre>
          material.calculate_linearized_stiffness = request_jacobian;
        
</pre>
</div>
<div class = "comment">
Get a reference to the auxiliary system
</div>

<div class ="fragment">
<pre>
          TransientExplicitSystem& aux_system = this-&gt;get_equation_systems().get_system&lt;
              TransientExplicitSystem&gt;("auxiliary");
          std::vector&lt;dof_id_type&gt; undefo_index;
        
</pre>
</div>
<div class = "comment">
Assume symmetry of local stiffness matrices
</div>

<div class ="fragment">
<pre>
          bool use_symmetry = args("assembly/use_symmetry", false);
        
</pre>
</div>
<div class = "comment">
Now we will build the element Jacobian and residual.
This must be calculated at each quadrature point by
summing the solution degree-of-freedom values by
the appropriate weight functions.
This class just takes care of the assembly. The matrix of
the jacobian and the residual vector are provided by the
constitutive formulation.


<br><br></div>

<div class ="fragment">
<pre>
          for (unsigned int qp = 0; qp != n_qpoints; qp++) {
</pre>
</div>
<div class = "comment">
Compute the displacement gradient
</div>

<div class ="fragment">
<pre>
            grad_u(0) = grad_u(1) = grad_u(2) = 0;
            for (unsigned int d = 0; d &lt; dim; ++d) {
              std::vector&lt;Number&gt; u_undefo;
              aux_system.get_dof_map().dof_indices(c.elem, undefo_index, undefo_var[d]);
              aux_system.current_local_solution-&gt;get(undefo_index, u_undefo);
              for (unsigned int l = 0; l != n_u_dofs; l++)
                grad_u(d).add_scaled(dphi[l][qp], u_undefo[l]); // u_current(l)); // -
            }
        
</pre>
</div>
<div class = "comment">
initialize the constitutive formulation with the current displacement
gradient
</div>

<div class ="fragment">
<pre>
            material.init_for_qp(grad_u, qp);
        
</pre>
</div>
<div class = "comment">
Aquire, scale and assemble residual and stiffness
</div>

<div class ="fragment">
<pre>
            for (unsigned int i = 0; i &lt; n_u_dofs; i++) {
              res.resize(dim);
              material.get_residual(res, i);
              res.scale(JxW[qp]);
              for (unsigned int ii = 0; ii &lt; dim; ++ii) {
                c.elem_subresiduals[ii]-&gt;operator ()(i) += res(ii);
              }
        
              if (request_jacobian && c.elem_solution_derivative) {
                libmesh_assert(c.elem_solution_derivative == 1.0);
                for (unsigned int j = (use_symmetry ? i : 0); j &lt; n_u_dofs; j++) {
                  material.get_linearized_stiffness(stiff, i, j);
                  stiff.scale(JxW[qp]);
                  for (unsigned int ii = 0; ii &lt; dim; ++ii) {
                    for (unsigned int jj = 0; jj &lt; dim; ++jj) {
                      c.elem_subjacobians[ii][jj]-&gt;operator ()(i, j) += stiff(ii, jj);
                      if (use_symmetry && i != j) {
                        c.elem_subjacobians[ii][jj]-&gt;operator ()(j, i) += stiff(jj, ii);
                      }
                    }
                  }
                }
              }
            }
          } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
        bool SolidSystem::side_time_derivative(bool request_jacobian,
            DiffContext &context) {
          FEMContext &c = libmesh_cast_ref&lt;FEMContext&&gt;(context);
        
</pre>
</div>
<div class = "comment">
Apply displacement boundary conditions with penalty method


<br><br>Get the current load step
</div>

<div class ="fragment">
<pre>
          Real ratio = this-&gt;get_equation_systems().parameters.get&lt;Real&gt;("progress")
              + 0.001;
        
</pre>
</div>
<div class = "comment">
The BC are stored in the simulation parameters as array containing sequences of
four numbers: Id of the side for the displacements and three values describing the
displacement. E.g.: bc/displacement = '5 nan nan -1.0'. This will move all nodes of
side 5 about 1.0 units down the z-axis while leaving all other directions unrestricted


<br><br>Get number of BCs to enforce
</div>

<div class ="fragment">
<pre>
          unsigned int num_bc = args.vector_variable_size("bc/displacement");
          if (num_bc % 4 != 0) {
            libMesh::err
                &lt;&lt; "ERROR, Odd number of values in displacement boundary condition.\n"
                &lt;&lt; std::endl;
            libmesh_error();
          }
          num_bc /= 4;
        
</pre>
</div>
<div class = "comment">
Loop over all BCs
</div>

<div class ="fragment">
<pre>
          for (unsigned int nbc = 0; nbc &lt; num_bc; nbc++) {
</pre>
</div>
<div class = "comment">
Get IDs of the side for this BC
</div>

<div class ="fragment">
<pre>
            short int positive_boundary_id = args("bc/displacement", 1, nbc * 4);
        
</pre>
</div>
<div class = "comment">
The current side may not be on the boundary to be restricted
</div>

<div class ="fragment">
<pre>
            if (!this-&gt;get_mesh().boundary_info-&gt;has_boundary_id
        	  (c.elem,c.side,positive_boundary_id))
              continue;
        
</pre>
</div>
<div class = "comment">
Read values from configuration file
</div>

<div class ="fragment">
<pre>
            Point diff_value;
            for (unsigned int d = 0; d &lt; c.dim; ++d) {
              diff_value(d) = args("bc/displacement", NAN, nbc * 4 + 1 + d);
            }
</pre>
</div>
<div class = "comment">
Scale according to current load step
</div>

<div class ="fragment">
<pre>
            diff_value *= ratio;
        
            Real penalty_number = args("bc/displacement_penalty", 1e7);
        
            FEBase * fe = c.side_fe_var[var[0]];
            const std::vector&lt;std::vector&lt;Real&gt; &gt; & phi = fe-&gt;get_phi();
            const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
            const std::vector&lt;Point&gt;& coords = fe-&gt;get_xyz();
        
            unsigned int n_x_dofs = c.dof_indices_var[this-&gt;var[0]].size();
        
</pre>
</div>
<div class = "comment">
get mappings for dofs for auxiliary system for original mesh positions
</div>

<div class ="fragment">
<pre>
            const System & auxsys = this-&gt;get_equation_systems().get_system(
                "auxiliary");
            const DofMap & auxmap = auxsys.get_dof_map();
            std::vector&lt;dof_id_type&gt; undefo_dofs[3];
            for (unsigned int d = 0; d &lt; c.dim; ++d) {
              auxmap.dof_indices(c.elem, undefo_dofs[d], undefo_var[d]);
            }
        
            for (unsigned int qp = 0; qp &lt; c.side_qrule-&gt;n_points(); ++qp) {
</pre>
</div>
<div class = "comment">
calculate coordinates of qp on undeformed mesh
</div>

<div class ="fragment">
<pre>
              Point orig_point;
              for (unsigned int i = 0; i &lt; n_x_dofs; ++i) {
                for (unsigned int d = 0; d &lt; c.dim; ++d) {
                  Real orig_val = auxsys.current_solution(undefo_dofs[d][i]);
                  orig_point(d) += phi[i][qp] * orig_val;
                }
              }
        
</pre>
</div>
<div class = "comment">
Calculate displacement to be enforced.
</div>

<div class ="fragment">
<pre>
              Point diff = coords[qp] - orig_point - diff_value;
        
</pre>
</div>
<div class = "comment">
Assemble
</div>

<div class ="fragment">
<pre>
              for (unsigned int i = 0; i &lt; n_x_dofs; ++i) {
                for (unsigned int d1 = 0; d1 &lt; c.dim; ++d1) {
                  if (libmesh_isnan(diff(d1)))
                    continue;
                  Real val = JxW[qp] * phi[i][qp] * diff(d1) * penalty_number;
                  c.elem_subresiduals[var[d1]]-&gt;operator ()(i) += val;
                }
                if (request_jacobian) {
                  for (unsigned int j = 0; j &lt; n_x_dofs; ++j) {
                    for (unsigned int d1 = 0; d1 &lt; c.dim; ++d1) {
                      if (libmesh_isnan(diff(d1)))
                        continue;
                      Real val = JxW[qp] * phi[i][qp] * phi[j][qp] * penalty_number;
                      c.elem_subjacobians[var[d1]][var[d1]]-&gt;operator ()(i, j) += val;
                    }
                  }
                }
              }
            }
          }
        
          return request_jacobian;
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file nonlinear_neohooke_cc.h without comments: </h1> 
<pre> 
  #ifndef NONLINEAR_NEOHOOKE_CC_H_
  #define NONLINEAR_NEOHOOKE_CC_H_
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tensor_value.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <I><FONT COLOR="#B22222">/**
   * This class implements a constitutive formulation for an Neo-Hookean elastic solid
   * in terms of the current configuration. This implementation is suitable for computing
   * large deformations. See e.g. &quot;Nonlinear finite element methods&quot; (P. Wriggers, 2008,
   * Springer) for details.
   */</FONT></I>
  <B><FONT COLOR="#228B22">class</FONT></B> NonlinearNeoHookeCurrentConfig {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    NonlinearNeoHookeCurrentConfig(
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi, GetPot&amp; args) :
        dphi(dphi) {
      E = args(<B><FONT COLOR="#BC8F8F">&quot;material/neohooke/e_modulus&quot;</FONT></B>, 10000.0);
      nu = args(<B><FONT COLOR="#BC8F8F">&quot;material/neohooke/nu&quot;</FONT></B>, 0.3);
    }
  
    <I><FONT COLOR="#B22222">/**
     * Initialize the class for the given displacement gradient at the
     * specified quadrature point.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> init_for_qp(VectorValue&lt;Gradient&gt; &amp; grad_u, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp);
  
    <I><FONT COLOR="#B22222">/**
     * Return the residual vector for the current state.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> get_residual(DenseVector&lt;Real&gt; &amp; residuum, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; i);
  
    <I><FONT COLOR="#B22222">/**
     * Return the stiffness matrix for the current state.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> get_linearized_stiffness(DenseMatrix&lt;Real&gt; &amp; stiffness,
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; i, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; j);
  
    <I><FONT COLOR="#B22222">/**
     * Flag to indicate if it is necessary to calculate values for stiffness
     * matrix during initialization.
     */</FONT></I>
    <B><FONT COLOR="#228B22">bool</FONT></B> calculate_linearized_stiffness;
  <B><FONT COLOR="#228B22">private</FONT></B>:
    <B><FONT COLOR="#228B22">void</FONT></B> build_b_0_mat(<B><FONT COLOR="#228B22">int</FONT></B> i, DenseMatrix&lt;Real&gt;&amp; b_l_mat);
    <B><FONT COLOR="#228B22">void</FONT></B> calculate_stress();
    <B><FONT COLOR="#228B22">void</FONT></B> calculate_tangent();
    <B><FONT COLOR="#228B22">static</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> tensor_to_voigt(<B><FONT COLOR="#228B22">const</FONT></B> RealTensor &amp;tensor, DenseVector&lt;Real&gt; &amp;vec);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> current_qp;
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi;
  
    DenseMatrix&lt;Real&gt; C_mat;
    Real E;
    Real nu;
    RealTensor F, S, tau, sigma;
    DenseMatrix&lt;Real&gt; B_L;
    DenseMatrix&lt;Real&gt; B_K;
  };
  
  #endif <I><FONT COLOR="#B22222">/* NONLINEAR_NEOHOOKE_CC_H_ */</FONT></I>
  
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file solid_system.h without comments: </h1> 
<pre> 
  #ifndef SOLID_SYSTEM_H_
  #define SOLID_SYSTEM_H_
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/auto_ptr.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">class</FONT></B> SolidSystem: <B><FONT COLOR="#228B22">public</FONT></B> FEMSystem {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    SolidSystem(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_data();
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> init_context(DiffContext &amp;context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> element_time_derivative(<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
        DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> side_time_derivative(<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
        DiffContext&amp; context);
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> eulerian_residual(<B><FONT COLOR="#228B22">bool</FONT></B>, DiffContext &amp;) {
      <B><FONT COLOR="#A020F0">return</FONT></B> false;
    }
  
    GetPot args;
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> std::string system_type() <B><FONT COLOR="#228B22">const</FONT></B> {
      <B><FONT COLOR="#A020F0">return</FONT></B> <B><FONT COLOR="#BC8F8F">&quot;Solid&quot;</FONT></B>;
    }
  
    <B><FONT COLOR="#228B22">virtual</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> update();
  
    <B><FONT COLOR="#228B22">void</FONT></B> save_initial_mesh();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> var[3];
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> undefo_var[3];
  };
  
  #endif <I><FONT COLOR="#B22222">/* SOLID_SYSTEM_H_ */</FONT></I>
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file fem_system_ex2.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_logging.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/time_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vtk_io.h&quot;</FONT></B>
  
  #include &lt;cstdio&gt;
  #include &lt;ctime&gt;
  #include &lt;fstream&gt;
  #include &lt;iomanip&gt;
  #include &lt;iostream&gt;
  #include &lt;sstream&gt;
  #include &lt;string&gt;
  #include &lt;unistd.h&gt;
  
  using namespace libMesh;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;solid_system.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">void</FONT></B> setup(EquationSystems&amp; systems, Mesh&amp; mesh, GetPot&amp; args)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    libmesh_assert (dim == 3); 
  
    ElemType eltype = Utility::string_to_enum&lt;ElemType&gt;(args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/element_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;hex8&quot;</FONT></B>));
    <B><FONT COLOR="#228B22">int</FONT></B> nx = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/num_elem&quot;</FONT></B>, 4, 0);
    <B><FONT COLOR="#228B22">int</FONT></B> ny = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/num_elem&quot;</FONT></B>, 4, 1);
    <B><FONT COLOR="#228B22">int</FONT></B> nz = dim &gt; 2 ? args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/num_elem&quot;</FONT></B>, 4, 2) : 0;
    <B><FONT COLOR="#228B22">double</FONT></B> origx = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/origin&quot;</FONT></B>, -1.0, 0);
    <B><FONT COLOR="#228B22">double</FONT></B> origy = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/origin&quot;</FONT></B>, -1.0, 1);
    <B><FONT COLOR="#228B22">double</FONT></B> origz = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/origin&quot;</FONT></B>, 0.0, 2);
    <B><FONT COLOR="#228B22">double</FONT></B> sizex = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/size&quot;</FONT></B>, 2.0, 0);
    <B><FONT COLOR="#228B22">double</FONT></B> sizey = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/size&quot;</FONT></B>, 2.0, 1);
    <B><FONT COLOR="#228B22">double</FONT></B> sizez = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/size&quot;</FONT></B>, 2.0, 2);
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube(mesh, nx, ny, nz,
        origx, origx+sizex, origy, origy+sizey, origz, origz+sizez, eltype);
  
    SolidSystem&amp; imms = systems.add_system&lt;SolidSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;solid&quot;</FONT></B>);
    imms.args = args;
  
    ExplicitSystem&amp; aux_sys = systems.add_system&lt;TransientExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
  
    systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;phase&quot;</FONT></B>) = 0;
    systems.init();
  
    imms.save_initial_mesh();
  
    aux_sys.reinit();
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> run_timestepping(EquationSystems&amp; systems, GetPot&amp; args)
  {
    TransientExplicitSystem&amp; aux_system = systems.get_system&lt;TransientExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
  
    SolidSystem&amp; solid_system = systems.get_system&lt;SolidSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;solid&quot;</FONT></B>);
  
    AutoPtr&lt;VTKIO&gt; io = AutoPtr&lt;VTKIO&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> VTKIO(systems.get_mesh()));
  
    Real duration = args(<B><FONT COLOR="#BC8F8F">&quot;duration&quot;</FONT></B>, 1.0);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step = 0; t_step &lt; duration/solid_system.deltat; t_step++) {
      Real progress = t_step * solid_system.deltat / duration;
      systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;progress&quot;</FONT></B>) = progress;
      systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;step&quot;</FONT></B>) = t_step;
  
  
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;===== Time Step &quot;</FONT></B> &lt;&lt; std::setw(4) &lt;&lt; t_step;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; (&quot;</FONT></B> &lt;&lt; std::fixed &lt;&lt; std::setprecision(2) &lt;&lt; std::setw(6) &lt;&lt; (progress*100.) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;%)&quot;</FONT></B>;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B> &lt;&lt; std::setw(7) &lt;&lt; solid_system.time;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; =====&quot;</FONT></B> &lt;&lt; std::endl;
  
      aux_system.current_local_solution-&gt;close();
      aux_system.old_local_solution-&gt;close();
      *aux_system.older_local_solution = *aux_system.old_local_solution;
      *aux_system.old_local_solution = *aux_system.current_local_solution;
  
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solving Solid&quot;</FONT></B> &lt;&lt; std::endl;
      solid_system.solve();
      aux_system.reinit();
  
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Doing a reinit of the equation systems&quot;</FONT></B> &lt;&lt; std::endl;
      systems.reinit();
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (t_step % args(<B><FONT COLOR="#BC8F8F">&quot;output/frequency&quot;</FONT></B>, 1) == 0) {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::string result;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::stringstream file_name;
        file_name &lt;&lt; args(<B><FONT COLOR="#BC8F8F">&quot;results_directory&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;./&quot;</FONT></B>) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;fem_&quot;</FONT></B>;
        file_name &lt;&lt; std::setw(6) &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>) &lt;&lt; t_step;
        file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.pvtu&quot;</FONT></B>;
  
  
        io-&gt;write_equation_systems(file_name.str(), systems);
      }
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Advancing to next step&quot;</FONT></B> &lt;&lt; std::endl;
      solid_system.time_solver-&gt;advance_timestep();
    }
  }
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
  #ifndef LIBMESH_HAVE_VTK
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-vtk&quot;</FONT></B>);
  #endif
  
    LibMeshInit init(argc, argv);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::n_threads() &gt; 1)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;We skip fem_system_ex2 when using threads.\n&quot;</FONT></B>
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
      }
  
    GetPot args = GetPot(<B><FONT COLOR="#BC8F8F">&quot;solid.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">int</FONT></B> dim = args(<B><FONT COLOR="#BC8F8F">&quot;mesh/generation/dimension&quot;</FONT></B>, 3);
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  
    Mesh mesh(dim);
    EquationSystems systems(mesh);
  
    setup(systems, mesh, args);
  
    run_timestepping(systems, args);
  
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Finished calculations&quot;</FONT></B> &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file nonlinear_neohooke_cc.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;nonlinear_neohooke_cc.h&quot;</FONT></B>
  
  <I><FONT COLOR="#B22222">/**
   * Return the inverse of the given TypeTensor. Algorithm taken from the tensor classes
   * of Ton van den Boogaard (http://tmku209.ctw.utwente.nl/~ton/tensor.html)
   */</FONT></I>
  <B><FONT COLOR="#228B22">template</FONT></B> &lt;typename T&gt; TypeTensor&lt;T&gt; inv(<B><FONT COLOR="#228B22">const</FONT></B> TypeTensor&lt;T&gt; &amp;A ) {
    <B><FONT COLOR="#228B22">double</FONT></B> Sub11, Sub12, Sub13;
    Sub11 = A._coords[4]*A._coords[8] - A._coords[5]*A._coords[7];
    Sub12 = A._coords[3]*A._coords[8] - A._coords[6]*A._coords[5];
    Sub13 = A._coords[3]*A._coords[7] - A._coords[6]*A._coords[4];
    <B><FONT COLOR="#228B22">double</FONT></B> detA = A._coords[0]*Sub11 - A._coords[1]*Sub12 + A._coords[2]*Sub13;
    libmesh_assert( std::fabs(detA)&gt;1.e-15 );
  
    TypeTensor&lt;T&gt; Ainv(A);
  
    Ainv._coords[0] =  Sub11/detA;
    Ainv._coords[1] = (-A._coords[1]*A._coords[8]+A._coords[2]*A._coords[7])/detA;
    Ainv._coords[2] = ( A._coords[1]*A._coords[5]-A._coords[2]*A._coords[4])/detA;
    Ainv._coords[3] = -Sub12/detA;
    Ainv._coords[4] = ( A._coords[0]*A._coords[8]-A._coords[2]*A._coords[6])/detA;
    Ainv._coords[5] = (-A._coords[0]*A._coords[5]+A._coords[2]*A._coords[3])/detA;
    Ainv._coords[6] =  Sub13/detA;
    Ainv._coords[7] = (-A._coords[0]*A._coords[7]+A._coords[1]*A._coords[6])/detA;
    Ainv._coords[8] = ( A._coords[0]*A._coords[4]-A._coords[1]*A._coords[3])/detA;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> Ainv;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::init_for_qp(VectorValue&lt;Gradient&gt; &amp; grad_u, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp) {
  	<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;current_qp = qp;
  	F.zero();
  	S.zero();
  
  	{
  	  RealTensor invF;
  	  invF.zero();
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; 3; ++i)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j &lt; 3; ++j) {
  	      invF(i, j) += grad_u(i)(j);
  	    }
      F.add(inv(invF));
  	}
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (F.det() &lt; -TOLERANCE) {
  		<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;detF &lt; 0&quot;</FONT></B> &lt;&lt; std::endl;
  		libmesh_error();
  	}
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;calculate_linearized_stiffness) {
  		<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;calculate_tangent();
  	}
  
  	<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;calculate_stress();
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::calculate_tangent() {
  	Real mu = E / (2 * (1 + nu));
  	Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  
  	Real detF = F.det();
  
  	C_mat.resize(6, 6);
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; 3; ++i) {
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j &lt; 3; ++j) {
  			<B><FONT COLOR="#A020F0">if</FONT></B> (i == j) {
  				C_mat(i, j) = 2 * mu + lambda;
  				C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
  			} <B><FONT COLOR="#A020F0">else</FONT></B> {
  				C_mat(i, j) = lambda * detF * detF;
  			}
  		}
  	}
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::calculate_stress() {
  
    <B><FONT COLOR="#228B22">double</FONT></B> mu = E / (2.0 * (1.0 + nu));
  	<B><FONT COLOR="#228B22">double</FONT></B> lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  
  	Real detF = F.det();
  	RealTensor Ft = F.transpose();
  
  	RealTensor C = Ft * F;
  	RealTensor b = F * Ft;
  	RealTensor identity;
  	identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
  	RealTensor invC = inv(C);
  
  	S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);
  
  	tau = (F * S) * Ft;
  	sigma = 1/detF * tau;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::get_residual(DenseVector&lt;Real&gt; &amp; residuum, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; i) {
  	B_L.resize(3, 6);
  	DenseVector&lt;Real&gt; sigma_voigt(6);
  
  	<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;build_b_0_mat(i, B_L);
  
  	tensor_to_voigt(sigma, sigma_voigt);
  
  	B_L.vector_mult(residuum, sigma_voigt);
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::tensor_to_voigt(<B><FONT COLOR="#228B22">const</FONT></B> RealTensor &amp;tensor, DenseVector&lt;Real&gt; &amp;vec) {
    vec(0) = tensor(0, 0);
    vec(1) = tensor(1, 1);
    vec(2) = tensor(2, 2);
    vec(3) = tensor(0, 1);
    vec(4) = tensor(1, 2);
    vec(5) = tensor(0, 2);
  
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::get_linearized_stiffness(DenseMatrix&lt;Real&gt; &amp; stiffness, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; i, <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> &amp; j) {
  	stiffness.resize(3, 3);
  
  	<B><FONT COLOR="#228B22">double</FONT></B> G_IK = (sigma * dphi[i][current_qp]) * dphi[j][current_qp];
  	stiffness(0, 0) += G_IK;
  	stiffness(1, 1) += G_IK;
  	stiffness(2, 2) += G_IK;
  
  	B_L.resize(3, 6);
  	<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;build_b_0_mat(i, B_L);
  	B_K.resize(3, 6);
  	<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;build_b_0_mat(j, B_K);
  
  	B_L.right_multiply(C_mat);
  	B_L.right_multiply_transpose(B_K);
  	B_L *= 1/F.det();
  
  	stiffness += B_L;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> NonlinearNeoHookeCurrentConfig::build_b_0_mat(<B><FONT COLOR="#228B22">int</FONT></B> i, DenseMatrix&lt;Real&gt;&amp; b_0_mat) {
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ii = 0; ii &lt; 3; ++ii) {
  		b_0_mat(ii, ii) = dphi[i][current_qp](ii);
  	}
  	b_0_mat(0, 3) = dphi[i][current_qp](1);
  	b_0_mat(1, 3) = dphi[i][current_qp](0);
  	b_0_mat(1, 4) = dphi[i][current_qp](2);
  	b_0_mat(2, 4) = dphi[i][current_qp](1);
  	b_0_mat(0, 5) = dphi[i][current_qp](2);
  	b_0_mat(2, 5) = dphi[i][current_qp](0);
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file solid_system.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/diff_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_base.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fem_context.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newton_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/steady_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;nonlinear_neohooke_cc.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;solid_system.h&quot;</FONT></B>
  
  #ifdef __SUNPRO_CC
    #define NAN (1.0/0.0)
  #endif
  
  using namespace libMesh;
  
  <B><FONT COLOR="#5F9EA0">SolidSystem</FONT></B>::SolidSystem(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name,
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number) :
      FEMSystem(es, name, number) {
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_solver = AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> SteadySolver(*<B><FONT COLOR="#A020F0">this</FONT></B>));
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> SolidSystem::save_initial_mesh() {
    System &amp; aux_sys = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().get_system(<B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_node_iterator nd_end =
        <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().local_nodes_end();
    <B><FONT COLOR="#A020F0">for</FONT></B> (MeshBase::const_node_iterator nd = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().local_nodes_begin();
        nd != nd_end; ++nd) {
      <B><FONT COLOR="#228B22">const</FONT></B> Node *node = *nd;
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d = 0; d &lt; dim; ++d) {
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> source_dof = node-&gt;dof_number(<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;number(), var[d], 0);
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dest_dof = node-&gt;dof_number(aux_sys.number(), undefo_var[d],
            0);
        Real value = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;current_local_solution-&gt;el(source_dof);
        aux_sys.current_local_solution-&gt;set(dest_dof, value);
      }
    }
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> SolidSystem::init_data() {
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    Order order = (*(<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().elements_begin()))-&gt;default_order();
  
    var[0] = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;x&quot;</FONT></B>, order);
    var[1] = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;y&quot;</FONT></B>, order);
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      var[2] = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable(<B><FONT COLOR="#BC8F8F">&quot;z&quot;</FONT></B>, order);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      var[2] = var[1];
  
    System&amp; aux_sys = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().get_system(<B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
    undefo_var[0] = aux_sys.add_variable(<B><FONT COLOR="#BC8F8F">&quot;undefo_x&quot;</FONT></B>, order);
    undefo_var[1] = aux_sys.add_variable(<B><FONT COLOR="#BC8F8F">&quot;undefo_y&quot;</FONT></B>, order);
    undefo_var[2] = aux_sys.add_variable(<B><FONT COLOR="#BC8F8F">&quot;undefo_z&quot;</FONT></B>, order);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;deltat = args(<B><FONT COLOR="#BC8F8F">&quot;schedule/dt&quot;</FONT></B>, 0.2);
  
    <B><FONT COLOR="#5F9EA0">FEMSystem</FONT></B>::init_data();
  
  
    set_mesh_system(<B><FONT COLOR="#A020F0">this</FONT></B>);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;set_mesh_x_var(var[0]);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;set_mesh_y_var(var[1]);
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;set_mesh_z_var(var[2]);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;mesh_position_get();
  
    <B><FONT COLOR="#5F9EA0">System</FONT></B>::reinit();
  
    DiffSolver &amp;solver = *(<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_solver-&gt;diff_solver().get());
    solver.quiet = args(<B><FONT COLOR="#BC8F8F">&quot;solver/quiet&quot;</FONT></B>, false);
    solver.max_nonlinear_iterations = args(
        <B><FONT COLOR="#BC8F8F">&quot;solver/nonlinear/max_nonlinear_iterations&quot;</FONT></B>, 100);
    solver.relative_step_tolerance = args(
        <B><FONT COLOR="#BC8F8F">&quot;solver/nonlinear/relative_step_tolerance&quot;</FONT></B>, 1.e-3);
    solver.relative_residual_tolerance = args(
        <B><FONT COLOR="#BC8F8F">&quot;solver/nonlinear/relative_residual_tolerance&quot;</FONT></B>, 1.e-8);
    solver.absolute_residual_tolerance = args(
        <B><FONT COLOR="#BC8F8F">&quot;solver/nonlinear/absolute_residual_tolerance&quot;</FONT></B>, 1.e-8);
    solver.verbose = !args(<B><FONT COLOR="#BC8F8F">&quot;solver/quiet&quot;</FONT></B>, false);
  
    ((NewtonSolver&amp;) solver).require_residual_reduction = args(
        <B><FONT COLOR="#BC8F8F">&quot;solver/nonlinear/require_reduction&quot;</FONT></B>, false);
  
    solver.max_linear_iterations = args(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 50000);
    solver.initial_linear_tolerance = args(<B><FONT COLOR="#BC8F8F">&quot;initial_linear_tolerance&quot;</FONT></B>, 1.e-3);
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> SolidSystem::update() {
    <B><FONT COLOR="#5F9EA0">System</FONT></B>::update();
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;mesh_position_set();
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> SolidSystem::init_context(DiffContext &amp;context) {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
    c.element_fe_var[var[0]]-&gt;get_JxW();
    c.element_fe_var[var[0]]-&gt;get_phi();
    c.element_fe_var[var[0]]-&gt;get_dphi();
    c.element_fe_var[var[0]]-&gt;get_xyz();
  
    c.side_fe_var[var[0]]-&gt;get_JxW();
    c.side_fe_var[var[0]]-&gt;get_phi();
    c.side_fe_var[var[0]]-&gt;get_xyz();
  }
  
  <I><FONT COLOR="#B22222">/**
   * Compute contribution from internal forces in elem_residual and contribution from
   * linearized internal forces (stiffness matrix) in elem_jacobian.
   */</FONT></I>
  <B><FONT COLOR="#228B22">bool</FONT></B> SolidSystem::element_time_derivative(<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
      DiffContext &amp;context) {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = c.element_fe_var[var[0]]-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi =
        c.element_fe_var[var[0]]-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().mesh_dimension();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = c.dof_indices_var[var[0]].size();
    libmesh_assert(n_u_dofs == c.dof_indices_var[var[1]].size());
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3) {
      libmesh_assert(n_u_dofs == c.dof_indices_var[var[2]].size());
    }
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_qpoints = c.element_qrule-&gt;n_points();
  
    DenseMatrix&lt;Real&gt; stiff;
    DenseVector&lt;Real&gt; res;
    VectorValue&lt;Gradient&gt; grad_u;
  
    NonlinearNeoHookeCurrentConfig material(dphi, args);
  
    material.calculate_linearized_stiffness = request_jacobian;
  
    TransientExplicitSystem&amp; aux_system = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().get_system&lt;
        TransientExplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; undefo_index;
  
    <B><FONT COLOR="#228B22">bool</FONT></B> use_symmetry = args(<B><FONT COLOR="#BC8F8F">&quot;assembly/use_symmetry&quot;</FONT></B>, false);
  
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp = 0; qp != n_qpoints; qp++) {
      grad_u(0) = grad_u(1) = grad_u(2) = 0;
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d = 0; d &lt; dim; ++d) {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; u_undefo;
        aux_system.get_dof_map().dof_indices(c.elem, undefo_index, undefo_var[d]);
        aux_system.current_local_solution-&gt;get(undefo_index, u_undefo);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l = 0; l != n_u_dofs; l++)
          grad_u(d).add_scaled(dphi[l][qp], u_undefo[l]); <I><FONT COLOR="#B22222">// u_current(l)); // -
</FONT></I>      }
  
      material.init_for_qp(grad_u, qp);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; n_u_dofs; i++) {
        res.resize(dim);
        material.get_residual(res, i);
        res.scale(JxW[qp]);
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ii = 0; ii &lt; dim; ++ii) {
          c.elem_subresiduals[ii]-&gt;<B><FONT COLOR="#A020F0">operator</FONT></B> ()(i) += res(ii);
        }
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian &amp;&amp; c.elem_solution_derivative) {
          libmesh_assert(c.elem_solution_derivative == 1.0);
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = (use_symmetry ? i : 0); j &lt; n_u_dofs; j++) {
            material.get_linearized_stiffness(stiff, i, j);
            stiff.scale(JxW[qp]);
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ii = 0; ii &lt; dim; ++ii) {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> jj = 0; jj &lt; dim; ++jj) {
                c.elem_subjacobians[ii][jj]-&gt;<B><FONT COLOR="#A020F0">operator</FONT></B> ()(i, j) += stiff(ii, jj);
                <B><FONT COLOR="#A020F0">if</FONT></B> (use_symmetry &amp;&amp; i != j) {
                  c.elem_subjacobians[ii][jj]-&gt;<B><FONT COLOR="#A020F0">operator</FONT></B> ()(j, i) += stiff(jj, ii);
                }
              }
            }
          }
        }
      }
    } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  <B><FONT COLOR="#228B22">bool</FONT></B> SolidSystem::side_time_derivative(<B><FONT COLOR="#228B22">bool</FONT></B> request_jacobian,
      DiffContext &amp;context) {
    FEMContext &amp;c = libmesh_cast_ref&lt;FEMContext&amp;&gt;(context);
  
  
    Real ratio = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;progress&quot;</FONT></B>)
        + 0.001;
  
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> num_bc = args.vector_variable_size(<B><FONT COLOR="#BC8F8F">&quot;bc/displacement&quot;</FONT></B>);
    <B><FONT COLOR="#A020F0">if</FONT></B> (num_bc % 4 != 0) {
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::err
          &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR, Odd number of values in displacement boundary condition.\n&quot;</FONT></B>
          &lt;&lt; std::endl;
      libmesh_error();
    }
    num_bc /= 4;
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nbc = 0; nbc &lt; num_bc; nbc++) {
      <B><FONT COLOR="#228B22">short</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> positive_boundary_id = args(<B><FONT COLOR="#BC8F8F">&quot;bc/displacement&quot;</FONT></B>, 1, nbc * 4);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (!<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_mesh().boundary_info-&gt;has_boundary_id
  	  (c.elem,c.side,positive_boundary_id))
        <B><FONT COLOR="#A020F0">continue</FONT></B>;
  
      Point diff_value;
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d = 0; d &lt; c.dim; ++d) {
        diff_value(d) = args(<B><FONT COLOR="#BC8F8F">&quot;bc/displacement&quot;</FONT></B>, NAN, nbc * 4 + 1 + d);
      }
      diff_value *= ratio;
  
      Real penalty_number = args(<B><FONT COLOR="#BC8F8F">&quot;bc/displacement_penalty&quot;</FONT></B>, 1e7);
  
      FEBase * fe = c.side_fe_var[var[0]];
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp; phi = fe-&gt;get_phi();
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
      <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; coords = fe-&gt;get_xyz();
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_x_dofs = c.dof_indices_var[<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;var[0]].size();
  
      <B><FONT COLOR="#228B22">const</FONT></B> System &amp; auxsys = <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;get_equation_systems().get_system(
          <B><FONT COLOR="#BC8F8F">&quot;auxiliary&quot;</FONT></B>);
      <B><FONT COLOR="#228B22">const</FONT></B> DofMap &amp; auxmap = auxsys.get_dof_map();
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; undefo_dofs[3];
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d = 0; d &lt; c.dim; ++d) {
        auxmap.dof_indices(c.elem, undefo_dofs[d], undefo_var[d]);
      }
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp = 0; qp &lt; c.side_qrule-&gt;n_points(); ++qp) {
        Point orig_point;
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; n_x_dofs; ++i) {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d = 0; d &lt; c.dim; ++d) {
            Real orig_val = auxsys.current_solution(undefo_dofs[d][i]);
            orig_point(d) += phi[i][qp] * orig_val;
          }
        }
  
        Point diff = coords[qp] - orig_point - diff_value;
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = 0; i &lt; n_x_dofs; ++i) {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d1 = 0; d1 &lt; c.dim; ++d1) {
            <B><FONT COLOR="#A020F0">if</FONT></B> (libmesh_isnan(diff(d1)))
              <B><FONT COLOR="#A020F0">continue</FONT></B>;
            Real val = JxW[qp] * phi[i][qp] * diff(d1) * penalty_number;
            c.elem_subresiduals[var[d1]]-&gt;<B><FONT COLOR="#A020F0">operator</FONT></B> ()(i) += val;
          }
          <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian) {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j &lt; n_x_dofs; ++j) {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> d1 = 0; d1 &lt; c.dim; ++d1) {
                <B><FONT COLOR="#A020F0">if</FONT></B> (libmesh_isnan(diff(d1)))
                  <B><FONT COLOR="#A020F0">continue</FONT></B>;
                Real val = JxW[qp] * phi[i][qp] * phi[j][qp] * penalty_number;
                c.elem_subjacobians[var[d1]][var[d1]]-&gt;<B><FONT COLOR="#A020F0">operator</FONT></B> ()(i, j) += val;
              }
            }
          }
        }
      }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example fem_system_ex2:
*  mpirun -np 2 example-dbg  
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/diff_physics.h, line 371, compiled Feb  1 2013 at 09:22:48 ***
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/vtk_io.h, line 178, compiled Feb  1 2013 at 09:22:42 ***
===== Time Step    0 (  0.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 3691.41
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 0.35
  Nonlinear solver current_residual 0.35 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 0, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    1 ( 20.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 738516.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 72.25
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.04
Assembling the System
Nonlinear Residual: 72.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 1.56
  Nonlinear solver current_residual 1.56 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    2 ( 40.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 788389.66
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 84.42
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.05
Assembling the System
Nonlinear Residual: 84.42
Linear solve starting, tolerance 0.00
Linear solve finished, step 12, residual 0.00
Trying full Newton step
  Current Residual: 3.97
  Nonlinear solver current_residual 3.97 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    3 ( 60.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 844777.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 124.74
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.06
Assembling the System
Nonlinear Residual: 124.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 8.64
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.02
Assembling the System
Nonlinear Residual: 8.64
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 0.03
  Nonlinear solver current_residual 0.03 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    4 ( 80.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 909688.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 256.39
  Nonlinear step: |du|/|u| = 0.07, |du| = 1.07
Assembling the System
Nonlinear Residual: 256.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 20.00
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.03
Assembling the System
Nonlinear Residual: 20.00
Linear solve starting, tolerance 0.00
Linear solve finished, step 17, residual 0.00
Trying full Newton step
  Current Residual: 0.12
  Nonlinear solver current_residual 0.12 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
Finished calculations

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh10DiffSolverE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh10FEAbstractE reference count information:
|  Creations:    534
|  Destructions: 534
| N7libMesh10Parameters5ValueE reference count information:
|  Creations:    5
|  Destructions: 5
| N7libMesh10TimeSolverE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh12LinearSolverIdEE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh12SparseMatrixIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh13NumericVectorIdEE reference count information:
|  Creations:    75
|  Destructions: 75
| N7libMesh15EquationSystemsE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh4ElemE reference count information:
|  Creations:    2465
|  Destructions: 2465
| N7libMesh4NodeE reference count information:
|  Creations:    125
|  Destructions: 125
| N7libMesh5QBaseE reference count information:
|  Creations:    465
|  Destructions: 465
| N7libMesh6DofMapE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh6SystemE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh9DofObjectE reference count information:
|  Creations:    4480
|  Destructions: 4480
 ---------------------------------------------------------------------------- 

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:31:49 2013                    |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=6.03736, Active time=5.89942                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     12        0.0633      0.005278    0.1407      0.011723    1.07     2.38     |
|   build_sparsity()                 6         0.0768      0.012795    0.1419      0.023651    1.30     2.41     |
|   create_dof_constraints()         12        0.0033      0.000272    0.0033      0.000272    0.06     0.06     |
|   distribute_dofs()                12        0.2208      0.018399    0.8698      0.072484    3.74     14.74    |
|   dof_indices()                    29856     2.4221      0.000081    2.4221      0.000081    41.06    41.06    |
|   old_dof_indices()                2560      0.2635      0.000103    0.2635      0.000103    4.47     4.47     |
|   prepare_send_list()              12        0.0020      0.000169    0.0020      0.000169    0.03     0.03     |
|   reinit()                         12        0.2184      0.018203    0.2184      0.018203    3.70     3.70     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          5         0.0055      0.001109    0.0743      0.014862    0.09     1.26     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        6272      0.5041      0.000080    0.5041      0.000080    8.55     8.55     |
|   init_shape_functions()           2242      0.0567      0.000025    0.0567      0.000025    0.96     0.96     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       11        0.9274      0.084313    2.1994      0.199943    15.72    37.28    |
|   assembly(get_residual)           11        0.1628      0.014801    1.4332      0.130294    2.76     24.29    |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             290       0.0045      0.000016    0.0045      0.000016    0.08     0.08     |
|   compute_face_map()               2112      0.0289      0.000014    0.0289      0.000014    0.49     0.49     |
|   compute_map()                    5982      0.2241      0.000037    0.2241      0.000037    3.80     3.80     |
|   init_face_shape_functions()      44        0.0012      0.000028    0.0012      0.000028    0.02     0.02     |
|   init_reference_to_physical_map() 2242      0.0796      0.000035    0.0796      0.000035    1.35     1.35     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   init()                           5         0.0013      0.000266    0.0013      0.000266    0.02     0.02     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       5         0.0003      0.000068    0.0007      0.000137    0.01     0.01     |
|   find_neighbors()                 1         0.0057      0.005743    0.0059      0.005915    0.10     0.10     |
|   renumber_nodes_and_elem()        7         0.0007      0.000099    0.0007      0.000099    0.01     0.01     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0024      0.001188    0.0024      0.001188    0.04     0.04     |
|   find_global_indices()            2         0.0011      0.000552    0.0053      0.002652    0.02     0.09     |
|   parallel_sort()                  2         0.0010      0.000523    0.0012      0.000619    0.02     0.02     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         5         0.0167      0.003331    0.0913      0.018268    0.28     1.55     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              5         0.0002      0.000036    0.0002      0.000040    0.00     0.00     |
|   _refine_elements()               5         0.0004      0.000082    0.0004      0.000087    0.01     0.01     |
|   make_coarsening_compatible()     10        0.0004      0.000039    0.0004      0.000039    0.01     0.01     |
|   make_refinement_compatible()     10        0.0000      0.000003    0.0001      0.000007    0.00     0.00     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0020      0.002030    0.0020      0.002030    0.03     0.03     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0045      0.004475    0.0068      0.006819    0.08     0.12     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          5         0.0639      0.012775    4.1008      0.820158    1.08     69.51    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      40        0.0005      0.000011    0.0005      0.000013    0.01     0.01     |
|   max(bool)                        15907     0.0504      0.000003    0.0504      0.000003    0.85     0.85     |
|   max(scalar)                      34154     0.1126      0.000003    0.1126      0.000003    1.91     1.91     |
|   max(vector)                      7131      0.0519      0.000007    0.1268      0.000018    0.88     2.15     |
|   max(vector<bool>)                25        0.0011      0.000043    0.0014      0.000057    0.02     0.02     |
|   min(bool)                        15316     0.0503      0.000003    0.0503      0.000003    0.85     0.85     |
|   min(scalar)                      27286     0.1475      0.000005    0.1475      0.000005    2.50     2.50     |
|   min(vector)                      7131      0.0522      0.000007    0.1287      0.000018    0.88     2.18     |
|   probe()                          152       0.0029      0.000019    0.0029      0.000019    0.05     0.05     |
|   receive()                        152       0.0012      0.000008    0.0042      0.000028    0.02     0.07     |
|   send()                           152       0.0007      0.000004    0.0007      0.000004    0.01     0.01     |
|   send_receive()                   156       0.0017      0.000011    0.0070      0.000045    0.03     0.12     |
|   sum()                            67        0.0008      0.000011    0.0015      0.000023    0.01     0.03     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           152       0.0003      0.000002    0.0003      0.000002    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0018      0.001837    0.0047      0.004694    0.03     0.08     |
|   set_parent_processor_ids()       1         0.0002      0.000220    0.0002      0.000220    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          11        0.0306      0.002786    0.0306      0.002786    0.52     0.52     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       20        0.0151      0.000756    0.2823      0.014114    0.26     4.79     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 20        0.0119      0.000597    0.4275      0.021376    0.20     7.25     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            159633    5.8994                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example fem_system_ex2:
*  mpirun -np 2 example-dbg  
***************************************************************
***************************************************************
* Running Example fem_system_ex2:
*  mpirun -np 2 example-devel  
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/diff_physics.h, line 371, compiled Feb  1 2013 at 09:23:13 ***
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/vtk_io.h, line 178, compiled Feb  1 2013 at 09:23:01 ***
===== Time Step    0 (  0.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 3691.41
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 0.35
  Nonlinear solver current_residual 0.35 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 0, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    1 ( 20.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 738516.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 72.25
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.04
Assembling the System
Nonlinear Residual: 72.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 1.56
  Nonlinear solver current_residual 1.56 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    2 ( 40.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 788389.66
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 84.42
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.05
Assembling the System
Nonlinear Residual: 84.42
Linear solve starting, tolerance 0.00
Linear solve finished, step 12, residual 0.00
Trying full Newton step
  Current Residual: 3.97
  Nonlinear solver current_residual 3.97 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    3 ( 60.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 844777.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 124.74
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.06
Assembling the System
Nonlinear Residual: 124.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 8.64
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.02
Assembling the System
Nonlinear Residual: 8.64
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 0.03
  Nonlinear solver current_residual 0.03 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    4 ( 80.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 909688.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 256.39
  Nonlinear step: |du|/|u| = 0.07, |du| = 1.07
Assembling the System
Nonlinear Residual: 256.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 20.00
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.03
Assembling the System
Nonlinear Residual: 20.00
Linear solve starting, tolerance 0.00
Linear solve finished, step 17, residual 0.00
Trying full Newton step
  Current Residual: 0.12
  Nonlinear solver current_residual 0.12 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
Finished calculations

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:31:50 2013                    |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.891264, Active time=0.832516                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     12        0.0045      0.000377    0.0093      0.000779    0.54     1.12     |
|   build_sparsity()                 6         0.0055      0.000923    0.0099      0.001648    0.66     1.19     |
|   create_dof_constraints()         12        0.0007      0.000059    0.0007      0.000059    0.08     0.08     |
|   distribute_dofs()                12        0.0031      0.000255    0.0086      0.000720    0.37     1.04     |
|   dof_indices()                    29856     0.1645      0.000006    0.1645      0.000006    19.75    19.75    |
|   old_dof_indices()                2560      0.0172      0.000007    0.0172      0.000007    2.07     2.07     |
|   prepare_send_list()              12        0.0003      0.000022    0.0003      0.000022    0.03     0.03     |
|   reinit()                         12        0.0036      0.000298    0.0036      0.000298    0.43     0.43     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          5         0.0019      0.000379    0.0078      0.001552    0.23     0.93     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        6272      0.0703      0.000011    0.0703      0.000011    8.44     8.44     |
|   init_shape_functions()           2242      0.0129      0.000006    0.0129      0.000006    1.54     1.54     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       11        0.3165      0.028769    0.4554      0.041398    38.01    54.70    |
|   assembly(get_residual)           11        0.0470      0.004271    0.1855      0.016864    5.64     22.28    |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             290       0.0012      0.000004    0.0012      0.000004    0.15     0.15     |
|   compute_face_map()               2112      0.0091      0.000004    0.0091      0.000004    1.09     1.09     |
|   compute_map()                    5982      0.0679      0.000011    0.0679      0.000011    8.15     8.15     |
|   init_face_shape_functions()      44        0.0003      0.000007    0.0003      0.000007    0.04     0.04     |
|   init_reference_to_physical_map() 2242      0.0179      0.000008    0.0179      0.000008    2.15     2.15     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   init()                           5         0.0004      0.000074    0.0004      0.000074    0.04     0.04     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       5         0.0001      0.000015    0.0002      0.000037    0.01     0.02     |
|   find_neighbors()                 1         0.0004      0.000379    0.0004      0.000442    0.05     0.05     |
|   renumber_nodes_and_elem()        7         0.0002      0.000023    0.0002      0.000023    0.02     0.02     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0006      0.000284    0.0006      0.000284    0.07     0.07     |
|   find_global_indices()            2         0.0001      0.000072    0.0011      0.000562    0.02     0.14     |
|   parallel_sort()                  2         0.0002      0.000082    0.0002      0.000115    0.02     0.03     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         5         0.0112      0.002241    0.0193      0.003851    1.35     2.31     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              5         0.0000      0.000010    0.0001      0.000012    0.01     0.01     |
|   _refine_elements()               5         0.0001      0.000024    0.0002      0.000031    0.01     0.02     |
|   make_coarsening_compatible()     10        0.0001      0.000015    0.0001      0.000015    0.02     0.02     |
|   make_refinement_compatible()     10        0.0000      0.000002    0.0000      0.000004    0.00     0.01     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000187    0.0002      0.000187    0.02     0.02     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0008      0.000788    0.0013      0.001289    0.09     0.15     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          5         0.0137      0.002734    0.7298      0.145954    1.64     87.66    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      40        0.0002      0.000006    0.0003      0.000006    0.03     0.03     |
|   max(bool)                        21        0.0001      0.000006    0.0001      0.000006    0.02     0.02     |
|   max(scalar)                      1185      0.0031      0.000003    0.0031      0.000003    0.37     0.37     |
|   max(vector)                      288       0.0014      0.000005    0.0037      0.000013    0.17     0.44     |
|   min(bool)                        1478      0.0037      0.000003    0.0037      0.000003    0.45     0.45     |
|   min(scalar)                      1172      0.0062      0.000005    0.0062      0.000005    0.74     0.74     |
|   min(vector)                      288       0.0015      0.000005    0.0038      0.000013    0.18     0.45     |
|   probe()                          152       0.0010      0.000007    0.0010      0.000007    0.12     0.12     |
|   receive()                        152       0.0006      0.000004    0.0016      0.000011    0.07     0.19     |
|   send()                           152       0.0005      0.000003    0.0005      0.000003    0.06     0.06     |
|   send_receive()                   156       0.0008      0.000005    0.0031      0.000020    0.09     0.37     |
|   sum()                            67        0.0005      0.000007    0.0006      0.000010    0.05     0.08     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           152       0.0001      0.000001    0.0001      0.000001    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000125    0.0003      0.000277    0.02     0.03     |
|   set_parent_processor_ids()       1         0.0000      0.000026    0.0000      0.000026    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          11        0.0296      0.002688    0.0296      0.002688    3.55     3.55     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       20        0.0052      0.000260    0.0243      0.001216    0.62     2.92     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 20        0.0057      0.000284    0.0400      0.001998    0.68     4.80     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            57115     0.8325                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example fem_system_ex2:
*  mpirun -np 2 example-devel  
***************************************************************
***************************************************************
* Running Example fem_system_ex2:
*  mpirun -np 2 example-opt  
***************************************************************
 
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/diff_physics.h, line 371, compiled Feb  1 2013 at 09:22:30 ***
*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/vtk_io.h, line 178, compiled Feb  1 2013 at 09:22:23 ***
===== Time Step    0 (  0.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 3691.41
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 0.35
  Nonlinear solver current_residual 0.35 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 0, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    1 ( 20.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 738516.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 8, residual 0.00
Trying full Newton step
  Current Residual: 72.25
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.04
Assembling the System
Nonlinear Residual: 72.25
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 1.56
  Nonlinear solver current_residual 1.56 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    2 ( 40.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 788389.66
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 84.42
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.05
Assembling the System
Nonlinear Residual: 84.42
Linear solve starting, tolerance 0.00
Linear solve finished, step 12, residual 0.00
Trying full Newton step
  Current Residual: 3.97
  Nonlinear solver current_residual 3.97 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 1, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    3 ( 60.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 844777.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 124.74
  Nonlinear step: |du|/|u| = 0.06, |du| = 1.06
Assembling the System
Nonlinear Residual: 124.74
Linear solve starting, tolerance 0.00
Linear solve finished, step 13, residual 0.00
Trying full Newton step
  Current Residual: 8.64
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.02
Assembling the System
Nonlinear Residual: 8.64
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 0.03
  Nonlinear solver current_residual 0.03 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
===== Time Step    4 ( 80.00%), time =    0.00 =====
Solving Solid
Assembling the System
Nonlinear Residual: 909688.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 9, residual 0.00
Trying full Newton step
  Current Residual: 256.39
  Nonlinear step: |du|/|u| = 0.07, |du| = 1.07
Assembling the System
Nonlinear Residual: 256.39
Linear solve starting, tolerance 0.00
Linear solve finished, step 14, residual 0.00
Trying full Newton step
  Current Residual: 20.00
  Nonlinear step: |du|/|u| = 0.00, |du| = 0.03
Assembling the System
Nonlinear Residual: 20.00
Linear solve starting, tolerance 0.00
Linear solve finished, step 17, residual 0.00
Trying full Newton step
  Current Residual: 0.12
  Nonlinear solver current_residual 0.12 > 0.00
  Nonlinear solver relative residual 0.00 > 0.00
  Nonlinear solver converged, step 2, relative step size 0.00 < 0.00
Doing a reinit of the equation systems
Advancing to next step
Finished calculations

 -------------------------------------------------------------
| Processor id:   0                                           |
| Num Processors: 2                                           |
| Time:           Fri Feb  1 09:31:54 2013                    |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.579611, Active time=0.511288                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     12        0.0012      0.000099    0.0018      0.000150    0.23     0.35     |
|   build_sparsity()                 6         0.0051      0.000847    0.0071      0.001187    0.99     1.39     |
|   create_dof_constraints()         12        0.0006      0.000052    0.0006      0.000052    0.12     0.12     |
|   distribute_dofs()                12        0.0012      0.000098    0.0037      0.000312    0.23     0.73     |
|   dof_indices()                    29856     0.0237      0.000001    0.0237      0.000001    4.63     4.63     |
|   old_dof_indices()                2560      0.0023      0.000001    0.0023      0.000001    0.45     0.45     |
|   prepare_send_list()              12        0.0003      0.000023    0.0003      0.000023    0.05     0.05     |
|   reinit()                         12        0.0018      0.000153    0.0018      0.000153    0.36     0.36     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          5         0.0015      0.000305    0.0026      0.000519    0.30     0.51     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        3520      0.0176      0.000005    0.0176      0.000005    3.44     3.44     |
|   init_shape_functions()           2156      0.0108      0.000005    0.0108      0.000005    2.11     2.11     |
|                                                                                                                |
| FEMSystem                                                                                                      |
|   assembly()                       11        0.2861      0.026005    0.3298      0.029985    55.95    64.51    |
|   assembly(get_residual)           11        0.0385      0.003498    0.0798      0.007258    7.52     15.62    |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             160       0.0004      0.000002    0.0004      0.000002    0.08     0.08     |
|   compute_face_map()               2112      0.0039      0.000002    0.0039      0.000002    0.77     0.77     |
|   compute_map()                    3360      0.0111      0.000003    0.0111      0.000003    2.17     2.17     |
|   init_face_shape_functions()      44        0.0003      0.000007    0.0003      0.000007    0.06     0.06     |
|   init_reference_to_physical_map() 2156      0.0158      0.000007    0.0158      0.000007    3.09     3.09     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   init()                           5         0.0002      0.000047    0.0002      0.000047    0.05     0.05     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       5         0.0000      0.000009    0.0001      0.000021    0.01     0.02     |
|   find_neighbors()                 1         0.0003      0.000311    0.0015      0.001484    0.06     0.29     |
|   renumber_nodes_and_elem()        7         0.0001      0.000013    0.0001      0.000013    0.02     0.02     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0006      0.000302    0.0006      0.000302    0.12     0.12     |
|   find_global_indices()            2         0.0001      0.000063    0.0016      0.000800    0.02     0.31     |
|   parallel_sort()                  2         0.0002      0.000100    0.0006      0.000319    0.04     0.12     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         5         0.0109      0.002182    0.0135      0.002701    2.13     2.64     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              5         0.0000      0.000005    0.0000      0.000007    0.00     0.01     |
|   _refine_elements()               5         0.0001      0.000012    0.0001      0.000019    0.01     0.02     |
|   make_coarsening_compatible()     5         0.0000      0.000007    0.0000      0.000007    0.01     0.01     |
|   make_refinement_compatible()     5         0.0000      0.000002    0.0000      0.000005    0.00     0.01     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000200    0.0002      0.000200    0.04     0.04     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0006      0.000613    0.0011      0.001053    0.12     0.21     |
|                                                                                                                |
| NewtonSolver                                                                                                   |
|   solve()                          5         0.0285      0.005695    0.4724      0.094477    5.57     92.39    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      40        0.0003      0.000008    0.0003      0.000008    0.06     0.06     |
|   max(bool)                        21        0.0001      0.000004    0.0001      0.000004    0.02     0.02     |
|   max(scalar)                      13        0.0012      0.000093    0.0012      0.000093    0.24     0.24     |
|   max(vector)                      2         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   min(bool)                        10        0.0000      0.000004    0.0000      0.000004    0.01     0.01     |
|   min(vector)                      2         0.0001      0.000055    0.0001      0.000055    0.02     0.02     |
|   probe()                          152       0.0025      0.000017    0.0025      0.000017    0.50     0.50     |
|   receive()                        152       0.0005      0.000003    0.0031      0.000020    0.10     0.61     |
|   send()                           152       0.0003      0.000002    0.0003      0.000002    0.07     0.07     |
|   send_receive()                   156       0.0007      0.000005    0.0044      0.000028    0.14     0.86     |
|   sum()                            67        0.0009      0.000013    0.0009      0.000013    0.17     0.17     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           152       0.0001      0.000001    0.0001      0.000001    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000072    0.0001      0.000105    0.01     0.02     |
|   set_parent_processor_ids()       1         0.0000      0.000026    0.0000      0.000026    0.01     0.01     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          11        0.0309      0.002809    0.0309      0.002809    6.04     6.04     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       20        0.0044      0.000218    0.0077      0.000386    0.85     1.51     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 20        0.0049      0.000247    0.0137      0.000685    0.97     2.68     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            47045     0.5113                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example fem_system_ex2:
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
