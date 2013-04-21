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
<br><br><br> <h1> The source file biharmonic.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __biharmonic_h__
        #define __biharmonic_h__
        
        #include "libmesh/equation_systems.h"
        #include "libmesh/serial_mesh.h"
        #include "libmesh/exodusII_io.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::EquationSystems;
        using libMesh::ExodusII_IO;
        using libMesh::MeshRefinement;
        using libMesh::Point;
        using libMesh::Real;
        using libMesh::UnstructuredMesh;
        
</pre>
</div>
<div class = "comment">
libmesh_error() and libmesh_assert() macros with a message
</div>

<div class ="fragment">
<pre>
        #define ERROR(message)                                                                                         \
          do {                                                                                                         \
            libMesh::err &lt;&lt; "Error: " &lt;&lt; message &lt;&lt; "\n";                                                              \
            libmesh_error();							                                       \
          } while(0)
        
        #define ASSERT(asserted, message) \
          do {                    \
            if(!(asserted)) {     \
              libMesh::err &lt;&lt; "Assertion '" #asserted "' violated: " #message; \
              libmesh_error();    \
            }                     \
          } while(0)
        
        
        /**
         * The Biharmonic class encapsulates most of the data structures
         * necessary to calculate the biharmonic residual and Jacobian,
         * auxiliary quantities, to take a timestep, and to output the state --
         * biharmonic solution and vectors of auxiliary quantities.
         *
         * The main reason for this design is to have a data structure that
         * has all of the necessary data in one place, where all of the
         * calculation subroutines can access these data. Currently these data
         * are split up among several interdependent objects with no clear
         * hierarchy between them: mesh, equation system, equation system
         * bundle, residual/Jacobian calculator.
         *
         * Since no object contains all others and the data are distributed
         * among many objects, the natural control and data flow resides outside
         * of these objects and is typically implemented in main().  We,
         * however, would like to split the calculation into natural chunks --
         * subroutines -- while retaining these subroutines access to the common
         * necessary data -- biharmonic parameters, mesh and time interval
         * sizes, etc. Thus, class Biharmonic.  Finally, making Biharmonic
         * inherit from EquationSystems makes it possible to include it in the
         * most common callbacks that do not pass back a user context, but only
         * an EquationSystems object.
         */
        class Biharmonic : public EquationSystems
        {
        public:
</pre>
</div>
<div class = "comment">
Initial state enumeration
</div>

<div class ="fragment">
<pre>
          enum InitialStateEnum {STRIP = 0,
        			 ROD   = 1,
        			 BALL  = 2};
        
</pre>
</div>
<div class = "comment">
Free energy enumeration
</div>

<div class ="fragment">
<pre>
          enum FreeEnergyEnum {DOUBLE_WELL         = 1,
        		       DOUBLE_OBSTACLE     = 2,
        		       LOG_DOUBLE_WELL     = 3,
        		       LOG_DOUBLE_OBSTACLE = 4};
        
          /**
           * Static creation/destruction routines.  FIXME - this looks like
           * object-oriented C, can we get rid of it?
           */
          static void Create(Biharmonic** b, const Parallel::Communicator &comm);
          static void Destroy(Biharmonic** b);
        
        
          /**
           * Constructor retrieves command-line options, setting  defaults, if necessary.
           * It then builds the mesh using these options, then the equations systems around it,
           * and, finally, sets up the output.
           * We recommend that this be used through the factory Create function, which allocates
           * the mesh. In that case don't forget to call Destroy at the end, to free the mesh up.
           */
          Biharmonic(UnstructuredMesh* m);
        
        
          /**
           * Destructor
           */
          ~Biharmonic()
          {
</pre>
</div>
<div class = "comment">
delete _meshRefinement;
</div>

<div class ="fragment">
<pre>
          }
        
        
</pre>
</div>
<div class = "comment">
Misc. getters
</div>

<div class ="fragment">
<pre>
          bool verbose()         { return _verbose; }
          Real dt0()             { return _dt0; }
          Real dt()              { return _dt; }
        
        
</pre>
</div>
<div class = "comment">
Public interface functions
</div>

<div class ="fragment">
<pre>
          void viewParameters();
          void init();
          void step(const Real& dt = -1.0);
          void output(int timestep, const Real& t, Real& o_t, bool force = false);
          void run();
        
        private:
          unsigned int  _dim, _N;
          Real _kappa, _theta, _theta_c;
          Real _tol;
          bool _growth, _degenerate, _cahn_hillard, _netforce;
          FreeEnergyEnum  _energy;
          int _log_truncation;
          bool _verbose;
          InitialStateEnum  _initialState;
          Point _initialCenter;
          Real _initialWidth;
          Real _dt0, _dt, _t0, _T, _t1;
          Real _cnWeight;
</pre>
</div>
<div class = "comment">

<br><br></div>

<div class ="fragment">
<pre>
          std::string  _ofile_base, _ofile;
          ExodusII_IO* _exio;
          Real    _o_dt;
          int     _o_count;
</pre>
</div>
<div class = "comment">

<br><br></div>

<div class ="fragment">
<pre>
          friend class JR;
          class JR;       // forward
          UnstructuredMesh*                       _mesh;
          MeshRefinement*                         _meshRefinement;
          JR*                                     _jr;
        };
        
        
        
        
        
        
        
        #endif // __biharmonic_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file biharmonic_jr.h with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #ifndef __biharmonic_jr_h__
        #define __biharmonic_jr_h__
        
</pre>
</div>
<div class = "comment">
LibMesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/transient_system.h"
        #include "libmesh/nonlinear_solver.h"
        
        
</pre>
</div>
<div class = "comment">
Example includes
</div>

<div class ="fragment">
<pre>
        #include "biharmonic.h"
        
</pre>
</div>
<div class = "comment">
Bring in bits from the libMesh namespace.
Just the bits we're using, since this is a header.
</div>

<div class ="fragment">
<pre>
        using libMesh::EquationSystems;
        using libMesh::Gradient;
        using libMesh::NonlinearImplicitSystem;
        using libMesh::Number;
        using libMesh::NumericVector;
        using libMesh::Parameters;
        using libMesh::Point;
        using libMesh::SparseMatrix;
        using libMesh::System;
        using libMesh::TransientNonlinearImplicitSystem;
        
        
        /**
         * Biharmonic's friend class definition
         */
        class Biharmonic::JR : public TransientNonlinearImplicitSystem,
        		       public NonlinearImplicitSystem::ComputeResidualandJacobian,
        		       public NonlinearImplicitSystem::ComputeBounds,
        		       public System::Initialization
        {
        public:
          /**
           * Constructor.
           */
          JR(EquationSystems& eqSys, const std::string& name, const unsigned int number);
        
          void initialize();
        
          /**
           * Static functions to be used for initialization
           */
          static Number InitialDensityBall(const Point& p, const Parameters& parameters, const std::string&, const std::string&);
          static Number InitialDensityRod(const Point& p, const Parameters& parameters, const std::string&, const std::string&);
          static Number InitialDensityStrip(const Point& p, const Parameters& parameters, const std::string&, const std::string&);
          static Gradient InitialGradientZero(const Point&, const Parameters&, const std::string&, const std::string&);
        
          /**
           * The residual and Jacobian assembly function for the Biharmonic system.
           */
          void residual_and_jacobian(const NumericVector&lt;Number&gt;& u,
        			     NumericVector&lt;Number&gt;* R,
        			     SparseMatrix&lt;Number&gt;* J,
        			     NonlinearImplicitSystem&);
        
        
          /**
           * Function defining the bounds of the Biharmonic system.
           */
          void bounds(NumericVector&lt;Number&gt;& XL,
        	      NumericVector&lt;Number&gt;& XU,
        	      NonlinearImplicitSystem&);
        
        private:
          Biharmonic& _biharmonic;
        };
        
        
        #endif // __biharmonic_jr_h__
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file biharmonic.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh_generation.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/serial_mesh.h"
        
</pre>
</div>
<div class = "comment">
Example includes
</div>

<div class ="fragment">
<pre>
        #include "biharmonic.h"
        #include "biharmonic_jr.h"
        
        using namespace libMesh;
        
        void Biharmonic::Create(Biharmonic** b, const Parallel::Communicator &comm)
        {
</pre>
</div>
<div class = "comment">
ParallelMesh doesn't yet understand periodic BCs
</div>

<div class ="fragment">
<pre>
          SerialMesh* mesh = new SerialMesh(comm);
          Biharmonic *biharmonic = new Biharmonic(mesh);
          *b = biharmonic;
        }
        
        
        
        
        
        void Biharmonic::Destroy(Biharmonic** b)
        {
          Biharmonic* biharmonic = *b;
          UnstructuredMesh* mesh = biharmonic-&gt;_mesh;
          delete biharmonic;
          delete mesh;
          *b = NULL;
        }
        
        
        
        void Biharmonic::viewParameters()
        {
          libMesh::out &lt;&lt; "Biharmonic parameters:\n";
        
</pre>
</div>
<div class = "comment">
Print verbosity status
</div>

<div class ="fragment">
<pre>
          if (_verbose)
            libMesh::out &lt;&lt; "verbose mode is on\n";
          else
            libMesh::out &lt;&lt; "verbose mode is off\n";
        
</pre>
</div>
<div class = "comment">
Print parameters
</div>

<div class ="fragment">
<pre>
          libMesh::out &lt;&lt; "mesh dimension           = " &lt;&lt; _dim &lt;&lt;               "\n";
          libMesh::out &lt;&lt; "initial linear mesh size = " &lt;&lt; _N   &lt;&lt;               "\n";
          libMesh::out &lt;&lt; "kappa                    = " &lt;&lt; _kappa &lt;&lt;             "\n";
          libMesh::out &lt;&lt; "growth                   = " &lt;&lt; (int)_growth &lt;&lt;       "\n";
          libMesh::out &lt;&lt; "degenerate               = " &lt;&lt; (int)_degenerate &lt;&lt;   "\n";
          libMesh::out &lt;&lt; "Cahn-Hillard             = " &lt;&lt; (int)_cahn_hillard &lt;&lt; "\n";
          libMesh::out &lt;&lt; "netforce                 = " &lt;&lt; (int)_netforce &lt;&lt;     "\n";
          libMesh::out &lt;&lt; "energy                   = " &lt;&lt; _energy        &lt;&lt;     "\n";
          libMesh::out &lt;&lt; "tol                      = " &lt;&lt; _tol           &lt;&lt;     "\n";
          libMesh::out &lt;&lt; "theta                    = " &lt;&lt; _theta         &lt;&lt;     "\n";
          libMesh::out &lt;&lt; "theta_c                  = " &lt;&lt; _theta_c       &lt;&lt;     "\n";
          libMesh::out &lt;&lt; "log truncation           = " &lt;&lt; _log_truncation &lt;&lt;    "\n";
          libMesh::out &lt;&lt; "initial timestep size    = " &lt;&lt; _dt0            &lt;&lt;    "\n";
        
          if (_initialState == STRIP)
            libMesh::out &lt;&lt; "initial state:             strip\n";
        
          if (_initialState == ROD)
            libMesh::out &lt;&lt; "initial state:             rod\n";
        
          if (_initialState == BALL)
            libMesh::out &lt;&lt; "initial state:             ball\n";
        
          libMesh::out &lt;&lt; "initial state center     = " &lt;&lt; _initialCenter(0) &lt;&lt; "\n";
          libMesh::out &lt;&lt; "initial state width      = " &lt;&lt; _initialWidth &lt;&lt; "\n";
          libMesh::out &lt;&lt; "initial time (min_time)  = " &lt;&lt; _t0 &lt;&lt; "\n";
          libMesh::out &lt;&lt; "integration time         = " &lt;&lt; _T  &lt;&lt; "\n";
          libMesh::out &lt;&lt; "final time   (max_time)  = " &lt;&lt; _t1 &lt;&lt; "\n";
          libMesh::out &lt;&lt; "Crank-Nicholson weight   = " &lt;&lt; _cnWeight &lt;&lt; "\n";
          libMesh::out &lt;&lt; "Output timestep          = " &lt;&lt; _o_dt &lt;&lt; "\n";
          libMesh::out &lt;&lt; "Output filename base:      " &lt;&lt;  _ofile_base &lt;&lt; "\n";
        }
        
        
        
        
        void Biharmonic::init()
        {
          if(_verbose)
            libMesh::out &lt;&lt; "&gt;&gt;&gt; Initializing Biharmonic\n";
        
          _dt  =  0;
          _o_count = 0;
          this-&gt;EquationSystems::init();
        
          if(_verbose)
            libMesh::out &lt;&lt; "&lt;&lt;&lt; Initializing Biharmonic\n";
        }
        
        
        
        
        
        void Biharmonic::step(const Real& dt_in)
        {
</pre>
</div>
<div class = "comment">
We need to update the old solution vector.
The old solution vector will be the current solution vector from the
previous time step. We use vector assignment.  Only \p TransientSystems
(and systems derived from them) contain old solutions.
</div>

<div class ="fragment">
<pre>
          if (dt_in &lt; 0)
            _dt = _dt0;
          else
            _dt = dt_in;
        
          *(_jr-&gt;old_local_solution) = *(_jr-&gt;current_local_solution);
        
</pre>
</div>
<div class = "comment">
this will localize the current solution, resulting in a current_local_solution with correct ghost values
</div>

<div class ="fragment">
<pre>
          _jr-&gt;solve();
        }
        
        
        
        void Biharmonic::output(int timestep, const Real& t, Real& o_t, bool force)
        {
        #ifdef LIBMESH_HAVE_EXODUS_API
          if (!force && t - o_t &lt; _o_dt)
            return;
        
          ++_o_count;
        
          if (_verbose)
            libMesh::out &lt;&lt; "Writing state " &lt;&lt; timestep &lt;&lt; " at time " &lt;&lt; t &lt;&lt; " to file " &lt;&lt; _ofile &lt;&lt; "; output a total of " &lt;&lt; _o_count &lt;&lt; " states so far\n";
        
          _exio-&gt;write_timestep(_ofile, *this, timestep, t);
        
          if (!force)
            o_t = t;
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        }
        
        
        
        void Biharmonic::run()
        {
          Real t = _t0, o_t = 0.0;
          int timestep = 1;
        
</pre>
</div>
<div class = "comment">
Force-write the initial timestep
</div>

<div class ="fragment">
<pre>
          output(timestep,t,o_t,true);
        
          while (t &lt; _t1)
            {
              ++timestep;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              if (_verbose)
        	libMesh::out &lt;&lt; "Solving for state " &lt;&lt; timestep &lt;&lt; ", time " &lt;&lt; t &lt;&lt; "\n";
        
</pre>
</div>
<div class = "comment">
Move biharmonic one timestep forward
</div>

<div class ="fragment">
<pre>
              step();
        
</pre>
</div>
<div class = "comment">
Keep track of time and timestep
</div>

<div class ="fragment">
<pre>
              t += _dt;
        
</pre>
</div>
<div class = "comment">
Output
</div>

<div class ="fragment">
<pre>
              output(timestep,t,o_t);
            } // while(t &lt; _t1)
        
</pre>
</div>
<div class = "comment">
Force-write the final timestep
</div>

<div class ="fragment">
<pre>
          output(timestep,t,o_t,true);
        }
        
        
        
        
        
        Biharmonic::Biharmonic(UnstructuredMesh* m) :
            EquationSystems(*m),
            _mesh(m)
          {
</pre>
</div>
<div class = "comment">
Retrieve parameters and set defaults
</div>

<div class ="fragment">
<pre>
            _verbose      = false; if(on_command_line("--verbose")) _verbose = true;
            _growth       = false; if(on_command_line("--growth"))       _growth = true;
            _degenerate   = false; if(on_command_line("--degenerate"))   _degenerate = true;
            _cahn_hillard = false; if(on_command_line("--cahn_hillard")) _cahn_hillard = true;
            _netforce     = false; if(on_command_line("--netforce"))     _netforce = true;
            _kappa = command_line_value("kappa", 1.0);
        
</pre>
</div>
<div class = "comment">
"type of energy (double well, double obstacle, logarithmic+double well, logarithmic+double obstacle)"
</div>

<div class ="fragment">
<pre>
            std::string energy = command_line_value("energy", std::string("double_well"));
            if (energy == "double_well")
              _energy = DOUBLE_WELL;
            else if (energy == "double_obstacle")
              _energy = DOUBLE_OBSTACLE;
            else if (energy == "log_double_well")
              _energy = LOG_DOUBLE_WELL;
            else if (energy == "log_double_obstacle")
              _energy = LOG_DOUBLE_OBSTACLE;
            else
              ERROR(std::string("Unknown energy type: ") + energy);
        
            _tol     = command_line_value("tol",1.0e-8);
            _theta   = command_line_value("theta", .001);
            _theta_c = command_line_value("theta_c",1.0);
        
</pre>
</div>
<div class = "comment">
"order of log truncation (0=none, 2=quadratic, 3=cubic)"
</div>

<div class ="fragment">
<pre>
            _log_truncation = command_line_value("log_truncation", 2);
        
            if (!_log_truncation)
              libMesh::out &lt;&lt; "WARNING: no truncation is being used for the logarithmic free energy term.\nWARNING: division by zero possible!\n";
        
        
</pre>
</div>
<div class = "comment">
Dimension
</div>

<div class ="fragment">
<pre>
            _dim = command_line_value("dim",1);
        
            ASSERT((_dim &lt;= 3) && (_dim &gt; 0), "Invalid mesh dimension");
        
</pre>
</div>
<div class = "comment">
Build the mesh
Yes, it's better to make a coarse mesh and then refine it. We'll get to it later.
</div>

<div class ="fragment">
<pre>
            _N = command_line_value("N", 8);
            ASSERT(_N &gt; 0, "Invalid mesh size");
        
            switch (_dim)
              {
              case 1:
        	MeshTools::Generation::build_line(*_mesh, _N, 0.0, 1.0, EDGE2);
        	break;
              case 2:
        	MeshTools::Generation::build_square(*_mesh, _N, _N, 0.0, 1.0, 0.0, 1.0, QUAD4);
        	break;
              case 3:
        	MeshTools::Generation::build_cube(*_mesh, _N, _N, _N, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8);
        	break;
              default:
        	ASSERT((_dim &lt;= 3) && (_dim &gt; 0), "Invalid mesh dimension");
        	break;
              }
        
</pre>
</div>
<div class = "comment">
Determine the initial timestep size
</div>

<div class ="fragment">
<pre>
            _dt0 = command_line_value("dt", 1.0/(10*_kappa*_N*_N*_N*_N));
            ASSERT(_dt0&gt;=0, "Negative initial timestep");
        
            _t0 = command_line_value("min_time", 0.0);
            _t1 = command_line_value("max_time", _t0 + 50.0*_dt0);
            ASSERT(_t1 &gt;= _t0, "Final time less than initial time");
            _T = _t1 - _t0;
        
            _cnWeight = command_line_value("crank_nicholson_weight", 1.0);
            ASSERT(_cnWeight &lt;= 1 && _cnWeight &gt;= 0, "Crank-Nicholson weight must be between 0 and 1");
        
</pre>
</div>
<div class = "comment">
Initial state
</div>

<div class ="fragment">
<pre>
            _initialState = STRIP;
            std::string initialState = command_line_value("initial_state", std::string("strip"));
            if (initialState == std::string("ball"))
              _initialState = BALL;
            else if (initialState == std::string("strip"))
              _initialState = STRIP;
            else if (initialState == std::string("rod"))
              _initialState = ROD;
            else
              ERROR("Unknown initial state: neither ball nor rod nor srip");
        
            std::vector&lt;Real&gt; icenter;
            command_line_vector("initial_center", icenter);
        
</pre>
</div>
<div class = "comment">
Check that the point defining the center was in the right spatial dimension
</div>

<div class ="fragment">
<pre>
            if (icenter.size() &gt; _dim)
              ASSERT(icenter.size() &gt; _dim, "Invalid dimension for the initial state center of mass");
        
</pre>
</div>
<div class = "comment">
Pad
</div>

<div class ="fragment">
<pre>
            icenter.resize(3);
            for (unsigned int i = icenter.size(); i &lt; _dim; ++i)
              icenter[i] = 0.5;
        
            for (unsigned int i = _dim; i &lt; 3; ++i)
              icenter[i] = 0.0;
        
            _initialCenter = Point(icenter[0],icenter[1], icenter[2]);
            _initialWidth = command_line_value("initial_width", 0.125);
        
</pre>
</div>
<div class = "comment">
Build the main equation encapsulated in the JR (Jacobian-Residual or J(R) "jet of R") object
</div>

<div class ="fragment">
<pre>
            _jr = &(add_system&lt;Biharmonic::JR&gt;(std::string("Biharmonic::JR")));
        
</pre>
</div>
<div class = "comment">
Output options
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_EXODUS_API
            if (on_command_line("output_base"))
              _ofile_base = command_line_value("output_base", std::string("bih"));
        
            else
              {
        	switch(_dim)
        	  {
        	  case 1:
        	    _ofile_base = std::string("bih.1");
        	    break;
        	  case 2:
        	    _ofile_base = std::string("bih.2");
        	    break;
        	  case 3:
        	    _ofile_base = std::string("bih.3");
        	    break;
        	  default:
        	    _ofile_base = std::string("bih");
        	    break;
        	  }
              }
            _ofile = _ofile_base + ".e";
            _exio = new ExodusII_IO(*_mesh);
            _o_dt = command_line_value("output_dt", 0.0);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
          } // constructor
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file biharmonic_jr.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh.h"
        #include "libmesh/quadrature.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/fourth_error_estimators.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/periodic_boundaries.h"
        #include "libmesh/periodic_boundary.h"
        
</pre>
</div>
<div class = "comment">
Example includes
</div>

<div class ="fragment">
<pre>
        #include "biharmonic_jr.h"
        
        using namespace libMesh;
        
        Biharmonic::JR::JR(EquationSystems& eqSys,
        		   const std::string& name_in,
        		   const unsigned int number_in) :
          TransientNonlinearImplicitSystem(eqSys,name_in,number_in),
          _biharmonic(dynamic_cast&lt;Biharmonic&&gt;(eqSys))
        {
</pre>
</div>
<div class = "comment">
Check that we can actually compute second derivatives
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
          ERROR("Must have second derivatives enabled");
        #endif
        
        #ifdef LIBMESH_ENABLE_PERIODIC
</pre>
</div>
<div class = "comment">
Add periodicity to the mesh
</div>

<div class ="fragment">
<pre>
          DofMap& dof_map = get_dof_map();
          PeriodicBoundary xbdry(RealVectorValue(1.0, 0.0, 0.0));
        #if LIBMESH_DIM &gt; 1
          PeriodicBoundary ybdry(RealVectorValue(0.0, 1.0, 0.0));
        #endif
        #if LIBMESH_DIM &gt; 2
          PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, 1.0));
        #endif
        
          switch(_biharmonic._dim)
            {
            case 1:
              xbdry.myboundary = 0;
              xbdry.pairedboundary = 1;
              dof_map.add_periodic_boundary(xbdry);
              break;
        #if LIBMESH_DIM &gt; 1
            case 2:
              xbdry.myboundary = 3;
              xbdry.pairedboundary = 1;
              dof_map.add_periodic_boundary(xbdry);
              ybdry.myboundary = 0;
              ybdry.pairedboundary = 2;
              dof_map.add_periodic_boundary(ybdry);
              break;
        #endif
        #if LIBMESH_DIM &gt; 2
            case 3:
              xbdry.myboundary = 4;
              xbdry.pairedboundary = 2;
              dof_map.add_periodic_boundary(xbdry);
              ybdry.myboundary = 1;
              ybdry.pairedboundary = 3;
              dof_map.add_periodic_boundary(ybdry);
              zbdry.myboundary = 0;
              zbdry.pairedboundary = 5;
              dof_map.add_periodic_boundary(zbdry);
              break;
        #endif
            default:
              libmesh_error();
            }
        #endif // LIBMESH_ENABLE_PERIODIC
        
</pre>
</div>
<div class = "comment">
Adaptivity stuff is commented out for now...
#ifndef   LIBMESH_ENABLE_AMR
libmesh_example_assert(false, "--enable-amr");
#else
// In case we ever get around to doing mesh refinement.
_biharmonic._meshRefinement = new MeshRefinement(_mesh);

<br><br>// Tell the MeshRefinement object about the periodic boundaries
// so that it can get heuristics like level-one conformity and unrefined
// island elimination right.
_biharmonic._mesh_refinement->set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());
#endif // LIBMESH_ENABLE_AMR


<br><br>Adds the variable "u" to the system.
u will be approximated using Hermite elements
</div>

<div class ="fragment">
<pre>
          add_variable("u", THIRD, HERMITE);
        
</pre>
</div>
<div class = "comment">
Give the system an object to compute the initial state.
</div>

<div class ="fragment">
<pre>
          attach_init_object(*this);
        
</pre>
</div>
<div class = "comment">
Attache the R & J calculation object
</div>

<div class ="fragment">
<pre>
          nonlinear_solver-&gt;residual_and_jacobian_object = this;
        
</pre>
</div>
<div class = "comment">
Attach the bounds calculation object
</div>

<div class ="fragment">
<pre>
          nonlinear_solver-&gt;bounds_object = this;
        }
        
        
        
        
        
        void Biharmonic::JR::initialize()
        {
          if (_biharmonic._verbose)
            libMesh::out &lt;&lt; "&gt;&gt;&gt; Initializing Biharmonic::JR\n";
        
          Parameters parameters;
          parameters.set&lt;Point&gt;("center") = _biharmonic._initialCenter;
          parameters.set&lt;Real&gt;("width")   = _biharmonic._initialWidth;
        
          if (_biharmonic._initialState == Biharmonic::BALL)
            project_solution(Biharmonic::JR::InitialDensityBall, Biharmonic::JR::InitialGradientZero, parameters);
        
          if (_biharmonic._initialState == Biharmonic::ROD)
            project_solution(Biharmonic::JR::InitialDensityRod, Biharmonic::JR::InitialGradientZero, parameters);
        
          if (_biharmonic._initialState == Biharmonic::STRIP)
            project_solution(Biharmonic::JR::InitialDensityStrip, Biharmonic::JR::InitialGradientZero, parameters);
        
</pre>
</div>
<div class = "comment">
both states are equal
</div>

<div class ="fragment">
<pre>
          *(old_local_solution) = *(current_local_solution);
        
          if (_biharmonic._verbose)
            libMesh::out &lt;&lt; "&lt;&lt;&lt; Initializing Biharmonic::JR\n";
        }
        
        
        
        
        
        
        Number Biharmonic::JR::InitialDensityBall(const Point& p,
        					  const Parameters& parameters,
        					  const std::string&,
        					  const std::string&)
        {
</pre>
</div>
<div class = "comment">
Initialize with a ball in the middle, which is a segment in 1D, a disk in 2D and a ball in 3D.
</div>

<div class ="fragment">
<pre>
          Point center = parameters.get&lt;Point&gt;("center");
          Real width = parameters.get&lt;Real&gt;("width");
          Point pc = p-center;
          Real r = pc.size();
          return (r &lt; width) ? 1.0 : -0.5;
        }
        
        
        
        
        Number Biharmonic::JR::InitialDensityRod(const Point& p,
        					 const Parameters& parameters,
        					 const std::string&,
        					 const std::string&)
        {
</pre>
</div>
<div class = "comment">
Initialize with a rod in the middle so that we have a z-homogeneous system to model the 2D disk.
</div>

<div class ="fragment">
<pre>
          Point center = parameters.get&lt;Point&gt;("center");
          Real width = parameters.get&lt;Real&gt;("width");
          Real r = sqrt((p(0)-center(0))*(p(0)-center(0)) + (p(1)-center(1))*(p(1)-center(1)));
          return (r &lt; width) ? 1.0 : -0.5;
        }
        
        
        
        
        
        Number Biharmonic::JR::InitialDensityStrip(const Point& p,
        					   const Parameters& parameters,
        					   const std::string&,
        					   const std::string&)
        {
</pre>
</div>
<div class = "comment">
Initialize with a wide strip in the middle so that we have a yz-homogeneous system to model the 1D.
</div>

<div class ="fragment">
<pre>
          Point center = parameters.get&lt;Point&gt;("center");
          Real width = parameters.get&lt;Real&gt;("width");
          Real r = sqrt((p(0)-center(0))*(p(0)-center(0)));
          return (r &lt; width) ? 1.0 : -0.5;
        }
        
        
        
        
        Gradient Biharmonic::JR::InitialGradientZero(const Point&,
        					     const Parameters&,
        					     const std::string&,
        					     const std::string&)
        {
          return Gradient(0.0,0.0,0.0);
        }
        
        
        
        
        void Biharmonic::JR::residual_and_jacobian(const NumericVector&lt;Number&gt; &u,
        					   NumericVector&lt;Number&gt; *R,
        					   SparseMatrix&lt;Number&gt; *J,
        					   NonlinearImplicitSystem&)
        {
        #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          if (!R && !J)
            return;
        
</pre>
</div>
<div class = "comment">
Declare a performance log.  Give it a descriptive
string to identify what part of the code we are
logging, since there may be many PerfLogs in an
application.
</div>

<div class ="fragment">
<pre>
          PerfLog perf_log ("Biharmonic Residual and Jacobian", false);
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = get_dof_map();
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = dof_map.variable_type(0);
        
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  This can be thought
of as a pointer that will clean up after itself.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(_biharmonic._dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Quadrature rule for numerical integration.
With 2D triangles, the Clough quadrature rule puts a Gaussian
quadrature rule on each of the 3 subelements
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(_biharmonic._dim));
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule (qrule.get());
        
</pre>
</div>
<div class = "comment">
Here we define some references to element-specific data that
will be used to assemble the linear system.
We begin with the element Jacobian * quadrature weight at each
integration point.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The element shape functions' derivatives evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The element shape functions'  second derivatives evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealTensor&gt; &gt;& d2phi = fe-&gt;get_d2phi();
        
</pre>
</div>
<div class = "comment">
For efficiency we will compute shape function laplacians n times,
not n^2
</div>

<div class ="fragment">
<pre>
          std::vector&lt;Real&gt; Laplacian_phi_qp;
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Je" and "Re". More detail is in example 3.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Je;
          DenseVector&lt;Number&gt; Re;
        
</pre>
</div>
<div class = "comment">
This vector will hold the degree of freedom indices for
the element.  These define where in the global system
the element degrees of freedom get mapped.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Old solution
</div>

<div class ="fragment">
<pre>
          const NumericVector&lt;Number&gt;& u_old = *old_local_solution;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh.  We will
compute the element matrix and right-hand-side contribution.  See
example 3 for a discussion of the element iterators.


<br><br></div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = _biharmonic._mesh-&gt;active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = _biharmonic._mesh-&gt;active_local_elements_end();
        
          for ( ; el != end_el; ++el) {
</pre>
</div>
<div class = "comment">
Store a pointer to the element we are currently
working on.  This allows for nicer syntax later.
</div>

<div class ="fragment">
<pre>
            const Elem* elem = *el;
        
</pre>
</div>
<div class = "comment">
Get the degree of freedom indices for the
current element.  These define where in the global
matrix and right-hand-side this element will
contribute to.
</div>

<div class ="fragment">
<pre>
            dof_map.dof_indices (elem, dof_indices);
        
</pre>
</div>
<div class = "comment">
Compute the element-specific data for the current
element.  This involves computing the location of the
quadrature points (q_point) and the shape function
values/derivatives (phi, dphi,d2phi) for the current element.
</div>

<div class ="fragment">
<pre>
            fe-&gt;reinit (elem);
        
</pre>
</div>
<div class = "comment">
Zero the element matrix, the right-hand side and the Laplacian matrix
before summing them.
</div>

<div class ="fragment">
<pre>
            if (J)
              Je.resize(dof_indices.size(), dof_indices.size());
        
            if (R)
              Re.resize(dof_indices.size());
        
            Laplacian_phi_qp.resize(dof_indices.size());
        
            for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
              {
</pre>
</div>
<div class = "comment">
AUXILIARY QUANTITIES:
Residual and Jacobian share a few calculations:
at the very least, in the case of interfacial energy only with a constant mobility,
both calculations use Laplacian_phi_qp; more is shared the case of a concentration-dependent
mobility and bulk potentials.
</div>

<div class ="fragment">
<pre>
                Number u_qp = 0.0, u_old_qp = 0.0, Laplacian_u_qp = 0.0, Laplacian_u_old_qp = 0.0;
        	Gradient grad_u_qp(0.0,0.0,0.0), grad_u_old_qp(0.0,0.0,0.0);
        	Number M_qp = 1.0, M_old_qp = 1.0, M_prime_qp = 0.0, M_prime_old_qp = 0.0;
        
        	for (unsigned int i=0; i&lt;phi.size(); i++)
        	  {
        	    Laplacian_phi_qp[i] = d2phi[i][qp](0,0);
        	    grad_u_qp(0) += u(dof_indices[i])*dphi[i][qp](0);
        	    grad_u_old_qp(0) += u_old(dof_indices[i])*dphi[i][qp](0);
        
        	    if (_biharmonic._dim &gt; 1)
        	      {
        		Laplacian_phi_qp[i] += d2phi[i][qp](1,1);
        		grad_u_qp(1) += u(dof_indices[i])*dphi[i][qp](1);
        		grad_u_old_qp(1) += u_old(dof_indices[i])*dphi[i][qp](1);
        	      }
        	    if (_biharmonic._dim &gt; 2)
        	      {
        		Laplacian_phi_qp[i] += d2phi[i][qp](2,2);
        		grad_u_qp(2) += u(dof_indices[i])*dphi[i][qp](2);
        		grad_u_old_qp(2) += u_old(dof_indices[i])*dphi[i][qp](2);
        	      }
        	    u_qp     += phi[i][qp]*u(dof_indices[i]);
        	    u_old_qp += phi[i][qp]*u_old(dof_indices[i]);
        	    Laplacian_u_qp     += Laplacian_phi_qp[i]*u(dof_indices[i]);
        	    Laplacian_u_old_qp += Laplacian_phi_qp[i]*u_old(dof_indices[i]);
        	  } // for i
        
        	if (_biharmonic._degenerate)
        	  {
        	    M_qp           = 1.0 - u_qp*u_qp;
        	    M_old_qp       = 1.0 - u_old_qp*u_old_qp;
        	    M_prime_qp     = -2.0*u_qp;
        	    M_prime_old_qp = -2.0*u_old_qp;
        	  }
        
</pre>
</div>
<div class = "comment">
ELEMENT RESIDUAL AND JACOBIAN
</div>

<div class ="fragment">
<pre>
                for (unsigned int i=0; i&lt;phi.size(); i++)
        	  {
</pre>
</div>
<div class = "comment">
RESIDUAL
</div>

<div class ="fragment">
<pre>
                    if (R)
        	      {
        		Number ri = 0.0, ri_old = 0.0;
        		ri     -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_u_qp;
        		ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._kappa*Laplacian_u_old_qp;
        
        		if (_biharmonic._degenerate)
        		  {
        		    ri       -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp);
        		    ri_old   -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*(_biharmonic._kappa*Laplacian_u_old_qp);
        		  }
        
        		if (_biharmonic._cahn_hillard)
        		  {
        		    if (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
        		      {
        			ri += Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
        			ri_old += Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
        			if (_biharmonic._degenerate)
        			  {
        			    ri     += (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
        			    ri_old += (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
        			  }
        		      }// if(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
        
        		    if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        		      {
        			ri -= Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*u_qp;
        			ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*u_old_qp;
        			if (_biharmonic._degenerate)
        			  {
        			    ri     -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*u_qp;
        			    ri_old -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*u_old_qp;
        			  }
        		      } // if(_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        
        		    if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        		      {
        			switch(_biharmonic._log_truncation)
        			  {
        			  case 2:
        			    break;
        			  case 3:
        			    break;
        			  default:
        			    break;
        			  }// switch(_biharmonic._log_truncation)
        		      }// if(_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        		  }// if(_biharmonic._cahn_hillard)
        		Re(i) += JxW[qp]*((u_qp-u_old_qp)*phi[i][qp]-_biharmonic._dt*0.5*((2.0-_biharmonic._cnWeight)*ri + _biharmonic._cnWeight*ri_old));
        	      } // if (R)
        
</pre>
</div>
<div class = "comment">
JACOBIAN
</div>

<div class ="fragment">
<pre>
                    if (J)
        	      {
        		Number M_prime_prime_qp = 0.0;
        		if(_biharmonic._degenerate) M_prime_prime_qp = -2.0;
        		for (unsigned int j=0; j&lt;phi.size(); j++)
        		  {
        		    Number ri_j = 0.0;
        		    ri_j -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_phi_qp[j];
        		    if (_biharmonic._degenerate)
        		      {
        			ri_j -=
        			  Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._kappa*Laplacian_u_qp               +
        			  (dphi[i][qp]*dphi[j][qp])*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp)                  +
        			  (dphi[i][qp]*grad_u_qp)*(M_prime_prime_qp*phi[j][qp])*(_biharmonic._kappa*Laplacian_u_qp) +
        			  (dphi[i][qp]*grad_u_qp)*(M_prime_qp)*(_biharmonic._kappa*Laplacian_phi_qp[j]);
        		      }
        
        		    if (_biharmonic._cahn_hillard)
        		      {
        			if(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
        			  {
        			    ri_j +=
        			      Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp +
        			      Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp]        +
        			      (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp      +
        			      (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp  +
        			      (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp];
        			  }// if(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
        
        			if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        			  {
        			    ri_j -=
        			      Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*u_qp                   +
        			      Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*phi[j][qp]                              +
        			      (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*u_qp                        +
        			      (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*u_qp                    +
        			      (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*phi[j][qp];
        			  } // if(_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        
        			if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        			  {
        			    switch(_biharmonic._log_truncation)
        			      {
        			      case 2:
        				break;
        			      case 3:
        				break;
        			      default:
        				break;
        			      }// switch(_biharmonic._log_truncation)
        			  }// if(_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
        		      }// if(_biharmonic._cahn_hillard)
        		    Je(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] - 0.5*_biharmonic._dt*(2.0-_biharmonic._cnWeight)*ri_j);
        		  } // for j
        	      } // if (J)
        	  } // for i
              } // for qp
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p SparseMatrix::add_matrix()
and \p NumericVector::add_vector() members do this for us.
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
            if (R)
              {
</pre>
</div>
<div class = "comment">
If the mesh has hanging nodes (e.g., as a result of refinement), those need to be constrained.
</div>

<div class ="fragment">
<pre>
                dof_map.constrain_element_vector(Re, dof_indices);
        	R-&gt;add_vector(Re, dof_indices);
              }
        
            if (J)
              {
</pre>
</div>
<div class = "comment">
If the mesh has hanging nodes (e.g., as a result of refinement), those need to be constrained.
</div>

<div class ="fragment">
<pre>
                dof_map.constrain_element_matrix(Je, dof_indices);
        	J-&gt;add_matrix(Je, dof_indices);
              }
          } // for el
        #endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
        }
        
        
        
        
        
        void Biharmonic::JR::bounds(NumericVector&lt;Number&gt; &XL, NumericVector&lt;Number&gt;& XU, NonlinearImplicitSystem&)
        {
</pre>
</div>
<div class = "comment">
sys is actually ignored, since it should be the same as *this.


<br><br>Declare a performance log.  Give it a descriptive
string to identify what part of the code we are
logging, since there may be many PerfLogs in an
application.
</div>

<div class ="fragment">
<pre>
          PerfLog perf_log ("Biharmonic bounds", false);
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = get_dof_map();
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = dof_map.variable_type(0);
        
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  This can be thought
of as a pointer that will clean up after itself.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(_biharmonic._dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the bound vectors contributions.
</div>

<div class ="fragment">
<pre>
          DenseVector&lt;Number&gt; XLe, XUe;
        
</pre>
</div>
<div class = "comment">
These vector will hold the degree of freedom indices for
the element.  These define where in the global system
the element degrees of freedom get mapped.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;dof_id_type&gt; dof_indices;
        
          MeshBase::const_element_iterator       el     = _biharmonic._mesh-&gt;active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = _biharmonic._mesh-&gt;active_local_elements_end();
        
          for ( ; el != end_el; ++el)
            {
</pre>
</div>
<div class = "comment">
Extract the shape function to be evaluated at the nodes
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
Get the degree of freedom indices for the current element.
They are in 1-1 correspondence with shape functions phi
and define where in the global vector this element will.
</div>

<div class ="fragment">
<pre>
              dof_map.dof_indices (*el, dof_indices);
        
</pre>
</div>
<div class = "comment">
Resize the local bounds vectors (zeroing them out in the process).
</div>

<div class ="fragment">
<pre>
              XLe.resize(dof_indices.size());
              XUe.resize(dof_indices.size());
        
</pre>
</div>
<div class = "comment">
Extract the element node coordinates in the reference frame
</div>

<div class ="fragment">
<pre>
              std::vector&lt;Point&gt; nodes;
              fe-&gt;get_refspace_nodes((*el)-&gt;type(), nodes);
        
</pre>
</div>
<div class = "comment">
Evaluate the shape functions at the nodes
</div>

<div class ="fragment">
<pre>
              fe-&gt;reinit(*el, &nodes);
        
</pre>
</div>
<div class = "comment">
Construct the bounds based on the value of the i-th phi at the nodes.
Observe that this doesn't really work in general: we rely on the fact
that for Hermite elements each shape function is nonzero at most at a
single node.
More generally the bounds must be constructed by inspecting a "mass-like"
matrix (m_{ij}) of the shape functions (i) evaluated at their corresponding nodes (j).
The constraints imposed on the dofs (d_i) are then are -1 \leq \sum_i d_i m_{ij} \leq 1,
since \sum_i d_i m_{ij} is the value of the solution at the j-th node.
Auxiliary variables will need to be introduced to reduce this to a "box" constraint.
Additional complications will arise since m might be singular (as is the case for Hermite,
which, however, is easily handled by inspection).
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i&lt;phi.size(); ++i)
        	{
</pre>
</div>
<div class = "comment">
FIXME: should be able to define INF and pass it to the solve
</div>

<div class ="fragment">
<pre>
                  Real infinity = 1.0e20;
        	  Real bound = infinity;
        	  for(unsigned int j = 0; j &lt; nodes.size(); ++j) {
        	    if(phi[i][j]) {
        	      bound = 1.0/fabs(phi[i][j]);
        	      break;
        	    }
        	  }
        
</pre>
</div>
<div class = "comment">
The value of the solution at this node must be between 1.0 and -1.0.
Based on the value of phi(i)(i) the nodal coordinate must be between 1.0/phi(i)(i) and its negative.
</div>

<div class ="fragment">
<pre>
                  XLe(i) = -bound;
        	  XUe(i) = bound;
        	}
</pre>
</div>
<div class = "comment">
The element bound vectors are now built for this element.
Insert them into the global vectors, potentially overwriting
the same dof contributions from other elements: no matter --
the bounds are always -1.0 and 1.0.
</div>

<div class ="fragment">
<pre>
              XL.insert(XLe, dof_indices);
              XU.insert(XUe, dof_indices);
            }
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex7.C with comments: </h1> 
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
              Biharmonic::Create(&biharmonic,init.comm());
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
<br><br><br> <h1> The source file biharmonic.h without comments: </h1> 
<pre> 
  #ifndef __biharmonic_h__
  #define __biharmonic_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  
  using libMesh::EquationSystems;
  using libMesh::ExodusII_IO;
  using libMesh::MeshRefinement;
  using libMesh::Point;
  using libMesh::Real;
  using libMesh::UnstructuredMesh;
  
  #define ERROR(message)                                                                                         \
    <B><FONT COLOR="#A020F0">do</FONT></B> {                                                                                                         \
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::err &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: &quot;</FONT></B> &lt;&lt; message &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;                                                              \
      libmesh_error();							                                       \
    } <B><FONT COLOR="#A020F0">while</FONT></B>(0)
  
  #define ASSERT(asserted, message) \
    <B><FONT COLOR="#A020F0">do</FONT></B> {                    \
      <B><FONT COLOR="#A020F0">if</FONT></B>(!(asserted)) {     \
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::err &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Assertion '&quot;</FONT></B> #asserted <B><FONT COLOR="#BC8F8F">&quot;' violated: &quot;</FONT></B> #message; \
        libmesh_error();    \
      }                     \
    } <B><FONT COLOR="#A020F0">while</FONT></B>(0)
  
  
  <I><FONT COLOR="#B22222">/**
   * The Biharmonic class encapsulates most of the data structures
   * necessary to calculate the biharmonic residual and Jacobian,
   * auxiliary quantities, to take a timestep, and to output the state --
   * biharmonic solution and vectors of auxiliary quantities.
   *
   * The main reason for this design is to have a data structure that
   * has all of the necessary data in one place, where all of the
   * calculation subroutines can access these data. Currently these data
   * are split up among several interdependent objects with no clear
   * hierarchy between them: mesh, equation system, equation system
   * bundle, residual/Jacobian calculator.
   *
   * Since no object contains all others and the data are distributed
   * among many objects, the natural control and data flow resides outside
   * of these objects and is typically implemented in main().  We,
   * however, would like to split the calculation into natural chunks --
   * subroutines -- while retaining these subroutines access to the common
   * necessary data -- biharmonic parameters, mesh and time interval
   * sizes, etc. Thus, class Biharmonic.  Finally, making Biharmonic
   * inherit from EquationSystems makes it possible to include it in the
   * most common callbacks that do not pass back a user context, but only
   * an EquationSystems object.
   */</FONT></I>
  <B><FONT COLOR="#228B22">class</FONT></B> Biharmonic : <B><FONT COLOR="#228B22">public</FONT></B> EquationSystems
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    <B><FONT COLOR="#228B22">enum</FONT></B> InitialStateEnum {STRIP = 0,
  			 ROD   = 1,
  			 BALL  = 2};
  
    <B><FONT COLOR="#228B22">enum</FONT></B> FreeEnergyEnum {DOUBLE_WELL         = 1,
  		       DOUBLE_OBSTACLE     = 2,
  		       LOG_DOUBLE_WELL     = 3,
  		       LOG_DOUBLE_OBSTACLE = 4};
  
    <I><FONT COLOR="#B22222">/**
     * Static creation/destruction routines.  FIXME - this looks like
     * object-oriented C, can we get rid of it?
     */</FONT></I>
    <B><FONT COLOR="#228B22">static</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> Create(Biharmonic** b, <B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator &amp;comm);
    <B><FONT COLOR="#228B22">static</FONT></B> <B><FONT COLOR="#228B22">void</FONT></B> Destroy(Biharmonic** b);
  
  
    <I><FONT COLOR="#B22222">/**
     * Constructor retrieves command-line options, setting  defaults, if necessary.
     * It then builds the mesh using these options, then the equations systems around it,
     * and, finally, sets up the output.
     * We recommend that this be used through the factory Create function, which allocates
     * the mesh. In that case don't forget to call Destroy at the end, to free the mesh up.
     */</FONT></I>
    Biharmonic(UnstructuredMesh* m);
  
  
    <I><FONT COLOR="#B22222">/**
     * Destructor
     */</FONT></I>
    ~Biharmonic()
    {
    }
  
  
    <B><FONT COLOR="#228B22">bool</FONT></B> verbose()         { <B><FONT COLOR="#A020F0">return</FONT></B> _verbose; }
    Real dt0()             { <B><FONT COLOR="#A020F0">return</FONT></B> _dt0; }
    Real dt()              { <B><FONT COLOR="#A020F0">return</FONT></B> _dt; }
  
  
    <B><FONT COLOR="#228B22">void</FONT></B> viewParameters();
    <B><FONT COLOR="#228B22">void</FONT></B> init();
    <B><FONT COLOR="#228B22">void</FONT></B> step(<B><FONT COLOR="#228B22">const</FONT></B> Real&amp; dt = -1.0);
    <B><FONT COLOR="#228B22">void</FONT></B> output(<B><FONT COLOR="#228B22">int</FONT></B> timestep, <B><FONT COLOR="#228B22">const</FONT></B> Real&amp; t, Real&amp; o_t, <B><FONT COLOR="#228B22">bool</FONT></B> force = false);
    <B><FONT COLOR="#228B22">void</FONT></B> run();
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>  _dim, _N;
    Real _kappa, _theta, _theta_c;
    Real _tol;
    <B><FONT COLOR="#228B22">bool</FONT></B> _growth, _degenerate, _cahn_hillard, _netforce;
    FreeEnergyEnum  _energy;
    <B><FONT COLOR="#228B22">int</FONT></B> _log_truncation;
    <B><FONT COLOR="#228B22">bool</FONT></B> _verbose;
    InitialStateEnum  _initialState;
    Point _initialCenter;
    Real _initialWidth;
    Real _dt0, _dt, _t0, _T, _t1;
    Real _cnWeight;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string  _ofile_base, _ofile;
    ExodusII_IO* _exio;
    Real    _o_dt;
    <B><FONT COLOR="#228B22">int</FONT></B>     _o_count;
    <B><FONT COLOR="#228B22">friend</FONT></B> <B><FONT COLOR="#228B22">class</FONT></B> JR;
    <B><FONT COLOR="#228B22">class</FONT></B> JR;       <I><FONT COLOR="#B22222">// forward
</FONT></I>    UnstructuredMesh*                       _mesh;
    MeshRefinement*                         _meshRefinement;
    JR*                                     _jr;
  };
  
  
  
  
  
  
  
  #endif <I><FONT COLOR="#B22222">// __biharmonic_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file biharmonic_jr.h without comments: </h1> 
<pre> 
  #ifndef __biharmonic_jr_h__
  #define __biharmonic_jr_h__
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/nonlinear_solver.h&quot;</FONT></B>
  
  
  #include <B><FONT COLOR="#BC8F8F">&quot;biharmonic.h&quot;</FONT></B>
  
  using libMesh::EquationSystems;
  using libMesh::Gradient;
  using libMesh::NonlinearImplicitSystem;
  using libMesh::Number;
  using libMesh::NumericVector;
  using libMesh::Parameters;
  using libMesh::Point;
  using libMesh::SparseMatrix;
  using libMesh::System;
  using libMesh::TransientNonlinearImplicitSystem;
  
  
  <I><FONT COLOR="#B22222">/**
   * Biharmonic's friend class definition
   */</FONT></I>
  <B><FONT COLOR="#228B22">class</FONT></B> Biharmonic::JR : <B><FONT COLOR="#228B22">public</FONT></B> TransientNonlinearImplicitSystem,
  		       <B><FONT COLOR="#228B22">public</FONT></B> NonlinearImplicitSystem::ComputeResidualandJacobian,
  		       <B><FONT COLOR="#228B22">public</FONT></B> NonlinearImplicitSystem::ComputeBounds,
  		       <B><FONT COLOR="#228B22">public</FONT></B> System::Initialization
  {
  <B><FONT COLOR="#228B22">public</FONT></B>:
    <I><FONT COLOR="#B22222">/**
     * Constructor.
     */</FONT></I>
    JR(EquationSystems&amp; eqSys, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name, <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number);
  
    <B><FONT COLOR="#228B22">void</FONT></B> initialize();
  
    <I><FONT COLOR="#B22222">/**
     * Static functions to be used for initialization
     */</FONT></I>
    <B><FONT COLOR="#228B22">static</FONT></B> Number InitialDensityBall(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p, <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
    <B><FONT COLOR="#228B22">static</FONT></B> Number InitialDensityRod(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p, <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
    <B><FONT COLOR="#228B22">static</FONT></B> Number InitialDensityStrip(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p, <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
    <B><FONT COLOR="#228B22">static</FONT></B> Gradient InitialGradientZero(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;, <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
    <I><FONT COLOR="#B22222">/**
     * The residual and Jacobian assembly function for the Biharmonic system.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> residual_and_jacobian(<B><FONT COLOR="#228B22">const</FONT></B> NumericVector&lt;Number&gt;&amp; u,
  			     NumericVector&lt;Number&gt;* R,
  			     SparseMatrix&lt;Number&gt;* J,
  			     NonlinearImplicitSystem&amp;);
  
  
    <I><FONT COLOR="#B22222">/**
     * Function defining the bounds of the Biharmonic system.
     */</FONT></I>
    <B><FONT COLOR="#228B22">void</FONT></B> bounds(NumericVector&lt;Number&gt;&amp; XL,
  	      NumericVector&lt;Number&gt;&amp; XU,
  	      NonlinearImplicitSystem&amp;);
  
  <B><FONT COLOR="#228B22">private</FONT></B>:
    Biharmonic&amp; _biharmonic;
  };
  
  
  #endif <I><FONT COLOR="#B22222">// __biharmonic_jr_h__
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file biharmonic.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;biharmonic.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;biharmonic_jr.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::Create(Biharmonic** b, <B><FONT COLOR="#228B22">const</FONT></B> Parallel::Communicator &amp;comm)
  {
    SerialMesh* mesh = <B><FONT COLOR="#A020F0">new</FONT></B> SerialMesh(comm);
    Biharmonic *biharmonic = <B><FONT COLOR="#A020F0">new</FONT></B> Biharmonic(mesh);
    *b = biharmonic;
  }
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::Destroy(Biharmonic** b)
  {
    Biharmonic* biharmonic = *b;
    UnstructuredMesh* mesh = biharmonic-&gt;_mesh;
    <B><FONT COLOR="#A020F0">delete</FONT></B> biharmonic;
    <B><FONT COLOR="#A020F0">delete</FONT></B> mesh;
    *b = NULL;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::viewParameters()
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Biharmonic parameters:\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;verbose mode is on\n&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B>
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;verbose mode is off\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;mesh dimension           = &quot;</FONT></B> &lt;&lt; _dim &lt;&lt;               <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial linear mesh size = &quot;</FONT></B> &lt;&lt; _N   &lt;&lt;               <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;kappa                    = &quot;</FONT></B> &lt;&lt; _kappa &lt;&lt;             <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;growth                   = &quot;</FONT></B> &lt;&lt; (<B><FONT COLOR="#228B22">int</FONT></B>)_growth &lt;&lt;       <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;degenerate               = &quot;</FONT></B> &lt;&lt; (<B><FONT COLOR="#228B22">int</FONT></B>)_degenerate &lt;&lt;   <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Cahn-Hillard             = &quot;</FONT></B> &lt;&lt; (<B><FONT COLOR="#228B22">int</FONT></B>)_cahn_hillard &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;netforce                 = &quot;</FONT></B> &lt;&lt; (<B><FONT COLOR="#228B22">int</FONT></B>)_netforce &lt;&lt;     <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;energy                   = &quot;</FONT></B> &lt;&lt; _energy        &lt;&lt;     <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;tol                      = &quot;</FONT></B> &lt;&lt; _tol           &lt;&lt;     <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;theta                    = &quot;</FONT></B> &lt;&lt; _theta         &lt;&lt;     <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;theta_c                  = &quot;</FONT></B> &lt;&lt; _theta_c       &lt;&lt;     <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;log truncation           = &quot;</FONT></B> &lt;&lt; _log_truncation &lt;&lt;    <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial timestep size    = &quot;</FONT></B> &lt;&lt; _dt0            &lt;&lt;    <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_initialState == STRIP)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial state:             strip\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_initialState == ROD)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial state:             rod\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_initialState == BALL)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial state:             ball\n&quot;</FONT></B>;
  
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial state center     = &quot;</FONT></B> &lt;&lt; _initialCenter(0) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial state width      = &quot;</FONT></B> &lt;&lt; _initialWidth &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;initial time (min_time)  = &quot;</FONT></B> &lt;&lt; _t0 &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;integration time         = &quot;</FONT></B> &lt;&lt; _T  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;final time   (max_time)  = &quot;</FONT></B> &lt;&lt; _t1 &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Crank-Nicholson weight   = &quot;</FONT></B> &lt;&lt; _cnWeight &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Output timestep          = &quot;</FONT></B> &lt;&lt; _o_dt &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Output filename base:      &quot;</FONT></B> &lt;&lt;  _ofile_base &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::init()
  {
    <B><FONT COLOR="#A020F0">if</FONT></B>(_verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&gt;&gt;&gt; Initializing Biharmonic\n&quot;</FONT></B>;
  
    _dt  =  0;
    _o_count = 0;
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;EquationSystems::init();
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(_verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&lt;&lt;&lt; Initializing Biharmonic\n&quot;</FONT></B>;
  }
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::step(<B><FONT COLOR="#228B22">const</FONT></B> Real&amp; dt_in)
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (dt_in &lt; 0)
      _dt = _dt0;
    <B><FONT COLOR="#A020F0">else</FONT></B>
      _dt = dt_in;
  
    *(_jr-&gt;old_local_solution) = *(_jr-&gt;current_local_solution);
  
    _jr-&gt;solve();
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::output(<B><FONT COLOR="#228B22">int</FONT></B> timestep, <B><FONT COLOR="#228B22">const</FONT></B> Real&amp; t, Real&amp; o_t, <B><FONT COLOR="#228B22">bool</FONT></B> force)
  {
  #ifdef LIBMESH_HAVE_EXODUS_API
    <B><FONT COLOR="#A020F0">if</FONT></B> (!force &amp;&amp; t - o_t &lt; _o_dt)
      <B><FONT COLOR="#A020F0">return</FONT></B>;
  
    ++_o_count;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Writing state &quot;</FONT></B> &lt;&lt; timestep &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; at time &quot;</FONT></B> &lt;&lt; t &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; to file &quot;</FONT></B> &lt;&lt; _ofile &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;; output a total of &quot;</FONT></B> &lt;&lt; _o_count &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; states so far\n&quot;</FONT></B>;
  
    _exio-&gt;write_timestep(_ofile, *<B><FONT COLOR="#A020F0">this</FONT></B>, timestep, t);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (!force)
      o_t = t;
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::run()
  {
    Real t = _t0, o_t = 0.0;
    <B><FONT COLOR="#228B22">int</FONT></B> timestep = 1;
  
    output(timestep,t,o_t,true);
  
    <B><FONT COLOR="#A020F0">while</FONT></B> (t &lt; _t1)
      {
        ++timestep;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (_verbose)
  	<B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solving for state &quot;</FONT></B> &lt;&lt; timestep &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time &quot;</FONT></B> &lt;&lt; t &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>;
  
        step();
  
        t += _dt;
  
        output(timestep,t,o_t);
      } <I><FONT COLOR="#B22222">// while(t &lt; _t1)
</FONT></I>  
    output(timestep,t,o_t,true);
  }
  
  
  
  
  
  <B><FONT COLOR="#5F9EA0">Biharmonic</FONT></B>::Biharmonic(UnstructuredMesh* m) :
      EquationSystems(*m),
      _mesh(m)
    {
      _verbose      = false; <B><FONT COLOR="#A020F0">if</FONT></B>(on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--verbose&quot;</FONT></B>)) _verbose = true;
      _growth       = false; <B><FONT COLOR="#A020F0">if</FONT></B>(on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--growth&quot;</FONT></B>))       _growth = true;
      _degenerate   = false; <B><FONT COLOR="#A020F0">if</FONT></B>(on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--degenerate&quot;</FONT></B>))   _degenerate = true;
      _cahn_hillard = false; <B><FONT COLOR="#A020F0">if</FONT></B>(on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--cahn_hillard&quot;</FONT></B>)) _cahn_hillard = true;
      _netforce     = false; <B><FONT COLOR="#A020F0">if</FONT></B>(on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--netforce&quot;</FONT></B>))     _netforce = true;
      _kappa = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;kappa&quot;</FONT></B>, 1.0);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string energy = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;energy&quot;</FONT></B>, std::string(<B><FONT COLOR="#BC8F8F">&quot;double_well&quot;</FONT></B>));
      <B><FONT COLOR="#A020F0">if</FONT></B> (energy == <B><FONT COLOR="#BC8F8F">&quot;double_well&quot;</FONT></B>)
        _energy = DOUBLE_WELL;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (energy == <B><FONT COLOR="#BC8F8F">&quot;double_obstacle&quot;</FONT></B>)
        _energy = DOUBLE_OBSTACLE;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (energy == <B><FONT COLOR="#BC8F8F">&quot;log_double_well&quot;</FONT></B>)
        _energy = LOG_DOUBLE_WELL;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (energy == <B><FONT COLOR="#BC8F8F">&quot;log_double_obstacle&quot;</FONT></B>)
        _energy = LOG_DOUBLE_OBSTACLE;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        ERROR(std::string(<B><FONT COLOR="#BC8F8F">&quot;Unknown energy type: &quot;</FONT></B>) + energy);
  
      _tol     = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;tol&quot;</FONT></B>,1.0e-8);
      _theta   = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;theta&quot;</FONT></B>, .001);
      _theta_c = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;theta_c&quot;</FONT></B>,1.0);
  
      _log_truncation = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;log_truncation&quot;</FONT></B>, 2);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (!_log_truncation)
        <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;WARNING: no truncation is being used for the logarithmic free energy term.\nWARNING: division by zero possible!\n&quot;</FONT></B>;
  
  
      _dim = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;dim&quot;</FONT></B>,1);
  
      ASSERT((_dim &lt;= 3) &amp;&amp; (_dim &gt; 0), <B><FONT COLOR="#BC8F8F">&quot;Invalid mesh dimension&quot;</FONT></B>);
  
      _N = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;N&quot;</FONT></B>, 8);
      ASSERT(_N &gt; 0, <B><FONT COLOR="#BC8F8F">&quot;Invalid mesh size&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">switch</FONT></B> (_dim)
        {
        <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
  	<B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_line(*_mesh, _N, 0.0, 1.0, EDGE2);
  	<B><FONT COLOR="#A020F0">break</FONT></B>;
        <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
  	<B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square(*_mesh, _N, _N, 0.0, 1.0, 0.0, 1.0, QUAD4);
  	<B><FONT COLOR="#A020F0">break</FONT></B>;
        <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
  	<B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube(*_mesh, _N, _N, _N, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8);
  	<B><FONT COLOR="#A020F0">break</FONT></B>;
        <B><FONT COLOR="#5F9EA0">default</FONT></B>:
  	ASSERT((_dim &lt;= 3) &amp;&amp; (_dim &gt; 0), <B><FONT COLOR="#BC8F8F">&quot;Invalid mesh dimension&quot;</FONT></B>);
  	<B><FONT COLOR="#A020F0">break</FONT></B>;
        }
  
      _dt0 = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>, 1.0/(10*_kappa*_N*_N*_N*_N));
      ASSERT(_dt0&gt;=0, <B><FONT COLOR="#BC8F8F">&quot;Negative initial timestep&quot;</FONT></B>);
  
      _t0 = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;min_time&quot;</FONT></B>, 0.0);
      _t1 = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;max_time&quot;</FONT></B>, _t0 + 50.0*_dt0);
      ASSERT(_t1 &gt;= _t0, <B><FONT COLOR="#BC8F8F">&quot;Final time less than initial time&quot;</FONT></B>);
      _T = _t1 - _t0;
  
      _cnWeight = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;crank_nicholson_weight&quot;</FONT></B>, 1.0);
      ASSERT(_cnWeight &lt;= 1 &amp;&amp; _cnWeight &gt;= 0, <B><FONT COLOR="#BC8F8F">&quot;Crank-Nicholson weight must be between 0 and 1&quot;</FONT></B>);
  
      _initialState = STRIP;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string initialState = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;initial_state&quot;</FONT></B>, std::string(<B><FONT COLOR="#BC8F8F">&quot;strip&quot;</FONT></B>));
      <B><FONT COLOR="#A020F0">if</FONT></B> (initialState == std::string(<B><FONT COLOR="#BC8F8F">&quot;ball&quot;</FONT></B>))
        _initialState = BALL;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (initialState == std::string(<B><FONT COLOR="#BC8F8F">&quot;strip&quot;</FONT></B>))
        _initialState = STRIP;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (initialState == std::string(<B><FONT COLOR="#BC8F8F">&quot;rod&quot;</FONT></B>))
        _initialState = ROD;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        ERROR(<B><FONT COLOR="#BC8F8F">&quot;Unknown initial state: neither ball nor rod nor srip&quot;</FONT></B>);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; icenter;
      command_line_vector(<B><FONT COLOR="#BC8F8F">&quot;initial_center&quot;</FONT></B>, icenter);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (icenter.size() &gt; _dim)
        ASSERT(icenter.size() &gt; _dim, <B><FONT COLOR="#BC8F8F">&quot;Invalid dimension for the initial state center of mass&quot;</FONT></B>);
  
      icenter.resize(3);
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = icenter.size(); i &lt; _dim; ++i)
        icenter[i] = 0.5;
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i = _dim; i &lt; 3; ++i)
        icenter[i] = 0.0;
  
      _initialCenter = Point(icenter[0],icenter[1], icenter[2]);
      _initialWidth = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;initial_width&quot;</FONT></B>, 0.125);
  
      _jr = &amp;(add_system&lt;Biharmonic::JR&gt;(std::string(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic::JR&quot;</FONT></B>)));
  
  #ifdef LIBMESH_HAVE_EXODUS_API
      <B><FONT COLOR="#A020F0">if</FONT></B> (on_command_line(<B><FONT COLOR="#BC8F8F">&quot;output_base&quot;</FONT></B>))
        _ofile_base = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;output_base&quot;</FONT></B>, std::string(<B><FONT COLOR="#BC8F8F">&quot;bih&quot;</FONT></B>));
  
      <B><FONT COLOR="#A020F0">else</FONT></B>
        {
  	<B><FONT COLOR="#A020F0">switch</FONT></B>(_dim)
  	  {
  	  <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
  	    _ofile_base = std::string(<B><FONT COLOR="#BC8F8F">&quot;bih.1&quot;</FONT></B>);
  	    <B><FONT COLOR="#A020F0">break</FONT></B>;
  	  <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
  	    _ofile_base = std::string(<B><FONT COLOR="#BC8F8F">&quot;bih.2&quot;</FONT></B>);
  	    <B><FONT COLOR="#A020F0">break</FONT></B>;
  	  <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
  	    _ofile_base = std::string(<B><FONT COLOR="#BC8F8F">&quot;bih.3&quot;</FONT></B>);
  	    <B><FONT COLOR="#A020F0">break</FONT></B>;
  	  <B><FONT COLOR="#5F9EA0">default</FONT></B>:
  	    _ofile_base = std::string(<B><FONT COLOR="#BC8F8F">&quot;bih&quot;</FONT></B>);
  	    <B><FONT COLOR="#A020F0">break</FONT></B>;
  	  }
        }
      _ofile = _ofile_base + <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
      _exio = <B><FONT COLOR="#A020F0">new</FONT></B> ExodusII_IO(*_mesh);
      _o_dt = command_line_value(<B><FONT COLOR="#BC8F8F">&quot;output_dt&quot;</FONT></B>, 0.0);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>    } <I><FONT COLOR="#B22222">// constructor
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file biharmonic_jr.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fourth_error_estimators.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundary.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;biharmonic_jr.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#5F9EA0">Biharmonic</FONT></B>::JR::JR(EquationSystems&amp; eqSys,
  		   <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; name_in,
  		   <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> number_in) :
    TransientNonlinearImplicitSystem(eqSys,name_in,number_in),
    _biharmonic(dynamic_cast&lt;Biharmonic&amp;&gt;(eqSys))
  {
  #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
    ERROR(<B><FONT COLOR="#BC8F8F">&quot;Must have second derivatives enabled&quot;</FONT></B>);
  #endif
  
  #ifdef LIBMESH_ENABLE_PERIODIC
    DofMap&amp; dof_map = get_dof_map();
    PeriodicBoundary xbdry(RealVectorValue(1.0, 0.0, 0.0));
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 1
    PeriodicBoundary ybdry(RealVectorValue(0.0, 1.0, 0.0));
  #endif
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 2
    PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, 1.0));
  #endif
  
    <B><FONT COLOR="#A020F0">switch</FONT></B>(_biharmonic._dim)
      {
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">1</FONT></B>:
        xbdry.myboundary = 0;
        xbdry.pairedboundary = 1;
        dof_map.add_periodic_boundary(xbdry);
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 1
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
        xbdry.myboundary = 3;
        xbdry.pairedboundary = 1;
        dof_map.add_periodic_boundary(xbdry);
        ybdry.myboundary = 0;
        ybdry.pairedboundary = 2;
        dof_map.add_periodic_boundary(ybdry);
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  #endif
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 2
      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
        xbdry.myboundary = 4;
        xbdry.pairedboundary = 2;
        dof_map.add_periodic_boundary(xbdry);
        ybdry.myboundary = 1;
        ybdry.pairedboundary = 3;
        dof_map.add_periodic_boundary(ybdry);
        zbdry.myboundary = 0;
        zbdry.pairedboundary = 5;
        dof_map.add_periodic_boundary(zbdry);
        <B><FONT COLOR="#A020F0">break</FONT></B>;
  #endif
      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
        libmesh_error();
      }
  #endif <I><FONT COLOR="#B22222">// LIBMESH_ENABLE_PERIODIC
</FONT></I>  
  
    add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, THIRD, HERMITE);
  
    attach_init_object(*<B><FONT COLOR="#A020F0">this</FONT></B>);
  
    nonlinear_solver-&gt;residual_and_jacobian_object = <B><FONT COLOR="#A020F0">this</FONT></B>;
  
    nonlinear_solver-&gt;bounds_object = <B><FONT COLOR="#A020F0">this</FONT></B>;
  }
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::JR::initialize()
  {
    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&gt;&gt;&gt; Initializing Biharmonic::JR\n&quot;</FONT></B>;
  
    Parameters parameters;
    parameters.set&lt;Point&gt;(<B><FONT COLOR="#BC8F8F">&quot;center&quot;</FONT></B>) = _biharmonic._initialCenter;
    parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;width&quot;</FONT></B>)   = _biharmonic._initialWidth;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._initialState == Biharmonic::BALL)
      project_solution(Biharmonic::JR::InitialDensityBall, Biharmonic::JR::InitialGradientZero, parameters);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._initialState == Biharmonic::ROD)
      project_solution(Biharmonic::JR::InitialDensityRod, Biharmonic::JR::InitialGradientZero, parameters);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._initialState == Biharmonic::STRIP)
      project_solution(Biharmonic::JR::InitialDensityStrip, Biharmonic::JR::InitialGradientZero, parameters);
  
    *(old_local_solution) = *(current_local_solution);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._verbose)
      <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;&lt;&lt;&lt; Initializing Biharmonic::JR\n&quot;</FONT></B>;
  }
  
  
  
  
  
  
  Number Biharmonic::JR::InitialDensityBall(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  					  <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
  					  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
  					  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Point center = parameters.get&lt;Point&gt;(<B><FONT COLOR="#BC8F8F">&quot;center&quot;</FONT></B>);
    Real width = parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;width&quot;</FONT></B>);
    Point pc = p-center;
    Real r = pc.size();
    <B><FONT COLOR="#A020F0">return</FONT></B> (r &lt; width) ? 1.0 : -0.5;
  }
  
  
  
  
  Number Biharmonic::JR::InitialDensityRod(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  					 <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
  					 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
  					 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Point center = parameters.get&lt;Point&gt;(<B><FONT COLOR="#BC8F8F">&quot;center&quot;</FONT></B>);
    Real width = parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;width&quot;</FONT></B>);
    Real r = sqrt((p(0)-center(0))*(p(0)-center(0)) + (p(1)-center(1))*(p(1)-center(1)));
    <B><FONT COLOR="#A020F0">return</FONT></B> (r &lt; width) ? 1.0 : -0.5;
  }
  
  
  
  
  
  Number Biharmonic::JR::InitialDensityStrip(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  					   <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
  					   <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
  					   <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    Point center = parameters.get&lt;Point&gt;(<B><FONT COLOR="#BC8F8F">&quot;center&quot;</FONT></B>);
    Real width = parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;width&quot;</FONT></B>);
    Real r = sqrt((p(0)-center(0))*(p(0)-center(0)));
    <B><FONT COLOR="#A020F0">return</FONT></B> (r &lt; width) ? 1.0 : -0.5;
  }
  
  
  
  
  Gradient Biharmonic::JR::InitialGradientZero(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp;,
  					     <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,
  					     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
  					     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> Gradient(0.0,0.0,0.0);
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::JR::residual_and_jacobian(<B><FONT COLOR="#228B22">const</FONT></B> NumericVector&lt;Number&gt; &amp;u,
  					   NumericVector&lt;Number&gt; *R,
  					   SparseMatrix&lt;Number&gt; *J,
  					   NonlinearImplicitSystem&amp;)
  {
  #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
    <B><FONT COLOR="#A020F0">if</FONT></B> (!R &amp;&amp; !J)
      <B><FONT COLOR="#A020F0">return</FONT></B>;
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Biharmonic Residual and Jacobian&quot;</FONT></B>, false);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(_biharmonic._dim, fe_type));
  
    AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(_biharmonic._dim));
  
    fe-&gt;attach_quadrature_rule (qrule.get());
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealTensor&gt; &gt;&amp; d2phi = fe-&gt;get_d2phi();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; Laplacian_phi_qp;
  
    DenseMatrix&lt;Number&gt; Je;
    DenseVector&lt;Number&gt; Re;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#228B22">const</FONT></B> NumericVector&lt;Number&gt;&amp; u_old = *old_local_solution;
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = _biharmonic._mesh-&gt;active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = _biharmonic._mesh-&gt;active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el) {
      <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
      dof_map.dof_indices (elem, dof_indices);
  
      fe-&gt;reinit (elem);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (J)
        Je.resize(dof_indices.size(), dof_indices.size());
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (R)
        Re.resize(dof_indices.size());
  
      Laplacian_phi_qp.resize(dof_indices.size());
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
        {
  	Number u_qp = 0.0, u_old_qp = 0.0, Laplacian_u_qp = 0.0, Laplacian_u_old_qp = 0.0;
  	Gradient grad_u_qp(0.0,0.0,0.0), grad_u_old_qp(0.0,0.0,0.0);
  	Number M_qp = 1.0, M_old_qp = 1.0, M_prime_qp = 0.0, M_prime_old_qp = 0.0;
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  {
  	    Laplacian_phi_qp[i] = d2phi[i][qp](0,0);
  	    grad_u_qp(0) += u(dof_indices[i])*dphi[i][qp](0);
  	    grad_u_old_qp(0) += u_old(dof_indices[i])*dphi[i][qp](0);
  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._dim &gt; 1)
  	      {
  		Laplacian_phi_qp[i] += d2phi[i][qp](1,1);
  		grad_u_qp(1) += u(dof_indices[i])*dphi[i][qp](1);
  		grad_u_old_qp(1) += u_old(dof_indices[i])*dphi[i][qp](1);
  	      }
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._dim &gt; 2)
  	      {
  		Laplacian_phi_qp[i] += d2phi[i][qp](2,2);
  		grad_u_qp(2) += u(dof_indices[i])*dphi[i][qp](2);
  		grad_u_old_qp(2) += u_old(dof_indices[i])*dphi[i][qp](2);
  	      }
  	    u_qp     += phi[i][qp]*u(dof_indices[i]);
  	    u_old_qp += phi[i][qp]*u_old(dof_indices[i]);
  	    Laplacian_u_qp     += Laplacian_phi_qp[i]*u(dof_indices[i]);
  	    Laplacian_u_old_qp += Laplacian_phi_qp[i]*u_old(dof_indices[i]);
  	  } <I><FONT COLOR="#B22222">// for i
</FONT></I>  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._degenerate)
  	  {
  	    M_qp           = 1.0 - u_qp*u_qp;
  	    M_old_qp       = 1.0 - u_old_qp*u_old_qp;
  	    M_prime_qp     = -2.0*u_qp;
  	    M_prime_old_qp = -2.0*u_old_qp;
  	  }
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  {
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (R)
  	      {
  		Number ri = 0.0, ri_old = 0.0;
  		ri     -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_u_qp;
  		ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._kappa*Laplacian_u_old_qp;
  
  		<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._degenerate)
  		  {
  		    ri       -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp);
  		    ri_old   -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*(_biharmonic._kappa*Laplacian_u_old_qp);
  		  }
  
  		<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._cahn_hillard)
  		  {
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
  		      {
  			ri += Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
  			ri_old += Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
  			<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._degenerate)
  			  {
  			    ri     += (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
  			    ri_old += (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
  			  }
  		      }<I><FONT COLOR="#B22222">// if(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
</FONT></I>  
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
  		      {
  			ri -= Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*u_qp;
  			ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*u_old_qp;
  			<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._degenerate)
  			  {
  			    ri     -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*u_qp;
  			    ri_old -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*u_old_qp;
  			  }
  		      } <I><FONT COLOR="#B22222">// if(_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
</FONT></I>  
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
  		      {
  			<B><FONT COLOR="#A020F0">switch</FONT></B>(_biharmonic._log_truncation)
  			  {
  			  <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
  			    <B><FONT COLOR="#A020F0">break</FONT></B>;
  			  <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
  			    <B><FONT COLOR="#A020F0">break</FONT></B>;
  			  <B><FONT COLOR="#5F9EA0">default</FONT></B>:
  			    <B><FONT COLOR="#A020F0">break</FONT></B>;
  			  }<I><FONT COLOR="#B22222">// switch(_biharmonic._log_truncation)
</FONT></I>  		      }<I><FONT COLOR="#B22222">// if(_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
</FONT></I>  		  }<I><FONT COLOR="#B22222">// if(_biharmonic._cahn_hillard)
</FONT></I>  		Re(i) += JxW[qp]*((u_qp-u_old_qp)*phi[i][qp]-_biharmonic._dt*0.5*((2.0-_biharmonic._cnWeight)*ri + _biharmonic._cnWeight*ri_old));
  	      } <I><FONT COLOR="#B22222">// if (R)
</FONT></I>  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (J)
  	      {
  		Number M_prime_prime_qp = 0.0;
  		<B><FONT COLOR="#A020F0">if</FONT></B>(_biharmonic._degenerate) M_prime_prime_qp = -2.0;
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
  		  {
  		    Number ri_j = 0.0;
  		    ri_j -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_phi_qp[j];
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._degenerate)
  		      {
  			ri_j -=
  			  Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._kappa*Laplacian_u_qp               +
  			  (dphi[i][qp]*dphi[j][qp])*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp)                  +
  			  (dphi[i][qp]*grad_u_qp)*(M_prime_prime_qp*phi[j][qp])*(_biharmonic._kappa*Laplacian_u_qp) +
  			  (dphi[i][qp]*grad_u_qp)*(M_prime_qp)*(_biharmonic._kappa*Laplacian_phi_qp[j]);
  		      }
  
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._cahn_hillard)
  		      {
  			<B><FONT COLOR="#A020F0">if</FONT></B>(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
  			  {
  			    ri_j +=
  			      Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp +
  			      Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp]        +
  			      (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp      +
  			      (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp  +
  			      (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp];
  			  }<I><FONT COLOR="#B22222">// if(_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
</FONT></I>  
  			<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
  			  {
  			    ri_j -=
  			      Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*u_qp                   +
  			      Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*phi[j][qp]                              +
  			      (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*u_qp                        +
  			      (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*u_qp                    +
  			      (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*phi[j][qp];
  			  } <I><FONT COLOR="#B22222">// if(_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
</FONT></I>  
  			<B><FONT COLOR="#A020F0">if</FONT></B> (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
  			  {
  			    <B><FONT COLOR="#A020F0">switch</FONT></B>(_biharmonic._log_truncation)
  			      {
  			      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">2</FONT></B>:
  				<B><FONT COLOR="#A020F0">break</FONT></B>;
  			      <B><FONT COLOR="#A020F0">case</FONT></B> <B><FONT COLOR="#5F9EA0">3</FONT></B>:
  				<B><FONT COLOR="#A020F0">break</FONT></B>;
  			      <B><FONT COLOR="#5F9EA0">default</FONT></B>:
  				<B><FONT COLOR="#A020F0">break</FONT></B>;
  			      }<I><FONT COLOR="#B22222">// switch(_biharmonic._log_truncation)
</FONT></I>  			  }<I><FONT COLOR="#B22222">// if(_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
</FONT></I>  		      }<I><FONT COLOR="#B22222">// if(_biharmonic._cahn_hillard)
</FONT></I>  		    Je(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] - 0.5*_biharmonic._dt*(2.0-_biharmonic._cnWeight)*ri_j);
  		  } <I><FONT COLOR="#B22222">// for j
</FONT></I>  	      } <I><FONT COLOR="#B22222">// if (J)
</FONT></I>  	  } <I><FONT COLOR="#B22222">// for i
</FONT></I>        } <I><FONT COLOR="#B22222">// for qp
</FONT></I>  
      <B><FONT COLOR="#A020F0">if</FONT></B> (R)
        {
  	dof_map.constrain_element_vector(Re, dof_indices);
  	R-&gt;add_vector(Re, dof_indices);
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (J)
        {
  	dof_map.constrain_element_matrix(Je, dof_indices);
  	J-&gt;add_matrix(Je, dof_indices);
        }
    } <I><FONT COLOR="#B22222">// for el
</FONT></I>  #endif <I><FONT COLOR="#B22222">// LIBMESH_ENABLE_SECOND_DERIVATIVES
</FONT></I>  }
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> Biharmonic::JR::bounds(NumericVector&lt;Number&gt; &amp;XL, NumericVector&lt;Number&gt;&amp; XU, NonlinearImplicitSystem&amp;)
  {
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Biharmonic bounds&quot;</FONT></B>, false);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(_biharmonic._dim, fe_type));
  
    DenseVector&lt;Number&gt; XLe, XUe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = _biharmonic._mesh-&gt;active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = _biharmonic._mesh-&gt;active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
        dof_map.dof_indices (*el, dof_indices);
  
        XLe.resize(dof_indices.size());
        XUe.resize(dof_indices.size());
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt; nodes;
        fe-&gt;get_refspace_nodes((*el)-&gt;type(), nodes);
  
        fe-&gt;reinit(*el, &amp;nodes);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); ++i)
  	{
  	  Real infinity = 1.0e20;
  	  Real bound = infinity;
  	  <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j = 0; j &lt; nodes.size(); ++j) {
  	    <B><FONT COLOR="#A020F0">if</FONT></B>(phi[i][j]) {
  	      bound = 1.0/fabs(phi[i][j]);
  	      <B><FONT COLOR="#A020F0">break</FONT></B>;
  	    }
  	  }
  
  	  XLe(i) = -bound;
  	  XUe(i) = bound;
  	}
        XL.insert(XLe, dof_indices);
        XU.insert(XUe, dof_indices);
      }
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex7.C without comments: </h1> 
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
        <B><FONT COLOR="#5F9EA0">Biharmonic</FONT></B>::Create(&amp;biharmonic,init.comm());
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
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex7'
***************************************************************
* Running Example miscellaneous_ex7:
*  mpirun -np 4 example-devel --verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-8 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
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
  NL step  1, |residual|_2 = 6.791353e-15
Writing state 2 at time 1.000000e-10 to file bih.1.e; output a total of 2 states so far
Solving for state 3, time 1.000000e-10
  NL step  0, |residual|_2 = 3.861870e+00
  NL step  1, |residual|_2 = 5.295328e-13
Writing state 3 at time 2.000000e-10 to file bih.1.e; output a total of 3 states so far
Solving for state 4, time 2.000000e-10
  NL step  0, |residual|_2 = 3.858634e+00
  NL step  1, |residual|_2 = 6.826469e-12
Writing state 4 at time 3.000000e-10 to file bih.1.e; output a total of 4 states so far
Solving for state 5, time 3.000000e-10
  NL step  0, |residual|_2 = 3.855483e+00
  NL step  1, |residual|_2 = 2.625789e-11
Writing state 5 at time 4.000000e-10 to file bih.1.e; output a total of 5 states so far
Solving for state 6, time 4.000000e-10
  NL step  0, |residual|_2 = 3.852464e+00
  NL step  1, |residual|_2 = 3.298126e-11
Writing state 6 at time 5.000000e-10 to file bih.1.e; output a total of 6 states so far
Solving for state 7, time 5.000000e-10
  NL step  0, |residual|_2 = 3.849501e+00
  NL step  1, |residual|_2 = 5.351723e-10
Writing state 7 at time 6.000000e-10 to file bih.1.e; output a total of 7 states so far
Solving for state 8, time 6.000000e-10
  NL step  0, |residual|_2 = 3.846621e+00
  NL step  1, |residual|_2 = 1.678735e-09
Writing state 8 at time 7.000000e-10 to file bih.1.e; output a total of 8 states so far
Solving for state 9, time 7.000000e-10
  NL step  0, |residual|_2 = 3.843784e+00
  NL step  1, |residual|_2 = 1.776743e-09
Writing state 9 at time 8.000000e-10 to file bih.1.e; output a total of 9 states so far
Solving for state 10, time 8.000000e-10
  NL step  0, |residual|_2 = 3.841009e+00
  NL step  1, |residual|_2 = 2.711450e-09
Writing state 10 at time 9.000000e-10 to file bih.1.e; output a total of 10 states so far
Solving for state 11, time 9.000000e-10
  NL step  0, |residual|_2 = 3.838268e+00
  NL step  1, |residual|_2 = 1.204330e-08
Writing state 11 at time 1.000000e-09 to file bih.1.e; output a total of 11 states so far
Solving for state 12, time 1.000000e-09
  NL step  0, |residual|_2 = 3.835577e+00
  NL step  1, |residual|_2 = 1.982070e-08
Writing state 12 at time 1.100000e-09 to file bih.1.e; output a total of 12 states so far
Solving for state 13, time 1.100000e-09
  NL step  0, |residual|_2 = 3.832915e+00
  NL step  1, |residual|_2 = 1.887360e-08
Writing state 13 at time 1.200000e-09 to file bih.1.e; output a total of 13 states so far
Solving for state 14, time 1.200000e-09
  NL step  0, |residual|_2 = 3.830295e+00
  NL step  1, |residual|_2 = 7.909848e-09
Writing state 14 at time 1.300000e-09 to file bih.1.e; output a total of 14 states so far
Solving for state 15, time 1.300000e-09
  NL step  0, |residual|_2 = 3.827700e+00
  NL step  1, |residual|_2 = 1.049956e-08
Writing state 15 at time 1.400000e-09 to file bih.1.e; output a total of 15 states so far
Solving for state 16, time 1.400000e-09
  NL step  0, |residual|_2 = 3.825140e+00
  NL step  1, |residual|_2 = 3.280321e-08
Writing state 16 at time 1.500000e-09 to file bih.1.e; output a total of 16 states so far
Solving for state 17, time 1.500000e-09
  NL step  0, |residual|_2 = 3.822604e+00
  NL step  1, |residual|_2 = 5.475702e-08
  NL step  2, |residual|_2 = 2.323455e-15
Writing state 17 at time 1.600000e-09 to file bih.1.e; output a total of 17 states so far
Solving for state 18, time 1.600000e-09
  NL step  0, |residual|_2 = 3.820098e+00
  NL step  1, |residual|_2 = 3.707085e-08
Writing state 18 at time 1.700000e-09 to file bih.1.e; output a total of 18 states so far
Solving for state 19, time 1.700000e-09
  NL step  0, |residual|_2 = 3.817613e+00
  NL step  1, |residual|_2 = 9.270300e-08
  NL step  2, |residual|_2 = 2.353413e-15
Writing state 19 at time 1.800000e-09 to file bih.1.e; output a total of 19 states so far
Solving for state 20, time 1.800000e-09
  NL step  0, |residual|_2 = 3.815156e+00
  NL step  1, |residual|_2 = 1.778004e-08
Writing state 20 at time 1.900000e-09 to file bih.1.e; output a total of 20 states so far
Solving for state 21, time 1.900000e-09
  NL step  0, |residual|_2 = 3.812718e+00
  NL step  1, |residual|_2 = 1.166306e-07
  NL step  2, |residual|_2 = 2.344579e-15
Writing state 21 at time 2.000000e-09 to file bih.1.e; output a total of 21 states so far
Solving for state 22, time 2.000000e-09
  NL step  0, |residual|_2 = 3.810305e+00
  NL step  1, |residual|_2 = 5.194027e-08
  NL step  2, |residual|_2 = 2.292010e-15
Writing state 22 at time 2.100000e-09 to file bih.1.e; output a total of 22 states so far
Solving for state 23, time 2.100000e-09
  NL step  0, |residual|_2 = 3.807909e+00
  NL step  1, |residual|_2 = 7.569577e-08
  NL step  2, |residual|_2 = 2.495208e-15
Writing state 23 at time 2.200000e-09 to file bih.1.e; output a total of 23 states so far
Solving for state 24, time 2.200000e-09
  NL step  0, |residual|_2 = 3.805536e+00
  NL step  1, |residual|_2 = 7.326139e-08
  NL step  2, |residual|_2 = 2.341764e-15
Writing state 24 at time 2.300000e-09 to file bih.1.e; output a total of 24 states so far
Solving for state 25, time 2.300000e-09
  NL step  0, |residual|_2 = 3.803179e+00
  NL step  1, |residual|_2 = 4.221623e-08
  NL step  2, |residual|_2 = 2.358325e-15
Writing state 25 at time 2.400000e-09 to file bih.1.e; output a total of 25 states so far
Solving for state 26, time 2.400000e-09
  NL step  0, |residual|_2 = 3.800842e+00
  NL step  1, |residual|_2 = 1.024370e-07
  NL step  2, |residual|_2 = 2.288452e-15
Writing state 26 at time 2.500000e-09 to file bih.1.e; output a total of 26 states so far
Solving for state 27, time 2.500000e-09
  NL step  0, |residual|_2 = 3.798522e+00
  NL step  1, |residual|_2 = 3.410027e-08
Writing state 27 at time 2.600000e-09 to file bih.1.e; output a total of 27 states so far
Solving for state 28, time 2.600000e-09
  NL step  0, |residual|_2 = 3.796220e+00
  NL step  1, |residual|_2 = 2.482465e-07
  NL step  2, |residual|_2 = 2.646697e-15
Writing state 28 at time 2.700000e-09 to file bih.1.e; output a total of 28 states so far
Solving for state 29, time 2.700000e-09
  NL step  0, |residual|_2 = 3.793932e+00
  NL step  1, |residual|_2 = 1.025109e-07
  NL step  2, |residual|_2 = 2.243753e-15
Writing state 29 at time 2.800000e-09 to file bih.1.e; output a total of 29 states so far
Solving for state 30, time 2.800000e-09
  NL step  0, |residual|_2 = 3.791663e+00
  NL step  1, |residual|_2 = 2.505431e-07
  NL step  2, |residual|_2 = 2.438573e-15
Writing state 30 at time 2.900000e-09 to file bih.1.e; output a total of 30 states so far
Solving for state 31, time 2.900000e-09
  NL step  0, |residual|_2 = 3.789407e+00
  NL step  1, |residual|_2 = 9.550927e-08
  NL step  2, |residual|_2 = 2.429171e-15
Writing state 31 at time 3.000000e-09 to file bih.1.e; output a total of 31 states so far
Solving for state 32, time 3.000000e-09
  NL step  0, |residual|_2 = 3.787167e+00
  NL step  1, |residual|_2 = 2.581924e-07
  NL step  2, |residual|_2 = 2.508660e-15
Writing state 32 at time 3.100000e-09 to file bih.1.e; output a total of 32 states so far
Solving for state 33, time 3.100000e-09
  NL step  0, |residual|_2 = 3.784941e+00
  NL step  1, |residual|_2 = 9.692862e-08
  NL step  2, |residual|_2 = 2.404583e-15
Writing state 33 at time 3.200000e-09 to file bih.1.e; output a total of 33 states so far
Solving for state 34, time 3.200000e-09
  NL step  0, |residual|_2 = 3.782730e+00
  NL step  1, |residual|_2 = 2.594179e-07
  NL step  2, |residual|_2 = 2.491096e-15
Writing state 34 at time 3.300000e-09 to file bih.1.e; output a total of 34 states so far
Solving for state 35, time 3.300000e-09
  NL step  0, |residual|_2 = 3.780531e+00
  NL step  1, |residual|_2 = 9.887434e-08
  NL step  2, |residual|_2 = 2.479266e-15
Writing state 35 at time 3.400000e-09 to file bih.1.e; output a total of 35 states so far
Solving for state 36, time 3.400000e-09
  NL step  0, |residual|_2 = 3.778347e+00
  NL step  1, |residual|_2 = 2.527746e-07
  NL step  2, |residual|_2 = 2.569841e-15
Writing state 36 at time 3.500000e-09 to file bih.1.e; output a total of 36 states so far
Solving for state 37, time 3.500000e-09
  NL step  0, |residual|_2 = 3.776174e+00
  NL step  1, |residual|_2 = 9.893108e-08
  NL step  2, |residual|_2 = 2.460867e-15
Writing state 37 at time 3.600000e-09 to file bih.1.e; output a total of 37 states so far
Solving for state 38, time 3.600000e-09
  NL step  0, |residual|_2 = 3.774015e+00
  NL step  1, |residual|_2 = 2.386953e-07
  NL step  2, |residual|_2 = 2.287197e-15
Writing state 38 at time 3.700000e-09 to file bih.1.e; output a total of 38 states so far
Solving for state 39, time 3.700000e-09
  NL step  0, |residual|_2 = 3.771868e+00
  NL step  1, |residual|_2 = 9.823210e-08
  NL step  2, |residual|_2 = 2.425269e-15
Writing state 39 at time 3.800000e-09 to file bih.1.e; output a total of 39 states so far
Solving for state 40, time 3.800000e-09
  NL step  0, |residual|_2 = 3.769733e+00
  NL step  1, |residual|_2 = 2.184248e-07
  NL step  2, |residual|_2 = 2.412187e-15
Writing state 40 at time 3.900000e-09 to file bih.1.e; output a total of 40 states so far
Solving for state 41, time 3.900000e-09
  NL step  0, |residual|_2 = 3.767609e+00
  NL step  1, |residual|_2 = 9.982364e-08
  NL step  2, |residual|_2 = 2.594741e-15
Writing state 41 at time 4.000000e-09 to file bih.1.e; output a total of 41 states so far
Solving for state 42, time 4.000000e-09
  NL step  0, |residual|_2 = 3.765498e+00
  NL step  1, |residual|_2 = 1.936364e-07
  NL step  2, |residual|_2 = 2.566911e-15
Writing state 42 at time 4.100000e-09 to file bih.1.e; output a total of 42 states so far
Solving for state 43, time 4.100000e-09
  NL step  0, |residual|_2 = 3.763397e+00
  NL step  1, |residual|_2 = 1.069263e-07
  NL step  2, |residual|_2 = 2.410036e-15
Writing state 43 at time 4.200000e-09 to file bih.1.e; output a total of 43 states so far
Solving for state 44, time 4.200000e-09
  NL step  0, |residual|_2 = 3.761307e+00
  NL step  1, |residual|_2 = 1.662948e-07
  NL step  2, |residual|_2 = 2.717666e-15
Writing state 44 at time 4.300000e-09 to file bih.1.e; output a total of 44 states so far
Solving for state 45, time 4.300000e-09
  NL step  0, |residual|_2 = 3.759228e+00
  NL step  1, |residual|_2 = 1.210237e-07
  NL step  2, |residual|_2 = 2.503640e-15
Writing state 45 at time 4.400000e-09 to file bih.1.e; output a total of 45 states so far
Solving for state 46, time 4.400000e-09
  NL step  0, |residual|_2 = 3.757160e+00
  NL step  1, |residual|_2 = 1.387212e-07
  NL step  2, |residual|_2 = 2.443778e-15
Writing state 46 at time 4.500000e-09 to file bih.1.e; output a total of 46 states so far
Solving for state 47, time 4.500000e-09
  NL step  0, |residual|_2 = 3.755101e+00
  NL step  1, |residual|_2 = 1.413200e-07
  NL step  2, |residual|_2 = 2.453897e-15
Writing state 47 at time 4.600000e-09 to file bih.1.e; output a total of 47 states so far
Solving for state 48, time 4.600000e-09
  NL step  0, |residual|_2 = 3.753053e+00
  NL step  1, |residual|_2 = 1.138969e-07
  NL step  2, |residual|_2 = 2.520689e-15
Writing state 48 at time 4.700000e-09 to file bih.1.e; output a total of 48 states so far
Solving for state 49, time 4.700000e-09
  NL step  0, |residual|_2 = 3.751014e+00
  NL step  1, |residual|_2 = 1.658550e-07
  NL step  2, |residual|_2 = 2.425168e-15
Writing state 49 at time 4.800000e-09 to file bih.1.e; output a total of 49 states so far
Solving for state 50, time 4.800000e-09
  NL step  0, |residual|_2 = 3.748985e+00
  NL step  1, |residual|_2 = 9.588393e-08
  NL step  2, |residual|_2 = 2.380361e-15
Writing state 50 at time 4.900000e-09 to file bih.1.e; output a total of 50 states so far
Solving for state 51, time 4.900000e-09
  NL step  0, |residual|_2 = 3.746966e+00
  NL step  1, |residual|_2 = 1.926354e-07
  NL step  2, |residual|_2 = 2.380779e-15
Writing state 51 at time 5.000000e-09 to file bih.1.e; output a total of 51 states so far
Solving for state 52, time 5.000000e-09
  NL step  0, |residual|_2 = 3.744956e+00
  NL step  1, |residual|_2 = 8.914698e-08
  NL step  2, |residual|_2 = 2.400894e-15
Writing state 52 at time 5.100000e-09 to file bih.1.e; output a total of 52 states so far
Solving for state 53, time 5.100000e-09
  NL step  0, |residual|_2 = 3.742955e+00
  NL step  1, |residual|_2 = 2.200581e-07
  NL step  2, |residual|_2 = 2.622921e-15
Writing state 53 at time 5.200000e-09 to file bih.1.e; output a total of 53 states so far
Solving for state 54, time 5.200000e-09
  NL step  0, |residual|_2 = 3.740963e+00
  NL step  1, |residual|_2 = 9.501731e-08
  NL step  2, |residual|_2 = 2.551815e-15
Writing state 54 at time 5.300000e-09 to file bih.1.e; output a total of 54 states so far
Solving for state 55, time 5.300000e-09
  NL step  0, |residual|_2 = 3.738979e+00
  NL step  1, |residual|_2 = 2.469333e-07
  NL step  2, |residual|_2 = 2.529382e-15
Writing state 55 at time 5.400000e-09 to file bih.1.e; output a total of 55 states so far
Solving for state 56, time 5.400000e-09
  NL step  0, |residual|_2 = 3.737005e+00
  NL step  1, |residual|_2 = 1.098357e-07
  NL step  2, |residual|_2 = 2.611139e-15
Writing state 56 at time 5.500000e-09 to file bih.1.e; output a total of 56 states so far
Solving for state 57, time 5.500000e-09
  NL step  0, |residual|_2 = 3.735039e+00
  NL step  1, |residual|_2 = 2.724056e-07
  NL step  2, |residual|_2 = 2.504157e-15
Writing state 57 at time 5.600000e-09 to file bih.1.e; output a total of 57 states so far
Solving for state 58, time 5.600000e-09
  NL step  0, |residual|_2 = 3.733082e+00
  NL step  1, |residual|_2 = 1.287789e-07
  NL step  2, |residual|_2 = 2.315606e-15
Writing state 58 at time 5.700000e-09 to file bih.1.e; output a total of 58 states so far
Solving for state 59, time 5.700000e-09
  NL step  0, |residual|_2 = 3.731132e+00
  NL step  1, |residual|_2 = 2.958737e-07
  NL step  2, |residual|_2 = 2.431640e-15
Writing state 59 at time 5.800000e-09 to file bih.1.e; output a total of 59 states so far
Solving for state 60, time 5.800000e-09
  NL step  0, |residual|_2 = 3.729191e+00
  NL step  1, |residual|_2 = 1.486029e-07
  NL step  2, |residual|_2 = 2.465803e-15
Writing state 60 at time 5.900000e-09 to file bih.1.e; output a total of 60 states so far
Solving for state 61, time 5.900000e-09
  NL step  0, |residual|_2 = 3.727257e+00
  NL step  1, |residual|_2 = 3.169299e-07
  NL step  2, |residual|_2 = 2.622898e-15
Writing state 61 at time 6.000000e-09 to file bih.1.e; output a total of 61 states so far
Solving for state 62, time 6.000000e-09
  NL step  0, |residual|_2 = 3.725332e+00
  NL step  1, |residual|_2 = 1.675296e-07
  NL step  2, |residual|_2 = 2.414417e-15
Writing state 62 at time 6.100000e-09 to file bih.1.e; output a total of 62 states so far
Solving for state 63, time 6.100000e-09
  NL step  0, |residual|_2 = 3.723414e+00
  NL step  1, |residual|_2 = 3.353144e-07
  NL step  2, |residual|_2 = 2.456989e-15
Writing state 63 at time 6.200000e-09 to file bih.1.e; output a total of 63 states so far
Solving for state 64, time 6.200000e-09
  NL step  0, |residual|_2 = 3.721504e+00
  NL step  1, |residual|_2 = 1.846217e-07
  NL step  2, |residual|_2 = 2.613437e-15
Writing state 64 at time 6.300000e-09 to file bih.1.e; output a total of 64 states so far
Solving for state 65, time 6.300000e-09
  NL step  0, |residual|_2 = 3.719601e+00
  NL step  1, |residual|_2 = 3.508814e-07
  NL step  2, |residual|_2 = 2.482544e-15
Writing state 65 at time 6.400000e-09 to file bih.1.e; output a total of 65 states so far
Solving for state 66, time 6.400000e-09
  NL step  0, |residual|_2 = 3.717706e+00
  NL step  1, |residual|_2 = 1.993937e-07
  NL step  2, |residual|_2 = 2.434944e-15
Writing state 66 at time 6.500000e-09 to file bih.1.e; output a total of 66 states so far
Solving for state 67, time 6.500000e-09
  NL step  0, |residual|_2 = 3.715818e+00
  NL step  1, |residual|_2 = 3.635735e-07
  NL step  2, |residual|_2 = 2.556709e-15
Writing state 67 at time 6.600000e-09 to file bih.1.e; output a total of 67 states so far
Solving for state 68, time 6.600000e-09
  NL step  0, |residual|_2 = 3.713937e+00
  NL step  1, |residual|_2 = 2.116130e-07
  NL step  2, |residual|_2 = 2.473761e-15
Writing state 68 at time 6.700000e-09 to file bih.1.e; output a total of 68 states so far
Solving for state 69, time 6.700000e-09
  NL step  0, |residual|_2 = 3.712063e+00
  NL step  1, |residual|_2 = 3.734010e-07
  NL step  2, |residual|_2 = 2.527880e-15
Writing state 69 at time 6.800000e-09 to file bih.1.e; output a total of 69 states so far
Solving for state 70, time 6.800000e-09
  NL step  0, |residual|_2 = 3.710197e+00
  NL step  1, |residual|_2 = 2.211980e-07
  NL step  2, |residual|_2 = 2.400609e-15
Writing state 70 at time 6.900000e-09 to file bih.1.e; output a total of 70 states so far
Solving for state 71, time 6.900000e-09
  NL step  0, |residual|_2 = 3.708336e+00
  NL step  1, |residual|_2 = 3.804259e-07
  NL step  2, |residual|_2 = 2.394281e-15
Writing state 71 at time 7.000000e-09 to file bih.1.e; output a total of 71 states so far
Solving for state 72, time 7.000000e-09
  NL step  0, |residual|_2 = 3.706484e+00
  NL step  1, |residual|_2 = 2.281627e-07
  NL step  2, |residual|_2 = 2.563391e-15
Writing state 72 at time 7.100000e-09 to file bih.1.e; output a total of 72 states so far
Solving for state 73, time 7.100000e-09
  NL step  0, |residual|_2 = 3.704637e+00
  NL step  1, |residual|_2 = 3.847491e-07
  NL step  2, |residual|_2 = 2.486401e-15
Writing state 73 at time 7.200000e-09 to file bih.1.e; output a total of 73 states so far
Solving for state 74, time 7.200000e-09
  NL step  0, |residual|_2 = 3.702797e+00
  NL step  1, |residual|_2 = 2.325837e-07
  NL step  2, |residual|_2 = 2.406717e-15
Writing state 74 at time 7.300000e-09 to file bih.1.e; output a total of 74 states so far
Solving for state 75, time 7.300000e-09
  NL step  0, |residual|_2 = 3.700963e+00
  NL step  1, |residual|_2 = 3.864999e-07
  NL step  2, |residual|_2 = 2.434656e-15
Writing state 75 at time 7.400000e-09 to file bih.1.e; output a total of 75 states so far
Solving for state 76, time 7.400000e-09
  NL step  0, |residual|_2 = 3.699137e+00
  NL step  1, |residual|_2 = 2.345800e-07
  NL step  2, |residual|_2 = 2.600773e-15
Writing state 76 at time 7.500000e-09 to file bih.1.e; output a total of 76 states so far
Solving for state 77, time 7.500000e-09
  NL step  0, |residual|_2 = 3.697316e+00
  NL step  1, |residual|_2 = 3.858276e-07
  NL step  2, |residual|_2 = 2.479296e-15
Writing state 77 at time 7.600000e-09 to file bih.1.e; output a total of 77 states so far
Solving for state 78, time 7.600000e-09
  NL step  0, |residual|_2 = 3.695502e+00
  NL step  1, |residual|_2 = 2.342992e-07
  NL step  2, |residual|_2 = 2.388925e-15
Writing state 78 at time 7.700000e-09 to file bih.1.e; output a total of 78 states so far
Solving for state 79, time 7.700000e-09
  NL step  0, |residual|_2 = 3.693694e+00
  NL step  1, |residual|_2 = 3.828952e-07
  NL step  2, |residual|_2 = 2.600408e-15
Writing state 79 at time 7.800000e-09 to file bih.1.e; output a total of 79 states so far
Solving for state 80, time 7.800000e-09
  NL step  0, |residual|_2 = 3.691892e+00
  NL step  1, |residual|_2 = 2.319080e-07
  NL step  2, |residual|_2 = 2.659483e-15
Writing state 80 at time 7.900000e-09 to file bih.1.e; output a total of 80 states so far
Solving for state 81, time 7.900000e-09
  NL step  0, |residual|_2 = 3.690096e+00
  NL step  1, |residual|_2 = 3.778738e-07
  NL step  2, |residual|_2 = 2.615159e-15
Writing state 81 at time 8.000000e-09 to file bih.1.e; output a total of 81 states so far
Solving for state 82, time 8.000000e-09
  NL step  0, |residual|_2 = 3.688306e+00
  NL step  1, |residual|_2 = 2.275868e-07
  NL step  2, |residual|_2 = 2.439835e-15
Writing state 82 at time 8.100000e-09 to file bih.1.e; output a total of 82 states so far
Solving for state 83, time 8.100000e-09
  NL step  0, |residual|_2 = 3.686522e+00
  NL step  1, |residual|_2 = 3.709386e-07
  NL step  2, |residual|_2 = 2.363289e-15
Writing state 83 at time 8.200000e-09 to file bih.1.e; output a total of 83 states so far
Solving for state 84, time 8.200000e-09
  NL step  0, |residual|_2 = 3.684744e+00
  NL step  1, |residual|_2 = 2.215246e-07
  NL step  2, |residual|_2 = 2.408869e-15
Writing state 84 at time 8.300000e-09 to file bih.1.e; output a total of 84 states so far
Solving for state 85, time 8.300000e-09
  NL step  0, |residual|_2 = 3.682971e+00
  NL step  1, |residual|_2 = 3.622657e-07
  NL step  2, |residual|_2 = 2.323864e-15
Writing state 85 at time 8.400000e-09 to file bih.1.e; output a total of 85 states so far
Solving for state 86, time 8.400000e-09
  NL step  0, |residual|_2 = 3.681205e+00
  NL step  1, |residual|_2 = 2.139174e-07
  NL step  2, |residual|_2 = 2.295697e-15
Writing state 86 at time 8.500000e-09 to file bih.1.e; output a total of 86 states so far
Solving for state 87, time 8.500000e-09
  NL step  0, |residual|_2 = 3.679443e+00
  NL step  1, |residual|_2 = 3.520296e-07
  NL step  2, |residual|_2 = 2.376511e-15
Writing state 87 at time 8.600000e-09 to file bih.1.e; output a total of 87 states so far
Solving for state 88, time 8.600000e-09
  NL step  0, |residual|_2 = 3.677688e+00
  NL step  1, |residual|_2 = 2.049672e-07
  NL step  2, |residual|_2 = 2.508000e-15
Writing state 88 at time 8.700000e-09 to file bih.1.e; output a total of 88 states so far
Solving for state 89, time 8.700000e-09
  NL step  0, |residual|_2 = 3.675938e+00
  NL step  1, |residual|_2 = 3.404019e-07
  NL step  2, |residual|_2 = 2.459481e-15
Writing state 89 at time 8.800000e-09 to file bih.1.e; output a total of 89 states so far
Solving for state 90, time 8.800000e-09
  NL step  0, |residual|_2 = 3.674193e+00
  NL step  1, |residual|_2 = 1.948828e-07
  NL step  2, |residual|_2 = 2.472973e-15
Writing state 90 at time 8.900000e-09 to file bih.1.e; output a total of 90 states so far
Solving for state 91, time 8.900000e-09
  NL step  0, |residual|_2 = 3.672454e+00
  NL step  1, |residual|_2 = 3.275496e-07
  NL step  2, |residual|_2 = 2.655535e-15
Writing state 91 at time 9.000000e-09 to file bih.1.e; output a total of 91 states so far
Solving for state 92, time 9.000000e-09
  NL step  0, |residual|_2 = 3.670720e+00
  NL step  1, |residual|_2 = 1.838827e-07
  NL step  2, |residual|_2 = 2.526653e-15
Writing state 92 at time 9.100000e-09 to file bih.1.e; output a total of 92 states so far
Solving for state 93, time 9.100000e-09
  NL step  0, |residual|_2 = 3.668992e+00
  NL step  1, |residual|_2 = 3.136350e-07
  NL step  2, |residual|_2 = 2.328935e-15
Writing state 93 at time 9.200000e-09 to file bih.1.e; output a total of 93 states so far
Solving for state 94, time 9.200000e-09
  NL step  0, |residual|_2 = 3.667269e+00
  NL step  1, |residual|_2 = 1.722002e-07
  NL step  2, |residual|_2 = 2.247037e-15
Writing state 94 at time 9.300000e-09 to file bih.1.e; output a total of 94 states so far
Solving for state 95, time 9.300000e-09
  NL step  0, |residual|_2 = 3.665550e+00
  NL step  1, |residual|_2 = 2.988151e-07
  NL step  2, |residual|_2 = 2.459460e-15
Writing state 95 at time 9.400000e-09 to file bih.1.e; output a total of 95 states so far
Solving for state 96, time 9.400000e-09
  NL step  0, |residual|_2 = 3.663837e+00
  NL step  1, |residual|_2 = 1.600913e-07
  NL step  2, |residual|_2 = 2.484313e-15
Writing state 96 at time 9.500000e-09 to file bih.1.e; output a total of 96 states so far
Solving for state 97, time 9.500000e-09
  NL step  0, |residual|_2 = 3.662129e+00
  NL step  1, |residual|_2 = 2.832423e-07
  NL step  2, |residual|_2 = 2.481569e-15
Writing state 97 at time 9.600000e-09 to file bih.1.e; output a total of 97 states so far
Solving for state 98, time 9.600000e-09
  NL step  0, |residual|_2 = 3.660426e+00
  NL step  1, |residual|_2 = 1.478468e-07
  NL step  2, |residual|_2 = 2.436396e-15
Writing state 98 at time 9.700000e-09 to file bih.1.e; output a total of 98 states so far
Solving for state 99, time 9.700000e-09
  NL step  0, |residual|_2 = 3.658728e+00
  NL step  1, |residual|_2 = 2.670648e-07
  NL step  2, |residual|_2 = 2.318336e-15
Writing state 99 at time 9.800000e-09 to file bih.1.e; output a total of 99 states so far
Solving for state 100, time 9.800000e-09
  NL step  0, |residual|_2 = 3.657035e+00
  NL step  1, |residual|_2 = 1.358089e-07
  NL step  2, |residual|_2 = 2.353539e-15
Writing state 100 at time 9.900000e-09 to file bih.1.e; output a total of 100 states so far
Solving for state 101, time 9.900000e-09
  NL step  0, |residual|_2 = 3.655347e+00
  NL step  1, |residual|_2 = 2.504283e-07
  NL step  2, |residual|_2 = 2.337847e-15
Writing state 101 at time 1.000000e-08 to file bih.1.e; output a total of 101 states so far
Solving for state 102, time 1.000000e-08
  NL step  0, |residual|_2 = 3.653664e+00
  NL step  1, |residual|_2 = 1.243929e-07
  NL step  2, |residual|_2 = 2.516733e-15
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 102 states so far
Writing state 102 at time 1.010000e-08 to file bih.1.e; output a total of 103 states so far

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:51:51 2013                                                                          |
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
| libMesh Performance: Alive time=12.2918, Active time=8.91372                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0010      0.001048    0.0011      0.001063    0.01     0.01     |
|   build_constraint_matrix()        120064    0.1420      0.000001    0.1420      0.000001    1.59     1.59     |
|   build_sparsity()                 1         0.0009      0.000905    0.0025      0.002509    0.01     0.03     |
|   constrain_elem_matrix()          47104     0.0514      0.000001    0.0514      0.000001    0.58     0.58     |
|   constrain_elem_vector()          72960     0.0755      0.000001    0.0755      0.000001    0.85     0.85     |
|   create_dof_constraints()         1         0.0121      0.012083    0.0191      0.019082    0.14     0.21     |
|   distribute_dofs()                1         0.0038      0.003789    0.0081      0.008060    0.04     0.09     |
|   dof_indices()                    146950    0.8103      0.000006    0.8103      0.000006    9.09     9.09     |
|   enforce_constraints_exactly()    470       0.1099      0.000234    0.1099      0.000234    1.23     1.23     |
|   prepare_send_list()              1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                         1         0.0039      0.003938    0.0039      0.003938    0.04     0.04     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          103       0.2317      0.002250    0.4340      0.004214    2.60     4.87     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               103       5.3463      0.051906    5.3463      0.051906    59.98    59.98    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        473       0.0079      0.000017    0.0079      0.000017    0.09     0.09     |
|   init_shape_functions()           473       0.0165      0.000035    0.0165      0.000035    0.19     0.19     |
|   inverse_map()                    17        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             120068    0.2075      0.000002    0.2075      0.000002    2.33     2.33     |
|   compute_face_map()               2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_face_shape_functions()      2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_reference_to_physical_map() 473       0.0023      0.000005    0.0023      0.000005    0.03     0.03     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0016      0.001646    0.0018      0.001840    0.02     0.02     |
|   renumber_nodes_and_elem()        2         0.0001      0.000062    0.0001      0.000062    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0085      0.004269    0.0085      0.004269    0.10     0.10     |
|   find_global_indices()            2         0.0013      0.000633    0.0261      0.013068    0.01     0.29     |
|   parallel_sort()                  2         0.0045      0.002233    0.0144      0.007190    0.05     0.16     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         103       0.0027      0.000026    5.7906      0.056219    0.03     64.96    |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0006      0.000617    0.0006      0.000617    0.01     0.01     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0035      0.003480    0.0176      0.017597    0.04     0.20     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      7         0.0010      0.000138    0.0013      0.000189    0.01     0.01     |
|   max(bool)                        1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   max(scalar)                      5391      0.0257      0.000005    0.0257      0.000005    0.29     0.29     |
|   max(vector)                      1271      0.0158      0.000012    0.0336      0.000026    0.18     0.38     |
|   min(bool)                        6656      0.0297      0.000004    0.0297      0.000004    0.33     0.33     |
|   min(scalar)                      5389      0.0416      0.000008    0.0416      0.000008    0.47     0.47     |
|   min(vector)                      1271      0.0123      0.000010    0.0299      0.000023    0.14     0.34     |
|   probe()                          36        0.0014      0.000040    0.0014      0.000040    0.02     0.02     |
|   receive()                        36        0.0001      0.000003    0.0016      0.000044    0.00     0.02     |
|   send()                           36        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   send_receive()                   40        0.0002      0.000005    0.0019      0.000048    0.00     0.02     |
|   sum()                            785       0.0225      0.000029    0.0401      0.000051    0.25     0.45     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           36        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0007      0.000737    0.0063      0.006303    0.01     0.07     |
|   set_parent_processor_ids()       1         0.0003      0.000294    0.0003      0.000294    0.00     0.00     |
|                                                                                                                |
| PetscNonlinearSolver                                                                                           |
|   jacobian()                       184       0.6562      0.003566    1.2088      0.006569    7.36     13.56    |
|   residual()                       285       0.8731      0.003064    1.7258      0.006055    9.80     19.36    |
|   solve()                          101       0.1751      0.001734    3.1098      0.030790    1.96     34.89    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  1         0.0064      0.006400    0.0066      0.006645    0.07     0.07     |
|   operator()                       4         0.0001      0.000022    0.0001      0.000027    0.00     0.00     |
|                                                                                                                |
| System                                                                                                         |
|   project_vector()                 1         0.0033      0.003339    0.0053      0.005274    0.04     0.06     |
|   solve()                          101       0.0016      0.000015    3.1113      0.030805    0.02     34.91    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            531016    8.9137                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex7:
*  mpirun -np 4 example-devel --verbose dim=1 N=1024 initial_state=strip initial_center=0.5 initial_width=0.1 dt=1e-10 max_time=1e-8 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex7'
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
