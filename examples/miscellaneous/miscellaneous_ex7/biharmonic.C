// Libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/replicated_mesh.h"

// Example includes
#include "biharmonic.h"
#include "biharmonic_jr.h"

using namespace libMesh;

// Constructor
Biharmonic::Biharmonic(ReplicatedMesh & mesh) :
  EquationSystems(mesh),
  _mesh(mesh)
{
  // Retrieve parameters and set defaults
  _verbose      = false;
  _growth       = false;
  _degenerate   = false;
  _cahn_hillard = false;
  _netforce     = false;

  if (on_command_line("--verbose"))
    _verbose = true;
  if (on_command_line("--growth"))
    _growth = true;
  if (on_command_line("--degenerate"))
    _degenerate = true;
  if (on_command_line("--cahn_hillard"))
    _cahn_hillard = true;
  if (on_command_line("--netforce"))
    _netforce = true;

  _kappa = command_line_value("kappa", 1.0);

  // "type of energy (double well, double obstacle, logarithmic+double well, logarithmic+double obstacle)"
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
    libmesh_error_msg("Unknown energy type: " << energy);

  _tol     = command_line_value("tol", 1.0e-8);
  _theta   = command_line_value("theta", .001);
  _theta_c = command_line_value("theta_c", 1.0);

  // "order of log truncation (0=none, 2=quadratic, 3=cubic)"
  _log_truncation = command_line_value("log_truncation", 2);

  if (!_log_truncation)
    libMesh::out << "WARNING: no truncation is being used for the logarithmic free energy term.\nWARNING: division by zero possible!\n";


  // Dimension
  _dim = command_line_value("dim", 1);

  libmesh_assert_msg((_dim <= 3) && (_dim > 0), "Invalid mesh dimension");

  // Build the mesh
  // Yes, it's better to make a coarse mesh and then refine it. We'll get to it later.
  _N = command_line_value("N", 8);
  libmesh_assert_msg(_N > 0, "Invalid mesh size");

  switch (_dim)
    {
    case 1:
      MeshTools::Generation::build_line(_mesh, _N, 0.0, 1.0, EDGE2);
      break;
    case 2:
      MeshTools::Generation::build_square(_mesh, _N, _N, 0.0, 1.0, 0.0, 1.0, QUAD4);
      break;
    case 3:
      MeshTools::Generation::build_cube(_mesh, _N, _N, _N, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, HEX8);
      break;
    default:
      libmesh_assert_msg((_dim <= 3) && (_dim > 0), "Invalid mesh dimension");
      break;
    }

  // Determine the initial timestep size
  _dt0 = command_line_value("dt", 1.0/(10*_kappa*_N*_N*_N*_N));
  libmesh_assert_msg(_dt0>=0, "Negative initial timestep");

  _t0 = command_line_value("min_time", 0.0);
  _t1 = command_line_value("max_time", _t0 + 50.0*_dt0);
  libmesh_assert_msg(_t1 >= _t0, "Final time less than initial time");
  _T = _t1 - _t0;

  _cnWeight = command_line_value("crank_nicholson_weight", 1.0);
  libmesh_assert_msg(_cnWeight <= 1 && _cnWeight >= 0, "Crank-Nicholson weight must be between 0 and 1");

  // Initial state
  _initialState = STRIP;
  std::string initialState = command_line_value("initial_state", std::string("strip"));

  if (initialState == std::string("ball"))
    _initialState = BALL;
  else if (initialState == std::string("strip"))
    _initialState = STRIP;
  else if (initialState == std::string("rod"))
    _initialState = ROD;
  else
    libmesh_error_msg("Unknown initial state: neither ball nor rod nor srip");

  std::vector<Real> icenter;
  command_line_vector("initial_center", icenter);

  // Check that the point defining the center was in the right spatial dimension
  if (icenter.size() > _dim)
    libmesh_assert_msg(icenter.size() > _dim, "Invalid dimension for the initial state center of mass");

  // Pad
  icenter.resize(3);
  for (std::size_t i = icenter.size(); i < _dim; ++i)
    icenter[i] = 0.5;

  for (unsigned int i = _dim; i < 3; ++i)
    icenter[i] = 0.0;

  _initialCenter = Point(icenter[0], icenter[1], icenter[2]);
  _initialWidth = command_line_value("initial_width", 0.125);

  // Build the main equation encapsulated in the JR (Jacobian-Residual or J(R) "jet of R") object
  _jr = &(add_system<Biharmonic::JR>(std::string("Biharmonic::JR")));

  // Output options
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
  _exio.reset(new ExodusII_IO(_mesh));
  _o_dt = command_line_value("output_dt", 0.0);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
} // constructor



void Biharmonic::viewParameters()
{
  libMesh::out << "Biharmonic parameters:\n";

  // Print verbosity status
  if (_verbose)
    libMesh::out << "verbose mode is on\n";
  else
    libMesh::out << "verbose mode is off\n";

  // Print parameters
  libMesh::out << "mesh dimension           = " << _dim <<               "\n";
  libMesh::out << "initial linear mesh size = " << _N   <<               "\n";
  libMesh::out << "kappa                    = " << _kappa <<             "\n";
  libMesh::out << "growth                   = " << (int)_growth <<       "\n";
  libMesh::out << "degenerate               = " << (int)_degenerate <<   "\n";
  libMesh::out << "Cahn-Hillard             = " << (int)_cahn_hillard << "\n";
  libMesh::out << "netforce                 = " << (int)_netforce <<     "\n";
  libMesh::out << "energy                   = " << _energy        <<     "\n";
  libMesh::out << "tol                      = " << _tol           <<     "\n";
  libMesh::out << "theta                    = " << _theta         <<     "\n";
  libMesh::out << "theta_c                  = " << _theta_c       <<     "\n";
  libMesh::out << "log truncation           = " << _log_truncation <<    "\n";
  libMesh::out << "initial timestep size    = " << _dt0            <<    "\n";

  if (_initialState == STRIP)
    libMesh::out << "initial state:             strip\n";

  if (_initialState == ROD)
    libMesh::out << "initial state:             rod\n";

  if (_initialState == BALL)
    libMesh::out << "initial state:             ball\n";

  libMesh::out << "initial state center     = " << _initialCenter(0) << "\n";
  libMesh::out << "initial state width      = " << _initialWidth << "\n";
  libMesh::out << "initial time (min_time)  = " << _t0 << "\n";
  libMesh::out << "integration time         = " << _T  << "\n";
  libMesh::out << "final time   (max_time)  = " << _t1 << "\n";
  libMesh::out << "Crank-Nicholson weight   = " << _cnWeight << "\n";
  libMesh::out << "Output timestep          = " << _o_dt << "\n";
  libMesh::out << "Output filename base:      " <<  _ofile_base << "\n";
}




void Biharmonic::init()
{
  if (_verbose)
    libMesh::out << ">>> Initializing Biharmonic\n";

  _dt  =  0;
  _o_count = 0;
  EquationSystems::init();

  if (_verbose)
    libMesh::out << "<<< Initializing Biharmonic\n";
}





void Biharmonic::step(const Real & dt_in)
{
  // We need to update the old solution vector.
  // The old solution vector will be the current solution vector from the
  // previous time step. We use vector assignment.  Only TransientSystems
  // (and systems derived from them) contain old solutions.
  if (dt_in < 0)
    _dt = _dt0;
  else
    _dt = dt_in;

  *(_jr->old_local_solution) = *(_jr->current_local_solution);

  // this will localize the current solution, resulting in a
  // current_local_solution with correct ghost values
  _jr->solve();
}



void Biharmonic::output(int timestep,
                        const Real & t,
                        Real & o_t,
                        bool force)
{
#ifdef LIBMESH_HAVE_EXODUS_API
  if (!force && t - o_t < _o_dt)
    return;

  ++_o_count;

  if (_verbose)
    libMesh::out << "Writing state "
                 << timestep
                 << " at time "
                 << t
                 << " to file "
                 << _ofile
                 << "; output a total of "
                 << _o_count
                 << " states so far\n";

  _exio->write_timestep(_ofile, *this, timestep, t);

  if (!force)
    o_t = t;
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
}



void Biharmonic::run()
{
  Real t = _t0, o_t = 0.0;
  int timestep = 1;

  // Force-write the initial timestep
  output(timestep, t, o_t, true);

  while (t < _t1)
    {
      ++timestep;

      // A pretty update message
      if (_verbose)
        libMesh::out << "Solving for state " << timestep << ", time " << t << "\n";

      // Move biharmonic one timestep forward
      step();

      // Keep track of time and timestep
      t += _dt;

      // Output
      output(timestep, t, o_t);
    } // while(t < _t1)

  // Force-write the final timestep
  output(timestep, t, o_t, true);
}
