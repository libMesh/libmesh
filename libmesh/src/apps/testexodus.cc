 // C++ includes
#include <iostream>
#include <map>



// Local Includes
#include "mesh_config.h"
#include "libmesh.h"
#include "mesh.h"
#include "quadrature_gauss.h"
#include "quadrature_trap.h"
#include "quadrature_simpson.h"
#include "fe.h"
#include "dof_map.h"
#include "boundary_info.h"
#include "elem.h"
#include "point.h"
#include "general_system.h"
#include "equation_systems.h"
#include "perfmon.h"

#if defined(HAVE_PETSC)
extern "C" {
#include <petsc.h>
}
#elif defined(HAVE_MPI)
extern "C" {
#include <mpi.h>
}
#endif



void assemble_primary(EquationSystems& es,
		      const std::string& system_name);

void assemble_secondary(EquationSystems& es,
			const std::string& system_name);




int main (int argc, char** argv)
{
  libMesh::init (argc, argv);

  {
    if (argc < 4)
      std::cout << "Usage: ./prog -d DIM filename" << std::endl;
    else
      {
	std::cout << " Running ";

	for (int i=0; i<argc; i++)
	  std::cout << argv[i] << " ";

	std::cout << " on " << libMesh::n_processors()
		  << " processors, I am processor "  << libMesh::processor_id()
		  << std::endl
		  << std::endl;	    
      }
    
    // Variables to get us started
    const unsigned int dim = atoi(argv[2]);

    std::cout << std::endl
	      << "Processor = " << libMesh::processor_id()
	      << ", dim = " << dim
	      << std::endl
	      << std::endl;
    
    std::string basename  (argv[3]);
    std::string in_name    = basename;
    std::string bndry_name = basename;
    std::string tec_name   = basename;
    std::string gmv_name   = basename;

    in_name    += ".exd";
    bndry_name += "_bndry.gmv";
    tec_name   += ".plt";
    gmv_name   += ".gmv";
    
    PerfMon perfmon("Code performance");

    // declare a mesh...
    Mesh mesh(dim, libMesh::processor_id());
  
    // Read an Exodus mesh
    //
    // Then partition the domain and find all the neighbor
    // information
    {
      {
	PerfMon pm("Mesh Input performance");
	mesh.read_exd(in_name);
      }
      {
	PerfMon pm("Find_neighbors performance");
	mesh.find_neighbors();
      }
      {
	PerfMon pm("Partitioner performance");
	mesh.metis_partition(libMesh::n_processors());
      }
      mesh.print_info();
    }

    // Set up the equation system(s)
    EquationSystems es(mesh);

    /*
    // read the system from disk
    {
      es.read("out.xdr", Xdr::DECODE);
      mesh.write_gmv_binary(gmv_name,
			    es,
			    true);	
      
      return;
    }
    */

    {
      PerfMon pm("EquationSystem Initialization");

      // Set up the primary system
      {  
	es.add_system("primary");
	es("primary").add_variable("U", SECOND);
	es("primary").add_variable("V", SECOND);
      
	es("primary").get_dof_map().dof_coupling.resize(2);      
	es("primary").get_dof_map().dof_coupling(0,0) = 1;
	es("primary").get_dof_map().dof_coupling(1,1) = 1;
	
	es("primary").attach_assemble_function(assemble_primary);
      };
      
      // Set up the secondary system
      {  
	es.add_system("secondary");
	FEType fe_type(SECOND, MONOMIAL);
	es("secondary").add_variable("w", fe_type);
	
	es("secondary").attach_assemble_function(assemble_secondary);
      };
      
      es.set_parameter("linear solver tolerance") = 1.e-6;
      
      es.init();

      es.print_info();
    };
    
    // assemble & solve the primary 
    {
      PerfMon pm ("Solver Performance");
      
      // call the solver.
      es("primary").solve ();

      // clear the unneded matrix and
      // Petsc interface
      es("primary").matrix->clear ();
      es("primary").linear_solver_interface->clear ();
    };
    
    // assemble & solve the secondary
    {
      PerfMon pm ("Solver Performance");
      
      // call the solver.
      es("secondary").solve ();
    };
    

    
    mesh.write_gmv_binary(gmv_name,
			  es,
			  true);	
       

    es.write("out.xdr", Xdr::ENCODE);
    es.write("out.xda", Xdr::WRITE);
  };

  
  return libMesh::close();
};
  


void assemble_primary(EquationSystems& es,
		      const std::string& system_name)
{
  assert (system_name == "primary");
  
  const Mesh& mesh       = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  const int proc_id      = mesh.processor_id();
  
  // In this section we assemble the matrix and rhs
  PerfMon pm("Matrix Assembly (primary)");
	
  // Also use a 3x3x3 quadrature rule (3D).  Then tell the FE
  // about the geometry of the problem and the quadrature rule
  
  AutoPtr<FEBase> fe(FEBase::build(dim, es("primary").get_dof_map().variable_type(0)));
  QGauss qrule(dim, SEVENTH);
  //QTrap qrule(dim);
  
  fe->attach_quadrature_rule (&qrule);
  
  AutoPtr<FEBase> fe_face(FEBase::build(dim, es("primary").get_dof_map().variable_type(0)));
  QGauss   qface0(dim, FIFTH);
  //QTrap   qface0(dim);
  //QSimpson   qface0(dim);
  
  QGauss   qface1(dim-1, FIFTH);
  //QTrap   qface1(dim-1);
  //QSimpson   qface1(dim-1);
  
  fe_face->attach_quadrature_rule(&qface0);
  
  
  // These are references to cell-specific data
  const std::vector<Real>& JxW_face            = fe_face->get_JxW();
  const std::vector<Real>& JxW                 = fe->get_JxW();
  const std::vector<Point>& q_point            = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi   = fe->get_phi();
  const std::vector<std::vector<Point> >& dphi = fe->get_dphi();
  
  std::vector<unsigned int> dof_indices_U;
  std::vector<unsigned int> dof_indices_V;
  const DofMap& dof_map = es(system_name).get_dof_map();
  
  RealDenseMatrix   Kuu;
  RealDenseMatrix   Kvv;
  std::vector<Real> Fu;
  std::vector<Real> Fv;
  
  Real vol=0., area=0.;

  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {
      const Elem* elem = mesh.elem(e);

      // Find the next active element on my processor
      if (elem->processor_id() != proc_id) continue;
      if (!elem->active())                 continue;

      // recompute the element-specific data for the current element
      fe->reinit (elem);

      
      //fe->print_info();

      dof_map.dof_indices(elem, dof_indices_U, 0);
      dof_map.dof_indices(elem, dof_indices_V, 1);
      
      // zero the element matrix and vector
      Kuu.resize (phi.size(),
		  phi.size());
	    
      Kvv.resize (phi.size(),
		  phi.size());
	    
      Fu.resize (phi.size());
      Fv.resize (phi.size());
	    
      for (unsigned int i=0; i<phi.size(); i++)
	Fu[i] = Fv[i] = 0.;
	    
      // standard stuff...  like in code 1.
      for (unsigned int gp=0; gp<qrule.n_points(); gp++)
	{
	  for (unsigned int i=0; i<phi.size(); ++i)
	    {
	      // this is tricky.  ig is the _global_ dof index corresponding
	      // to the _global_ vertex number elem->node(i).  Note that
	      // in general these numbers will not be the same (except for
	      // the case of one unknown per node on one subdomain) so
	      // we need to go through the dof_map
		  
	      const Real f = q_point[gp]*q_point[gp];
	      //		    const Real f = (q_point[gp](0) +
	      //				    q_point[gp](1) +
	      //				    q_point[gp](2));
		    
	      // add jac*weight*f*phi to the RHS in position ig
		    
	      Fu[i] += JxW[gp]*f*phi[i][gp];
	      Fv[i] += JxW[gp]*f*phi[i][gp];
		    
	      for (unsigned int j=0; j<phi.size(); ++j)
		{
			
		  Kuu(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]));
			
		  Kvv(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]) +
				       1.*((dphi[i][gp])*(dphi[j][gp])));
		};
	    };
	  vol += JxW[gp];
	};


      // You can't compute "area" (perimeter) if you are in 2D
      if (dim == 3)
	{
	  for (unsigned int side=0; side<elem->n_sides(); side++)
	    if (elem->neighbor(side) == NULL)
	      {
		fe_face->reinit (&qface1, elem, side);
	    
		//fe_face->print_info();
	    
		for (unsigned int gp=0; gp<JxW_face.size(); gp++)
		  area += JxW_face[gp];
	      }
	}

      es("primary").rhs->add_vector(Fu,
				    dof_indices_U);
      es("primary").rhs->add_vector(Fv,
				    dof_indices_V);

      es("primary").matrix->add_matrix(Kuu,
				       dof_indices_U);
      es("primary").matrix->add_matrix(Kvv,
				       dof_indices_V);
    };
  
  std::cout << "Vol="  << vol << std::endl;

  if (dim == 3)
    std::cout << "Area=" << area << std::endl;
};
  


void assemble_secondary(EquationSystems& es,
			const std::string& system_name)
{
  assert (system_name == "secondary");
  
  const Mesh& mesh       = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  const int proc_id      = mesh.processor_id();
  
  // In this section we assemble the matrix and rhs
  PerfMon pm("Matrix Assembly (secondary)");

  AutoPtr<FEBase> fe(FEBase::build(dim, es("secondary").get_dof_map().variable_type(0)));
  QGauss qrule(dim, FIFTH);
  
  fe->attach_quadrature_rule (&qrule);
  
  // These are references to cell-specific data
  const std::vector<Real>& JxW                 = fe->get_JxW();
  const std::vector<Point>& q_point            = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi   = fe->get_phi();
  //const std::vector<std::vector<Point> >& dphi = fe->get_dphi();
  
  std::vector<unsigned int> dof_indices;
  const DofMap& dof_map = es(system_name).get_dof_map();
  
  RealDenseMatrix   Kww;
  std::vector<Real> Fw;
  
	
  for (unsigned int e=0; e<mesh.n_elem(); e++)
    {
      const Elem* elem = mesh.elem(e);

      // Find the next active element on my processor
      if (elem->processor_id() != proc_id) continue;
      if (!elem->active())                 continue;
	    
      dof_map.dof_indices(elem, dof_indices);
      
      // recompute the element-specific data for the current element
      fe->reinit (elem);
	    
      // zero the element matrix and vector
      Kww.resize (phi.size(),
		  phi.size());
	    
      Fw.resize (phi.size());
	    
      for (unsigned int i=0; i<phi.size(); i++)
	Fw[i] = 0.;
	    
      // standard stuff...  like in code 1.
      for (unsigned int gp=0; gp<qrule.n_points(); gp++)
	for (unsigned int i=0; i<phi.size(); ++i)
	  {
	    // this is tricky.  ig is the _global_ dof index corresponding
	    // to the _global_ vertex number elem->node(i).  Note that
	    // in general these numbers will not be the same (except for
	    // the case of one unknown per node on one subdomain) so
	    // we need to go through the dof_map
	    
	    const Real f = q_point[gp]*q_point[gp];
	    
	    Fw[i] += JxW[gp]*f*phi[i][gp];
	    
	    for (unsigned int j=0; j<phi.size(); ++j)
	      Kww(i,j) += JxW[gp]*(phi[i][gp])*(phi[j][gp]);
	  };

      es("secondary").matrix->add_matrix(Kww,
					 dof_indices);
      
      es("secondary").rhs->add_vector(Fw,
				      dof_indices);
      
    };
};
