 // C++ includes



// Local Includes
#include "libmesh_config.h"
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
#include "equation_systems.h"
#include "perfmon.h"
#include "steady_system.h"
#include "enum_xdr_mode.h"
#include "sparse_matrix.h"
#include "gmv_io.h"
#include "dense_matrix.h"
#include "dense_vector.h"

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
    Mesh mesh(dim);
  
    // Read an Exodus mesh
    //
    // Then partition the domain and find all the neighbor
    // information
    {
      {
	PerfMon pm("Mesh Input performance");
	mesh.read(in_name);
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



    // Set up the primary system
    SteadySystem& system1 = es.add_system<SteadySystem>("primary");
    system1.add_variable("U", SECOND);
    system1.add_variable("V", SECOND);
    system1.attach_assemble_function(assemble_primary);
      
    // Set up the secondary system
    SteadySystem& system2 = es.add_system<SteadySystem>("secondary");
    FEType fe_type(SECOND, MONOMIAL);
    system2.add_variable("w", fe_type);
    system2.attach_assemble_function(assemble_secondary);
      
    es.set_parameter("linear solver tolerance") = 1.e-6;
      
    es.init();

    es.print_info();
    
    // assemble & solve the primary
    system1.solve ();
      
    // assemble & solve the primary
    system2.solve ();

    // Write solution and mesh to file.
    GMVIO(mesh).write_equation_systems(gmv_name, es);	

    es.write("out.xdr", libMeshEnums::ENCODE);
    es.write("out.xda", libMeshEnums::WRITE);
  }

  
  return libMesh::close();
}
  


void assemble_primary(EquationSystems& es,
		      const std::string& system_name)
{
  assert (system_name == "primary");

  SteadySystem& system1 = es.get_system<SteadySystem>(system_name);
  
  const Mesh& mesh       = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  const int proc_id      = mesh.processor_id();
  
  // In this section we assemble the matrix and rhs
  PerfMon pm("Matrix Assembly (primary)");
	
  // Also use a 3x3x3 quadrature rule (3D).  Then tell the FE
  // about the geometry of the problem and the quadrature rule
  
  AutoPtr<FEBase> fe(FEBase::build(dim, system1.get_dof_map().variable_type(0)));
  QGauss qrule(dim, SEVENTH);
  //QTrap qrule(dim);
  
  fe->attach_quadrature_rule (&qrule);
  
  AutoPtr<FEBase> fe_face(FEBase::build(dim, system1.get_dof_map().variable_type(0)));
  
  QGauss   qface(dim-1, FIFTH);
  //QTrap   qface(dim-1);
  //QSimpson   qface(dim-1);
  
  fe_face->attach_quadrature_rule(&qface);
  
  
  // These are references to cell-specific data
  const std::vector<Real>& JxW_face            = fe_face->get_JxW();
  const std::vector<Real>& JxW                 = fe->get_JxW();
  const std::vector<Point>& q_point            = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi   = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  
  std::vector<unsigned int> dof_indices_U;
  std::vector<unsigned int> dof_indices_V;
  const DofMap& dof_map = system1.get_dof_map();
  
  DenseMatrix<Number>   Kuu;
  DenseMatrix<Number>   Kvv;
  DenseVector<Number> Fu;
  DenseVector<Number> Fv;
  
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
	    
      Fu.zero();
      Fv.zero();
      
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
		    
	      Fu(i) += JxW[gp]*f*phi[i][gp];
	      Fv(i) += JxW[gp]*f*phi[i][gp];
		    
	      for (unsigned int j=0; j<phi.size(); ++j)
		{
			
		  Kuu(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]));
			
		  Kvv(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]) +
				       1.*((dphi[i][gp])*(dphi[j][gp])));
		}
	    }
	  vol += JxW[gp];
	}


      // Compute surface area (perimiter in 2D)
      {
	for (unsigned int side=0; side<elem->n_sides(); side++)
	  if (elem->neighbor(side) == NULL)
	    {
	      fe_face->reinit (elem, side);
	      
	      //fe_face->print_info();
	      
	      for (unsigned int gp=0; gp<JxW_face.size(); gp++)
		area += JxW_face[gp];
	    }
      }

      system1.rhs->add_vector(Fu,
				    dof_indices_U);
      system1.rhs->add_vector(Fv,
				    dof_indices_V);

      system1.matrix->add_matrix(Kuu,
				       dof_indices_U);
      system1.matrix->add_matrix(Kvv,
				       dof_indices_V);
    }

    if (dim == 3)
      {
      	std::cout << "Vol="  << vol << std::endl;	
	std::cout << "Area=" << area << std::endl;
      }
    else if (dim == 2)
      {
      	std::cout << "Area="  << vol << std::endl;	
	std::cout << "Perimeter=" << area << std::endl;
      }
}
  


void assemble_secondary(EquationSystems& es,
			const std::string& system_name)
{
  assert (system_name == "secondary");

  SteadySystem& system2 = es.get_system<SteadySystem>(system_name);

  const Mesh& mesh       = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  const int proc_id      = mesh.processor_id();
  
  // In this section we assemble the matrix and rhs
  PerfMon pm("Matrix Assembly (secondary)");

  AutoPtr<FEBase> fe(FEBase::build(dim, system2.get_dof_map().variable_type(0)));
  QGauss qrule(dim, FIFTH);
  
  fe->attach_quadrature_rule (&qrule);
  
  // These are references to cell-specific data
  const std::vector<Real>& JxW                 = fe->get_JxW();
  const std::vector<Point>& q_point            = fe->get_xyz();
  const std::vector<std::vector<Real> >& phi   = fe->get_phi();
  //const std::vector<std::vector<Point> >& dphi = fe->get_dphi();
  
  std::vector<unsigned int> dof_indices;
  const DofMap& dof_map = es(system_name).get_dof_map();
  
  DenseMatrix<Number> Kww;
  DenseVector<Number> Fw;
  
	
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
	    
      Fw.zero();
      
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
	    
	    Fw(i) += JxW[gp]*f*phi[i][gp];
	    
	    for (unsigned int j=0; j<phi.size(); ++j)
	      Kww(i,j) += JxW[gp]*(phi[i][gp])*(phi[j][gp]);
	  }

      system2.matrix->add_matrix(Kww,
					 dof_indices);
      
      system2.rhs->add_vector(Fw,
				      dof_indices);
      
    }
}
