// $Id: ex5.C,v 1.16 2003-03-03 02:15:57 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2003  Benjamin S. Kirk
  
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
  
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA




/**
 * C++ include files that we need
 */
#include <iostream>
#include <sstream> 
#include <algorithm>
#include <math.h>

/**
 * Basic include file needed for the mesh functionality.
 */
#include "libmesh.h"
#include "mesh.h"
#include "general_system.h"
#include "equation_systems.h"

/**
 * Define the Finite Element object.
 */
#include "fe.h"

/**
 * Define the base quadrature class, with which
 * specialized quadrature rules will be built.
 */
#include "quadrature_gauss.h"

/**
 * Include the namespace \p QuadratureRules for
 * some handy descriptions.
 */
#include "quadrature_rules.h"

/**
 * Define useful datatypes for finite element
 * matrix and vector components.
 */
#include "dense_matrix.h"
#include "dense_vector.h"

/**
 * Define the DofMap, which handles degree of freedom
 * indexing.
 */
#include "dof_map.h"






/**
 * \mainpage Example 5
 *
 * \section Introduction
 *
 * This is the fifth example program.  It builds on
 * the previous two examples, and extends the use
 * of the \p AutoPtr<> as a convenient build method to
 * determine the quadrature rule at run time.
 */


/**
 * Function prototype, as before.
 */
void assemble_poisson(EquationSystems& es,
                      const std::string& system_name);



/**
 * Exact solution function prototype, as before.
 */
Real exact_solution (const Real x,
		     const Real y,
		     const Real z = 0.);


/**
 * The quadrature type the user requests.
 */
QuadratureType quad_type=INVALID_Q_RULE;




int main (int argc, char** argv)
{
  /**
   * Initialize Petsc, like in example 2.
   */
  libMesh::init (argc, argv);
  
  /**
   * Braces are used to force object scope, like in example 2
   */   
  {
    /**
     * Check for proper usage.  The quadrature rule
     * must be given at run time.
     */
    if (argc < 3)
      {
	std::cerr << "Usage: " << argv[0] << " -q n"
		  << std::endl;
	std::cerr << "  where n stands for:" << std::endl;

	/**
	 * Note that only some of all quadrature rules are
	 * valid choices.  For example, the Jacobi quadrature
	 * is actually a "helper" for higher-order rules,
	 * included in QGauss.
	 */
	for (unsigned int n=0; n<QuadratureRules::num_valid_elem_rules; n++)
	  std::cerr << "  " << QuadratureRules::valid_elem_rules[n] << "    " 
		    << QuadratureRules::name(QuadratureRules::valid_elem_rules[n])
		    << std::endl;

	std::cerr << std::endl;
	
	error();
      }
    
    /**
     * Tell the user what we are doing.
     */
    else 
      {
	std::cout << "Running " << argv[0];
	
	for (int i=1; i<argc; i++)
	  std::cout << " " << argv[i];
	
	std::cout << std::endl << std::endl;
      };
    

    /**
     * Set the quadrature rule type that the user wants from argv[2]
     */
    quad_type = static_cast<QuadratureType>(atoi(argv[2]));


    /**
     * Independence of dimension has already been shown in
     * example 4.  For the time being, restrict to 3 dimensions.
     */
    const unsigned int dim=3;
    
    /**
     * The following is identical to example 4, and therefore
     * not commented.  Differences are mentioned when present.
     */
    Mesh mesh (dim);
    
    mesh.build_cube (16, 16, 16,
		     -1., 1.,
		     -1., 1.,
		     -1., 1.,
		     HEX8);
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    {
      equation_systems.add_system("Poisson");
      
      equation_systems("Poisson").add_variable("u", FIRST);

      equation_systems("Poisson").attach_assemble_function (assemble_poisson);
      
      equation_systems.init();
      
      equation_systems.print_info();
    };


    equation_systems("Poisson").solve();

    /**
     * "Personalize" the output, with the
     * number of the quadrature rule appended.
     */
    std::ostringstream f_name;
    f_name << "out_" << quad_type << ".gmv";

    mesh.write_gmv (f_name.str(), equation_systems);
  };


  /**
   * All done.
   */
  return libMesh::close ();
};




void assemble_poisson(EquationSystems& es,
                      const std::string& system_name)
{
  assert (system_name == "Poisson");

  const Mesh& mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  FEType fe_type = es("Poisson").get_dof_map().variable_type(0);

  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  Below, the
   * functionality of \p AutoPtr's is described more detailed in 
   * the context of building quadrature rules.
   */
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
  
  /**
   * Now this deviates from example 4.  we create a 
   * 5th order quadrature rule of user-specified type
   * for numerical integration.  Note that not all
   * quadrature rules support this order.
   */
  AutoPtr<QBase> qrule(QBase::build(quad_type, dim, THIRD));


  /**
   * Tell the finte element object to use our
   * quadrature rule.  Note that a \p AutoPtr<QBase> returns
   * a \p QBase* pointer to the object it handles with \p get().  
   * However, using \p get(), the \p AutoPtr<QBase> \p qface0 is 
   * still in charge of this pointer. I.e., when \p qface0 goes 
   * out of scope, it will safely delete the \p QBase object it 
   * points to.  This behavior may be overridden using
   * \p AutoPtr<Xyz>::release(), but is currently not
   * recommended.
   */
  fe->attach_quadrature_rule (qrule.get());

  

  /**
   *--------------------------------------------------------------------
   * This is again identical to example 4, and not commented.
   */

  const std::vector<Real>& JxW = fe->get_JxW();

  const std::vector<Point>& q_point = fe->get_xyz();

  const std::vector<std::vector<Real> >& phi = fe->get_phi();

  const std::vector<std::vector<Point> >& dphi = fe->get_dphi();

  const DofMap& dof_map = es("Poisson").get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<unsigned int> dof_indices;




  /**
   *--------------------------------------------------------------------
   * Now we will loop over all the elements in the mesh.
   * See example 3 for details.
   */
  
  const_elem_iterator           el (mesh.elements_begin());
  const const_elem_iterator end_el (mesh.elements_end());
  
  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);

      fe->reinit (elem);

      Ke.resize (dof_indices.size(),
		 dof_indices.size());

      Fe.resize (dof_indices.size());



      /**
       *----------------------------------------------------------------
       * Now loop over the quadrature points.  This handles
       * the numeric integration.  Note the slightly different
       * access to the QBase members!
       */
      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
	{

	  for (unsigned int i=0; i<phi.size(); i++)
	    for (unsigned int j=0; j<phi.size(); j++)
	      {
		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
	      };


	  for (unsigned int i=0; i<phi.size(); i++)
	    {
	      const Real x = q_point[qp](0);
	      const Real y = q_point[qp](1);
	      const Real z = q_point[qp](2);
	      const Real eps = 1.e-3;
	       
	      const Real uxx = (exact_solution(x-eps,y,z) +
				exact_solution(x+eps,y,z) +
				-2.*exact_solution(x,y,z))/eps/eps;
	      
	      const Real uyy = (exact_solution(x,y-eps,z) +
				exact_solution(x,y+eps,z) +
				-2.*exact_solution(x,y,z))/eps/eps;
	      
	      const Real uzz = (exact_solution(x,y,z-eps) +
				exact_solution(x,y,z+eps) +
				-2.*exact_solution(x,y,z))/eps/eps;

	      const Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
	      
	      Fe(i) += JxW[qp]*fxy*phi[i][qp];
	    }; // end of the RHS summation loop
	  
	}; // end of quadrature point loop




      
      /**
       *----------------------------------------------------------------
       * Most of this has already been seen before, except
       * for the build routines of QBase, described below
       */
      {
	for (unsigned int side=0; side<elem->n_sides(); side++)
	  if (elem->neighbor(side) == NULL)
	    {
	      AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
	      /**
	       * As already seen in example 3, boundary integration 
	       * requires a quadraure rule.  Here, however,
	       * we use the more convenient way of building this
	       * rule.  Note that one could also have initialized
	       * the face quadrature rules with the type directly
	       * determined from \p qrule, namely through:
	       * \verbatim
	         AutoPtr<QBase>  qface1 (QBase::build(qrule->type(),
		                                      dim-1, 
						      THIRD));
		 \endverbatim
	       * And again: using the \p AutoPtr<QBase> relaxes
	       * the need to delete the object afterwards,
	       * they clean up themselves.
	       */
	      AutoPtr<QBase>  qface (QBase::build(quad_type,
						  dim-1, 
						  THIRD));
	      
	      /**
	       * Tell the finte element object to use our
	       * quadrature rule.  Note that a \p AutoPtr<QBase> returns
	       * a \p QBase* pointer to the object it handles with \p get().  
	       * However, using \p get(), the \p AutoPtr<QBase> \p qface is 
	       * still in charge of this pointer. I.e., when \p qface goes 
	       * out of scope, it will safely delete the \p QBase object it 
	       * points to.  This behavior may be overridden using
	       * \p AutoPtr<Xyz>::release(), but is not recommended.
	       */
	      fe_face->attach_quadrature_rule (qface.get());
	      
	      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();

	      const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      
	      const std::vector<Point >& qface_point = fe_face->get_xyz();
	      
	      /**
	       * Compute the shape function values on the element
	       * face.
	       */
	      fe_face->reinit(elem, side);
	      
	      /**
	       * Loop over the face quagrature points for integration.
	       * Note that the \p AutoPtr<QBase> overloaded the operator->,
	       * so that QBase methods may safely be accessed.  It may
	       * be said: accessing an \p AutoPtr<Xyz> through the
	       * "." operator returns \p AutoPtr methods, while access
	       * through the "->" operator returns Xyz methods.
	       * This allows almost no change in syntax when switching
	       * to "safe pointers".
	       */
	      for (unsigned int qp=0; qp<qface->n_points(); qp++)
		{
		  const Real xf = qface_point[qp](0);
		  const Real yf = qface_point[qp](1);
		  const Real zf = qface_point[qp](2);
		  
		  const Real penalty = 1.e10;
		  
		  const Real value = exact_solution(xf, yf, zf);
		  
		  for (unsigned int i=0; i<phi_face.size(); i++)
		    for (unsigned int j=0; j<phi_face.size(); j++)
		      {
			Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
		      };

		  for (unsigned int i=0; i<phi_face.size(); i++)
		    {
		      Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
		    };
		  
		}; // end face quadrature point loop	  
	    }; // end if (elem->neighbor(side) == NULL)
      }; // end boundary condition section	  


      
      
      /**
       *----------------------------------------------------------------
       * The element matrix and right-hand-side are now built
       * for this element.  Add them to the global matrix and
       * right-hand-side vector.  The \p PetscMatrix::add_matrix()
       * and \p PetscVector::add_vector() members do this for us.
       */
      
      es("Poisson").matrix->add_matrix (Ke, dof_indices);
      es("Poisson").rhs->add_vector    (Fe, dof_indices);
      
    }; // end of element loop


  
  /**
   * All done!
   */
  return;
};
