// $Id: ex6.C,v 1.21 2003-04-03 14:17:18 ddreyer Exp $

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


//#define _p(arg)  std::cout << "Here: " << #arg << std::endl
#define _p(arg) 



/**
 * C++ include files that we need
 */
#include <iostream>
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
 * Define the Finite and Infinite Element object.
 */
#include "fe.h"
#include "inf_fe.h"

/**
 * Define Gauss quadrature rules.
 */
#include "quadrature_gauss.h"

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
 * \mainpage Example 6
 *
 * \section Introduction
 *
 * WARNING! This example is under development.
 * Better do NOT use it!
 *
 * This is the sixth example program.  It builds on
 * the previous examples, and introduces the Infinite
 * Element class.  Note that the library must be compiled
 * with Infinite Elements enabled.  Otherwise, this
 * example will abort.
 */



/**
 * Function prototype.  This is similar to the Poisson
 * assemble function of example 4. 
 */
void assemble_wave(EquationSystems<GeneralSystem>& es,
		   const std::string& system_name);




int main (int argc, char** argv)
{

  /**
   * Initialize Petsc, like in example 2.
   */
  libMesh::init (argc, argv);


  /**
   * This short nice macro indicates the user
   * that the following code may be neither
   *  stable nor correct.
   */
  untested();

  /**
   * This example requires Infinite Elements
   */
#ifndef ENABLE_INFINITE_ELEMENTS

  std::cerr << "ERROR: This example requires the library to be " << std::endl
	    << " compiled with Infinite Element support!" << std::endl;
  here();

  return 0;

#else


  /**
   * This example is designed for real numbers only.
   */
# ifdef USE_COMPLEX_NUMBERS

  std::cerr << "ERROR: This example is not intended for " << std::endl
	    << " use with complex numbers." << std::endl;
  here();

  return 0;

# endif


  /**
   * Braces are used to force object scope, like in example 2
   */   
  {    
    /**
     * For the moment, only allow 3D
     */
    const unsigned int dim = 3; 

    /**
     * Tell the user what we are doing.
     */
    std::cout << "Running ex6 with dim = " << dim << std::endl << std::endl;        
    
    /**
     * Create a mesh with user-defined dimension 
     */
    Mesh mesh (dim);

    /**
     * Use the internal mesh generator to create 8 elements
     * on the square [-1,1]^D, of type \p Quad9 or \p Hex27.
     */
    mesh.build_cube (2, 2, 2,
//    mesh.build_cube (1, 1, 1,
		     -1., 1.,
		     -1., 1.,
		     -1., 1.,
// 		     (dim == 2) ? QUAD8 : HEX27);
		     (dim == 2) ? QUAD4 : HEX8); //HEX20); //HEX27);

    /**
     * Print information about the mesh to the screen.
     */
    mesh.print_info();

    /**
     * Build infinite elements on the outer boundary of the existing 
     * volume mesh.  This method automatically determines the
     * origin of the infinite elements.  The \p bool determines
     * whether to be verbose.
     */
    mesh.build_inf_elem(true);

    /**
     * Print information about the mesh to the screen.
     */
    mesh.print_info();

    /**
     * Save only the mesh, with infinite elements added.
     */
    mesh.write_gmv ("ifems_added.gmv");

    /**
     * After building finite elements, we have to let 
     * the elements find their neighbors.
     */
//TODO:[DD] Does this work with infinite elements?
    mesh.find_neighbors();
    


    /**
     * Create an equation systems object.
     */
    EquationSystems<GeneralSystem> equation_systems (mesh);

    /**
     * Declare the system and its variables.
     */
    {
      /**
       * Create a system named "Wave"
       */
      equation_systems.add_system("Wave");
      
      /**
       * Create an FEType describing the approximation
       * characteristics of the InfFE object.  Note that
       * the constructor automatically defaults to some
       * sensible values.  But use \p FIRST order 
       * approximation.
       */
      FEType fe_type(FIRST);

      /**
       * Add the variable "p" to "Wave".  Note that there exist
       * various approaches in adding variables.  In example 3, 
       * \p add_variable took the order of approximation and used
       * default values for the \p FEFamily, while here the \p FEType 
       * is used.
       */
      equation_systems("Wave").add_variable("p", fe_type);

      /**
       * Give the system a pointer to the matrix assembly
       * function.
       */
      equation_systems("Wave").attach_assemble_function (assemble_wave);
      

      _p(kurz_vor_eqnsys_init);
      /**
       * Initialize the data structures for the equation system.
       */
      equation_systems.init();
      
      /**
       * Prints information about the system to the screen.
       */
      equation_systems.print_info();
    }



    /**
     * Solve the system "Wave".
     */
    equation_systems("Wave").solve();


    _p(Write results);
    /**
     * After solving the system write the solution
     * to a GMV-formatted plot file.
     */
    mesh.write_gmv ("out.gmv", equation_systems);
  }

  
  /**
   * All done.  
   */
  return libMesh::close ();


#endif // else part of ifndef ENABLE_INFINITE_ELEMENTS
}




void assemble_wave(EquationSystems<GeneralSystem>& es,
		   const std::string& system_name)
{
  _p(starting assemble_wave);

  /**
   * It is a good idea to make sure we are assembling
   * the proper system.
   */
  assert (system_name == "Wave");


#ifdef ENABLE_INFINITE_ELEMENTS


  /**
   * Get a constant reference to the mesh object.
   */
  const Mesh& mesh = es.get_mesh();

  /**
   * The dimension that we are running
   */
  const unsigned int dim = mesh.mesh_dimension();

  /**
   * Get a constant reference to the Finite Element type
   * for the first (and only) variable in the system.
   */
  FEType fe_type = es("Wave").get_dof_map().variable_type(0);

  /**
   * Build a Finite Element object of the specified type.  Since the
   * \p FEBase::build() member dynamically creates memory we will
   * store the object as an \p AutoPtr<FEBase>.  This can be thought
   * of as a pointer that will clean up after itself.
   */
  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

  /**
   * Do the same for an infinite element.
   */
  AutoPtr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));
  


  /**
   * A 2nd order Gauss quadrature rule for numerical integration.
   */
  QGauss qrule (dim, SECOND);

  /**
   * Tell the finite element object to use our quadrature rule.
   */
  fe->attach_quadrature_rule (&qrule);
  _p(fe->attach just finished.);


  /**
   * Due to its internal structure, the infinite element handles 
   * quadrature rules differently.  It takes the quadrature
   * rule which has been initialized for the FE object, but
   * creates suitable quadrature rules by itself.  The user
   * need not worry about this.
   */
  inf_fe->attach_quadrature_rule (&qrule);
  _p(inf_fe->attach just finished.);





  /**
   * A reference to the \p DofMap object for this system.  The \p DofMap
   * object handles the index translation from node and element numbers
   * to degree of freedom numbers.  We will talk more about the \p DofMap
   * in future examples.
   */
  const DofMap& dof_map = es("Wave").get_dof_map();

  /**
   * Define data structures to contain the element matrix
   * and right-hand-side vector contribution.  Following
   * basic finite element terminology we will denote these
   * "Ke",  "Ce", "Me", and "Fe" for the stiffness, damping
   * and mass matrices, and the load vector.  Note that in Acoustics,
   * these descriptors do not match the true physical meaning
   * of the projectors.  The final overall system, however, 
   * resembles the conventional notation, again.
   */
  DenseMatrix<Number> Ke;
  DenseMatrix<Number> Ce;
  DenseMatrix<Number> Me;

  DenseVector<Number> Fe;

  /**
   * This vector will hold the degree of freedom indices for
   * the element.  These define where in the global system
   * the element degrees of freedom get mapped.
   */
  std::vector<unsigned int> dof_indices;



  _p(Start element loop);
  
  /**
   *--------------------------------------------------------------------
   * Now we will loop over all the elements in the mesh.
   * We will compute the element matrix and right-hand-side
   * contribution.
   */
  
  const_elem_iterator           el (mesh.elements_begin());
  const const_elem_iterator end_el (mesh.elements_end());
  
  for ( ; el != end_el; ++el)
    {
      /**
       * Store a pointer to the element we are currently
       * working on.  This allows for nicer syntax later.
       */
      const Elem* elem = *el;

      /**
       * Get the degree of freedom indices for the
       * current element.  These define where in the global
       * matrix and right-hand-side this element will
       * contribute to.
       */
      dof_map.dof_indices (elem, dof_indices);


      /**
       * The mesh contains both finite and infinite elements.  These
       * elements are handled through different classes, namely
       * \p FE and \p InfFE, respectively.  However, since both
       * are derived from \p FEBase, they share the same interface,
       * and overall burden of coding is greatly reduced through
       * using a pointer, which is adjusted appropriately to the
       * current element type.
       */
      FEBase* cfe=NULL;


      /**
       *----------------------------------------------------------
       * This here is almost the only place where we need to
       * distinguish between finite and infinite elements.
       * For faster computation, however, different approaches
       * may be feasible.
       *
       * Up to now, we do not know what kind of element we
       * have.  Aske the element of what type it is:
        */
      if (elem->infinite())
        {
	  /** 
	   * We have an infinite element.  Let \p cfe point
	   * to our \p InfFE object.  This is handled through
	   * an AutoPtr.  Through the \p AutoPtr::get() we "borrow"
	   * the pointer, while the \p  AutoPtr \p inf_fe is
	   * still in charge of memory management.
	   */
	  cfe = inf_fe.get(); 
	  _p(InfFE);
	}
      else
        {
	  /** 
	   * This is a conventional finite element.  Let \p fe handle it.
	   */
  	  cfe = fe.get();
	  _p(FE);


	  /**
	   *----------------------------------------------------------------
	   * Boundary conditions.  Currently, only  natural boundary
	   * conditions on faces of finite elements are supported (this
	   * is likely to change in the future!).  Therefore, we do this 
	   * only for finite elements. 
	   */
	    {
	      /**
	       * Zero the RHS for this element. 
	       */
	      Fe.resize (dof_indices.size());


	      /**
	       * See previous example(s) for details
	       */
	      for (unsigned int side=0; side<elem->n_sides(); side++)
		  if (elem->neighbor(side) == NULL)
		    {
		      AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
	      
		      QGauss qface (dim, SECOND);
	      
		      fe_face->attach_quadrature_rule (&qface);
	      
		      /**
		       * The value of the shape functions at the quadrature
		       * points.
		       */
		      const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
	      
		      /**
		       * The Jacobian * Quadrature Weight at the quadrature
		       * points on the face.
		       */
		      const std::vector<Real>& JxW_face = fe_face->get_JxW();
	      
		      /**
		       * The XYZ locations (in physical space) of the
		       * quadrature points on the face.  This is where
		       * we will interpolate the boundary value function.
		       */
		      //const std::vector<Point >& qface_point = fe_face->get_xyz();
	      
		      /**
		       * Compute the shape function values on the element
		       * face.
		       */
		      fe_face->reinit(elem, side);

	      
		      /**
		       * Loop over the face quagrature points for integration.
		       */
		      for (unsigned int qp=0; qp<fe_face->n_quadrature_points(); qp++)
		        {
			  /**
			   * The location on the boundary of the current
			   * face quadrature point.
			   */
			  //const Real xf = qface_point[qp](0);
			  //const Real yf = qface_point[qp](1);
			  //const Real zf = qface_point[qp](2);
		  
			  /**
			   * The boundary value.
			   */
			  const Real value = 1.;
		  		  
			  /**
			   * Right-hand-side contribution.
			   */
			  for (unsigned int i=0; i<fe_face->n_shape_functions(); i++)
			    {
			      Fe(i) += JxW_face[qp]*value*phi_face[i][qp];
			    }
		  
			} // end face quadrature point loop	  


		      es("Wave").rhs->add_vector    (Fe, dof_indices);

		    } // end if (elem->neighbor(side) == NULL)

	    } // end boundary condition section	     

	} // else ( if (elem->infinite())) )



      /**
       *----------------------------------------------------------
       * Now this is all independent of whether we use an \p FE
       * or an \p InfFE.  Nice, hm? ;-)
       *
       * Compute the element-specific data, as described
       * in previous examples.
       */
      cfe->reinit (elem);

      /**
       * This is slightly different from the Poisson solver:
       * Since the finite element object may change, we have to
       * initialize the constant references to the data fields
       * each time again, when a new element is processed.
       *
       * The element Jacobian * quadrature weight at each integration point.   
       */
      const std::vector<Real>& JxW = cfe->get_JxW();

      /**
       * The element shape functions evaluated at the quadrature points.
       */
      const std::vector<std::vector<Real> >& phi = cfe->get_phi();

      /**
       * The element shape function gradients evaluated at the quadrature
       * points.
       */
      const std::vector<std::vector<Point> >& dphi = cfe->get_dphi();

      /**
       * The infinite elements need more data fields than conventional FE.  
       * These are the gradients of the phase term \p dphase, an additional 
       * radial weight for the test functions \p Sobolev_weight, and its
       * gradient.
       * 
       * Note that these data fields are also initialized appropriately by
       * the \p FE method, so that the weak form (below) is valid for both
       * finite and infinite elements.
       */
      const std::vector<Point>& dphase  = cfe->get_dphase();
      const std::vector<Real>&  weight  = cfe->get_Sobolev_weight();
      const std::vector<Point>& dweight = cfe->get_Sobolev_dweight();


      /**
       * Zero the element matrices.  Boundary conditions were already
       * processed in the \p FE-only section, see above.
       */
      Ke.resize (dof_indices.size(), dof_indices.size());
      Ce.resize (dof_indices.size(), dof_indices.size());
      Me.resize (dof_indices.size(), dof_indices.size());

      /**
       * The total number of quadrature points for infinite elements
       * @e has to be determined in a different way, compared to
       * conventional finite elements.  This type of access is also
       * valid for finite elements, so this can safely be used
       * anytime, instead of asking the quadrature rule, as
       * seen in previous examples.
       */
      unsigned int max_qp = cfe->n_quadrature_points();


      /**
       * Loop over the quadrature points. 
       */
      for (unsigned int qp=0; qp<max_qp; qp++)
        {
	  /**
	   * Similar to the modified access to the number of quadrature 
	   * points, the number of shape functions may also be obtained
	   * in a different manner.  This offers the great advantage
	   * of being valid for both finite and infinite elements.
	   */
	  unsigned int n_sf = cfe->n_shape_functions();

	  /**
	   * Now we will build the element matrices.  Since the infinite
	   * elements are based on a Petrov-Galerkin scheme, the
	   * resulting system matrices are non-symmetric. The additional
	   * weight, described before, is part of the trial space.
	   *
	   * For the finite elements, though, these matrices are symmetric
	   * just as we know them, since the additional fields \p dphase,
	   * \p weight, and \p dweight are initialized appropriately.
	   *
	   * test functions:    weight[qp]*phi[i][qp]
	   * trial functions:   phi[j][qp]
	   * phase term:        phase[qp]
	   * 
	   * derivatives are similar, but note that these are of type
	   * Point, not of type Real.
	   */
	  for (unsigned int i=0; i<n_sf; i++)
	    for (unsigned int j=0; j<n_sf; j++)
	      {
		//         (ndt*Ht + nHt*d) * nH 
		Ke(i,j) +=
		    (                                /*    (                         */
			(                            /*      (                       */
			  dweight[qp] * phi[i][qp]   /*        Point * Real  = Point */
			  +                          /*        +                     */
			  dphi[i][qp] * weight[qp]   /*        Point * Real  = Point */
			) * dphi[j][qp]              /*      )       * Point = Real  */
		    ) * JxW[qp];                     /*    )         * Real  = Real  */


		// (d*Ht*nmut*nH - ndt*nmu*Ht*H - d*nHt*nmu*H)
		Ce(i,j) +=
		    (                                /*    (                         */
			(dphase[qp] * dphi[j][qp])   /*      (Point * Point) = Real  */
			* weight[qp] * phi[i][qp]    /*      * Real * Real   = Real  */
			-                            /*      -                       */
			(dweight[qp] * dphase[qp])   /*      (Point * Point) = Real  */
			* phi[i][qp] * phi[j][qp]    /*      * Real * Real   = Real  */
			-                            /*      -                       */
			(dphi[i][qp] * dphase[qp])   /*      (Point * Point) = Real  */
			* weight[qp] * phi[j][qp]    /*      * Real * Real   = Real  */
		    ) * JxW[qp];                     /*    )         * Real  = Real  */


		// (d*Ht*H * (1 - nmut*nmu))
		Me(i,j) +=
		    (                                          /*    (                                  */
			(1. - (dphase[qp] * dphase[qp]))       /*      (Real  - (Point * Point)) = Real */
			* phi[i][qp] * phi[j][qp] * weight[qp] /*      * Real *  Real  * Real    = Real */
		    ) * JxW[qp];                               /*    ) * Real                    = Real */

	      } // end of the matrix summation loop

	} // end of quadrature point loop

          
      /**
       *----------------------------------------------------------------
       * The element matrices are now built for this element.  
       * Collect them in Ke, and then add them to the global matrix.  
       * The \p PetscMatrix::add_matrix() member does this for us.
       */
      //Ke.print();

// Doesn't make sense!
      Ke.add(1., Ce);
      Ke.add(1., Me);
      es("Wave").matrix->add_matrix (Ke, dof_indices);

      
    } // end of element loop

  _p(Finished element loop);

  
#else

  /* dummy assert */
  assert(es.get_mesh().mesh_dimension() != 1);

#endif //ifdef ENABLE_INFINITE_ELEMENTS


  /**
   * All done!
   */
  return;
}

