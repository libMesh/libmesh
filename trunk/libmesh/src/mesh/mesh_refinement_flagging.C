// $Id: mesh_refinement_flagging.C,v 1.12 2004-11-08 00:11:05 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "libmesh_config.h"

// only compile these functions if the user requests AMR support
#ifdef ENABLE_AMR

// C++ includes
#include <algorithm> // for std::sort

// Local includes
#include "error_vector.h"
#include "mesh_refinement.h"
#include "mesh_base.h"
#include "elem.h"



//-----------------------------------------------------------------
// Mesh refinement methods
void MeshRefinement::flag_elements_by_error_fraction (const ErrorVector& error_per_cell_in,
						      const Real refine_fraction,
						      const Real coarsen_fraction,
						      const unsigned int max_level)
{
  // Copy the input error_per_cell so that we can modify it
  ErrorVector error_per_cell (error_per_cell_in);
  
  // Check for valid fractions..
  // The fraction values must be in [0,1]
  assert (refine_fraction  >= 0.);  assert (refine_fraction  <= 1.);
  assert (coarsen_fraction >= 0.);  assert (coarsen_fraction <= 1.);

  // Clean up the refinement flags.  These could be left
  // over from previous refinement steps.
  this->clean_refinement_flags();
  

//   // Set the error_per_cell to zero for elements at the
//   // maximum allowable refinement level
//   {
//     active_elem_iterator       elem_it (_mesh.elements_begin());
//     const active_elem_iterator elem_end(_mesh.elements_end());

//     for (; elem_it != elem_end; ++elem_it)
//       {
// 	Elem* elem               = *elem_it;
// 	const unsigned int id    = elem->id();
// 	const unsigned int level = elem->level();

// 	assert (id < error_per_cell.size());
      
// 	if (level >= max_level)
// 	  error_per_cell[id] = 0.;
//       }
//   }
  
  // Get the minimum, maximum, and delta error values
  // for the elements
  const Real error_max   = error_per_cell.maximum();
  const Real error_min   = error_per_cell.minimum();
  const Real error_delta = (error_max - error_min);


  // Compute the cutoff values for coarsening and refinement
  const Real refine_cutoff  = (1.- refine_fraction)*error_max;
  const Real coarsen_cutoff = coarsen_fraction*error_delta + error_min;

//   // Print information about the error
//   std::cout << " Error Information:"                     << std::endl
// 	    << " ------------------"                     << std::endl
// 	    << "   min:              " << error_min      << std::endl
// 	    << "   max:              " << error_max      << std::endl
// 	    << "   delta:            " << error_delta    << std::endl
// 	    << "     refine_cutoff:  " << refine_cutoff  << std::endl
// 	    << "     coarsen_cutoff: " << coarsen_cutoff << std::endl;
  
  

  // Loop over the elements and flag them for coarsening or
  // refinement based on the element error
//   active_elem_iterator       elem_it (_mesh.elements_begin());
//   const active_elem_iterator elem_end(_mesh.elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem             = *elem_it;
      const unsigned int id  = elem->id();

      assert (id < error_per_cell.size());
      
      const float elem_error = error_per_cell[id];

      // Flag the element for coarsening if its error
      // is <= coarsen_fraction*delta + error_min
      if (elem_error <= coarsen_cutoff)
	{
	  elem->set_refinement_flag(Elem::COARSEN);
	}
      
      // Flag the element for refinement if its error
      // is >= refinement_cutoff.
      if (elem_error >= refine_cutoff)
	if (elem->level() < max_level)
	  elem->set_refinement_flag(Elem::REFINE);
    }
}



void MeshRefinement::flag_elements_by_elem_fraction (const ErrorVector& error_per_cell_in,
						     const Real refine_fraction,
						     const Real coarsen_fraction,
						     const unsigned int max_level)
{
  // Copy the input error_per_cell so that we can modify it
  ErrorVector error_per_cell (error_per_cell_in);
  
  // Check for valid fractions..
  // The fraction values must be in [0,1]
  assert (refine_fraction  >= 0.);  assert (refine_fraction  <= 1.);
  assert (coarsen_fraction >= 0.);  assert (coarsen_fraction <= 1.);

  // The number of active elements in the mesh
  const unsigned int n_active_elem  = _mesh.n_elem();

  // The number of elements to flag for coarsening
  const unsigned int n_elem_coarsen = static_cast<unsigned int>(coarsen_fraction * n_active_elem);

  // The number of elements to flag for refinement
  const unsigned int n_elem_refine  = static_cast<unsigned int>(refine_fraction  * n_active_elem);


  
  // Clean up the refinement flags.  These could be left
  // over from previous refinement steps.
  this->clean_refinement_flags();



//   // Set the error_per_cell to zero for elements at the
//   // maximum allowable refinement level
//   {
//     active_elem_iterator       elem_it (_mesh.elements_begin());
//     const active_elem_iterator elem_end(_mesh.elements_end());

//     for (; elem_it != elem_end; ++elem_it)
//       {
// 	Elem* elem               = *elem_it;
// 	const unsigned int id    = elem->id();
// 	const unsigned int level = elem->level();

// 	assert (id < error_per_cell.size());
      
// 	if (level >= max_level)
// 	  error_per_cell[id] = 0.;
//       }
//   }
  
  // This vector stores the error and element number for all the
  // active elements.  It will be sorted and the top & bottom
  // elements will then be flagged for coarsening & refinement
  std::vector<float> sorted_error;

  sorted_error.reserve (n_active_elem);

  // Loop over the active elements and create the entry
  // in the sorted_error vector
//   const_active_elem_iterator       elem_it (_mesh.const_elements_begin());
//   const const_active_elem_iterator elem_end(_mesh.const_elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;
      const unsigned int elem_number = elem->id();
      const float        elem_error  = error_per_cell[elem_number];

      sorted_error.push_back (elem_error);
    }

  // Now sort the sorted_error vector
  std::sort (sorted_error.begin(), sorted_error.end());
  

  float
    top_error= 0., bottom_error = 0.,
    min_error= 0., max_error    = 0.;

  // Get the maximum error value corresponding to the
  // bottom n_elem_coarsen elements
  {
    std::vector<float>::iterator
      it  = sorted_error.begin();

    min_error = *it;
    
    std::advance (it, n_elem_coarsen);

    bottom_error = *it;
  }

  
  // Get the minimum error value corresponding to the
  // top n_elem_refine elements
  {
    std::vector<float>::reverse_iterator
      it  = sorted_error.rbegin();

    max_error = *it;
    
    std::advance (it, n_elem_refine);

    top_error = *it;
  }

  assert (max_error >= min_error);
  const Real delta       = max_error - min_error;  
  const Real top_frac    = (max_error    - top_error)/max_error;
  const Real bottom_frac = (bottom_error - min_error)/delta;
    
  

  // Call the other refinement scheme
  this->flag_elements_by_error_fraction  (error_per_cell, top_frac,
					  bottom_frac, max_level);
}



void MeshRefinement::flag_elements_by_mean_stddev (const ErrorVector& error_per_cell_in,
						   const Real refine_fraction,
						   const Real coarsen_fraction,
						   const unsigned int max_level)
{
  // Copy the input error_per_cell so that we can modify it
  ErrorVector error_per_cell (error_per_cell_in);
  
//   // Set the error_per_cell to zero for elements at the
//   // maximum allowable refinement level
//   {
//     active_elem_iterator       elem_it (_mesh.elements_begin());
//     const active_elem_iterator elem_end(_mesh.elements_end());

//     for (; elem_it != elem_end; ++elem_it)
//       {
// 	Elem* elem               = *elem_it;
// 	const unsigned int id    = elem->id();
// 	const unsigned int level = elem->level();

// 	assert (id < error_per_cell.size());
      
// 	if (level >= max_level)
// 	  error_per_cell[id] = 0.;
//       }
//   }
  
  // Get the mean value from the error vector
  const Real mean = error_per_cell.mean();
  
  // Get the standard deviation.  This equals the
  // square-root of the variance
  const Real stddev = sqrt (error_per_cell.variance());
  
  // Check for valid fractions
  assert (refine_fraction  >= 0.);
  assert (coarsen_fraction >= 0.);

  // The refine and coarsen cutoff
  const Real refine_cutoff  =  mean + refine_fraction * stddev;
  const Real coarsen_cutoff = ((coarsen_fraction == 0.) ? 0. :
			       mean - coarsen_fraction * stddev);

  

  // Loop over the elements and flag them for coarsening or
  // refinement based on the element error
//   active_elem_iterator       elem_it (_mesh.elements_begin());
//   const active_elem_iterator elem_end(_mesh.elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.active_elements_end(); 

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem             = *elem_it;
      const unsigned int id  = elem->id();

      assert (id < error_per_cell.size());
      
      const float elem_error = error_per_cell[id];

      // Flag the element for coarsening if its error
      // is <= coarsen_fraction*delta + error_min
      if (elem_error <= coarsen_cutoff)
	{
	  elem->set_refinement_flag(Elem::COARSEN);
	}
      
      // Flag the element for refinement if its error
      // is >= refinement_cutoff.
      if (elem_error >= refine_cutoff)
	if (elem->level() < max_level)
	  elem->set_refinement_flag(Elem::REFINE);
    }
}



void MeshRefinement::clean_refinement_flags ()
{
  // Possibly clean up the refinement flags from
  // a previous step
//   elem_iterator       elem_it (_mesh.elements_begin());
//   const elem_iterator elem_end(_mesh.elements_end());

  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end(); 

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->set_refinement_flag(Elem::DO_NOTHING);
}

#endif
