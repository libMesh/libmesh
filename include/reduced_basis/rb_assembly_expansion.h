// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_RB_ASSEMBLY_EXPANSION_H
#define LIBMESH_RB_ASSEMBLY_EXPANSION_H

// libMesh includes
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <vector>


namespace libMesh
{

// Forward declarations
class ElemAssembly;
class FEMContext;

/**
 * This class stores the set of ElemAssembly functor objects that define
 * the "parameter-independent expansion" of a PDE.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBAssemblyExpansion : public ReferenceCountedObject<RBAssemblyExpansion>
{
public:

  /**
   * Constructor.
   */
  RBAssemblyExpansion();

  /**
   * Destructor.
   */
  virtual ~RBAssemblyExpansion() {}

  /**
   * Perform the specified A interior assembly.
   */
  void perform_A_interior_assembly(unsigned int q,
                                   FEMContext & context);

  /**
   * Perform the specified A boundary assembly.
   */
  void perform_A_boundary_assembly(unsigned int q,
                                   FEMContext & context);

  /**
   * Perform the specified F interior assembly.
   */
  void perform_F_interior_assembly(unsigned int q,
                                   FEMContext & context);

  /**
   * Perform the specified F boundary assembly.
   */
  void perform_F_boundary_assembly(unsigned int q,
                                   FEMContext & context);

  /**
   * Perform the specified output assembly.
   */
  void perform_output_interior_assembly(unsigned int output_index,
                                        unsigned int q_l,
                                        FEMContext & context);

  /**
   * Perform the specified output assembly.
   */
  void perform_output_boundary_assembly(unsigned int output_index,
                                        unsigned int q_l,
                                        FEMContext & context);

  /**
   * Get Q_a, the number of terms in the affine
   * expansion for the bilinear form.
   */
  unsigned int get_n_A_terms() const;

  /**
   * Get Q_f, the number of terms in the affine
   * expansion for the right-hand side.
   */
  unsigned int get_n_F_terms() const;

  /**
   * Get n_outputs, the number output functionals.
   */
  unsigned int get_n_outputs() const;

  /**
   * Get the number of affine terms associated with the specified output.
   */
  unsigned int get_n_output_terms(unsigned int output_index) const;

  /**
   * Attach ElemAssembly object for the left-hand side
   * (both interior and boundary assembly).
   */
  void attach_A_assembly(ElemAssembly * Aq_assembly);

  /**
   * Attach multiple ElemAssembly objects for the left-hand side
   * (both interior and boundary assembly).
   */
  void attach_multiple_A_assembly(std::vector<ElemAssembly *> Aq_assembly);

  /**
   * Attach ElemAssembly object for the right-hand side
   * (both interior and boundary assembly).
   */
  void attach_F_assembly(ElemAssembly * Fq_assembly);

  /**
   * Attach multiple ElemAssembly objects for the right-hand side
   * (both interior and boundary assembly).
   */
  void attach_multiple_F_assembly(std::vector<ElemAssembly *> Fq_assembly);

  /**
   * Attach ElemAssembly object for an output
   * (both interior and boundary assembly).
   * In this case we pass in vector arguments to allow for Q_l > 1.
   */
  virtual void attach_output_assembly(std::vector<ElemAssembly *> output_assembly);

  /**
   * Attach ElemAssembly object for an output
   * (both interior and boundary assembly).
   * This function provides simpler syntax in the case that Q_l = 1; we
   * do not need to use a vector in this case.
   */
  virtual void attach_output_assembly(ElemAssembly * output_assembly);

  /**
   * Return a reference to the specified A_assembly object.
   */
  ElemAssembly & get_A_assembly(unsigned int q);

  /**
   * Return a reference to the specified F_assembly object.
   */
  ElemAssembly & get_F_assembly(unsigned int q);

  /**
   * Return a reference to the specified output assembly object.
   */
  ElemAssembly & get_output_assembly(unsigned int output_index, unsigned int q_l);

private:

  /**
   * Vectors storing the function pointers to the assembly
   * routines for the affine operators, both interior and boundary
   * assembly.
   */
  std::vector<ElemAssembly *> _A_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the rhs affine vectors.
   */
  std::vector<ElemAssembly *> _F_assembly_vector;

  /**
   * Vector storing the function pointers to the assembly
   * routines for the outputs. Element interior part.
   */
  std::vector< std::vector<ElemAssembly *> > _output_assembly_vector;
};

}

#endif // LIBMESH_RB_ASSEMBLY_EXPANSION_H
