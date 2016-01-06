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

#ifndef LIBMESH_RB_THETA_EXPANSION_H
#define LIBMESH_RB_THETA_EXPANSION_H

// libMesh includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"

// C++ includes
#include <vector>


namespace libMesh
{

class RBTheta;
class RBParameters;

/**
 * This class stores the set of RBTheta functor objects that define
 * the "parameter-dependent expansion" of a PDE.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class RBThetaExpansion : public ReferenceCountedObject<RBThetaExpansion>
{
public:

  /**
   * Constructor.
   */
  RBThetaExpansion();

  /**
   * Destructor.
   */
  virtual ~RBThetaExpansion() {}

  /**
   * Evaluate theta_q_a at the current parameter. Overload
   * if the theta functions need to be treated differently
   * in subclasses.
   */
  virtual Number eval_A_theta(unsigned int q,
                              const RBParameters & mu);

  /**
   * Evaluate theta_q_f at the current parameter.
   */
  virtual Number eval_F_theta(unsigned int q,
                              const RBParameters & mu);

  /**
   * Evaluate theta_q_l at the current parameter.
   */
  virtual Number eval_output_theta(unsigned int output_index,
                                   unsigned int q_l,
                                   const RBParameters & mu);

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
   * Attach a pointer to a functor object that defines one
   * of the theta_q_a terms.
   */
  virtual void attach_A_theta(RBTheta * theta_q_a);

  /**
   * Attach a vector of pointers to functor objects that each define one
   * of the theta_q_a terms.
   */
  virtual void attach_multiple_A_theta(std::vector<RBTheta *> theta_q_a);

  /**
   * Attach a pointer to a functor object that defines one
   * of the theta_q_a terms.
   */
  virtual void attach_F_theta(RBTheta * theta_q_f);

  /**
   * Attach a vector of pointers to functor objects that each define one
   * of the theta_q_f terms.
   */
  virtual void attach_multiple_F_theta(std::vector<RBTheta *> theta_q_f);

  /**
   * Attach a vector of pointers to functor objects that define one
   * of the outputs.
   */
  virtual void attach_output_theta(std::vector<RBTheta *> theta_q_l);

  /**
   * Attach a pointer to a functor object that defines one
   * of the outputs.
   */
  virtual void attach_output_theta(RBTheta * theta_q_l);


private:

  /**
   * Vector storing the pointers to the RBTheta functors for A.
   */
  std::vector<RBTheta *> _A_theta_vector;

  /**
   * Vector storing the RBTheta functors for the affine expansion of the rhs.
   */
  std::vector<RBTheta *> _F_theta_vector;

  /**
   * Vector storing the RBTheta functors for the affine expansion of the outputs.
   */
  std::vector< std::vector<RBTheta *> > _output_theta_vector;

};

}

#endif // LIBMESH_RB_THETA_EXPANSION_H
