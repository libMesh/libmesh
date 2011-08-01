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

#ifndef __rb_theta_expansion_h__
#define __rb_theta_expansion_h__

// libMesh includes
#include "libmesh_common.h"
#include "reference_counted_object.h"

// misc includes
#include <vector>


namespace libMesh
{

class RBTheta;

/**
 * This class stores the set of RBTheta functor objects that define
 * the "parameter-dependent expansion" of a PDE.
 *
 * @author David J. Knezevic, 2011
 */

// ------------------------------------------------------------
// RBThetaExpansion class definition
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
  virtual Number eval_theta_q_a(unsigned int q,
                                const std::vector<Real>& mu);

  /**
   * Evaluate theta_q_f at the current parameter.
   */
  virtual Number eval_theta_q_f(unsigned int q,
                                const std::vector<Real>& mu);

  /**
   * Evaluate theta_q_l at the current parameter.
   */
  Number eval_theta_q_l(unsigned int output_index,
                        unsigned int q_l,
                        const std::vector<Real>& mu);

  /**
   * Get Q_a, the number of terms in the affine
   * expansion for the bilinear form.
   */
  virtual unsigned int get_Q_a();

  /**
   * Get Q_f, the number of terms in the affine
   * expansion for the right-hand side.
   */
  unsigned int get_Q_f() const;

  /**
   * Get n_outputs, the number output functionals.
   */
  unsigned int get_n_outputs() const;
  
  /**
   * Get the number of affine terms associated with the specified output.
   */
  unsigned int get_Q_l(unsigned int output_index) const;

  /**
   * Attach a pointer to a functor object that defines one
   * of the theta_q_a terms.
   */
  virtual void attach_theta_q_a(RBTheta* theta_q_a);

  /**
   * Attach a pointer to a functor object that defines one
   * of the theta_q_a terms.
   */
  virtual void attach_theta_q_f(RBTheta* theta_q_f);

  /**
   * Attach a vector of pointers to functor objects that define one
   * of the outputs.
   */
  virtual void attach_output_theta(std::vector<RBTheta*> theta_q_l);

  /**
   * Attach a pointer to a functor object that defines one
   * of the outputs.
   */
  virtual void attach_output_theta(RBTheta* theta_q_l);


  // -------- Data members --------


  /**
   * Vector storing the pointers to the RBTheta functors.
   */
  std::vector<RBTheta*> theta_q_a_vector;

  /**
   * Vector storing the RBTheta functors for the theta_q_f (affine expansion of the rhs).
   */
  std::vector<RBTheta*> theta_q_f_vector;

  /**
   * Vector storing the RBTheta functors for the theta_q_l (affine expansion of the outputs).
   */
  std::vector< std::vector<RBTheta*> > theta_q_l_vector;

};

}

#endif