// $Id: laspack_vector.C,v 1.6 2003-02-13 22:56:12 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#include "mesh_config.h"

#if defined(HAVE_LASPACK) && !defined(USE_COMPLEX_NUMBERS)


// C++ includes
#include <math.h>

// Local Includes
#include "laspack_vector.h"




// void LaspackVector::init (const NumericVector& v, const bool fast)
// {
//   error();
  
//   init (v.local_size(), v.size(), fast);

//   vec = reinterpret_cast<const LaspackVector&>(v).vec;
// }



Real LaspackVector::l1_norm () const
{
  assert(closed());
  
  return static_cast<Real>(Laspack::l1Norm_V(const_cast<Laspack::QVector*>(&_vec)));
}



Real LaspackVector::l2_norm () const
{
  assert(closed());
  
  return static_cast<Real>(Laspack::l2Norm_V(const_cast<Laspack::QVector*>(&_vec)));
}



Real LaspackVector::linfty_norm () const
{
  assert(closed());
  
  return static_cast<Real>(Laspack::MaxNorm_V(const_cast<Laspack::QVector*>(&_vec)));
}




NumericVector& LaspackVector::operator += (const NumericVector& v)
{
  assert(closed());
  
  add(1., v);
  
  return *this;
}




NumericVector& LaspackVector::operator -= (const NumericVector& v)
{
  assert(closed());
  
  add(-1., v);
  
  return *this;
}



void LaspackVector::add (const Complex v)
{
  for (unsigned int i=0; i<size(); i++)
    add (i, v);
}




void LaspackVector::add (const NumericVector& v)
{
  add (1., v);
}



void LaspackVector::add (const Complex a, const NumericVector& v_in)
{
  const LaspackVector& v = reinterpret_cast<const LaspackVector&>(v_in);
  
  assert(size() == v.size());

  for (unsigned int i=0; i<v.size(); i++)
    add (i, a*v(i));
}


void LaspackVector::scale (const Complex factor)
{
  assert (initialized());
  
  Laspack::Mul_SV (factor, &_vec);
}



NumericVector& 
LaspackVector::operator = (const Complex s)
{
  assert (initialized());

  Laspack::V_SetAllCmp (&_vec, s);
  
  return *this;
}



NumericVector&
LaspackVector::operator = (const NumericVector& v_in)
{
  const LaspackVector& v = reinterpret_cast<const LaspackVector&>(v_in);

  *this = v;

  return *this;
}



LaspackVector&
LaspackVector::operator = (const LaspackVector& v)
{
  assert (size() == v.size());

  if (size() != 0)    
    Laspack::Asgn_VV (const_cast<Laspack::QVector*>(&_vec),
		      const_cast<Laspack::QVector*>(&v._vec)
		      );
  return *this;
}



NumericVector&
LaspackVector::operator = (const std::vector<Complex>& v)
{
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if (size() == v.size())      
    for (unsigned int i=0; i<v.size(); i++)
      set (i, v[i]);
  
  else
    error();

  return *this;
}



void LaspackVector::localize (NumericVector& v_local_in) const
{
  LaspackVector& v_local =
    reinterpret_cast<LaspackVector&>(v_local_in);

  v_local = *this;
}



void LaspackVector::localize (NumericVector& v_local_in,
			      const std::vector<unsigned int>& send_list) const
{
  LaspackVector& v_local =
    reinterpret_cast<LaspackVector&>(v_local_in);

  assert (send_list.size() == v_local.size());

  v_local = *this;
}



void LaspackVector::localize (std::vector<Complex>& v_local) const

{
  v_local.resize(size());

  for (unsigned int i=0; i<v_local.size(); i++)
    v_local[i] = (*this)(i);  
}



void LaspackVector::localize_to_one (std::vector<Complex>& v_local,
				     const unsigned int pid) const
{
  assert (pid == 0);

  localize (v_local);
}



#endif // #ifdef HAVE_LASPACK
