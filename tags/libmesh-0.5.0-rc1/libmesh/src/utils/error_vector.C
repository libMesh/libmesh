// $Id: error_vector.C,v 1.9 2005-02-22 22:17:43 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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


// C++ includes

// Local includes
#include "error_vector.h"
#include "libmesh_logging.h"



// ------------------------------------------------------------
// ErrorVector class member functions
float ErrorVector::minimum() const
{
  START_LOG ("minimum()", "ErrorVector");

  const unsigned int n = this->size();
  float min = 1.e30;
  
  for (unsigned int i=0; i<n; i++)
    {
      // Only positive (or zero) values in the error vector
      assert((*this)[i] >= 0.);
      min = std::min (min, (*this)[i]);
    }
  STOP_LOG ("minimum()", "ErrorVector");

  // ErrorVectors are for positive values
  assert (min >= 0.);
  
  return min;
}



Real ErrorVector::mean() const
{
  START_LOG ("mean()", "ErrorVector");
  
  const unsigned int n = this->size();

  Real mean  = 0;
  unsigned int nnz = 0;
  
  for (unsigned int i=0; i<n; i++)
    if ((*this)[i] != 0.)
      {
	mean += ( static_cast<Real>((*this)[i]) - mean ) / (nnz + 1);

	nnz++;
      }
  
  STOP_LOG ("mean()", "ErrorVector");
  
  return mean;
}




Real ErrorVector::median() 
{
  const unsigned int n   = this->size();
  
  if (n == 0)
    return 0.;
  

  // Build a StatisticsVector<float> containing
  // only our nonzero entries and take its mean
  StatisticsVector<float> sv;

  sv.reserve (n);

  for (unsigned int i=0; i<n; i++)
    if((*this)[i] != 0.)
      sv.push_back((*this)[i]);
  
  return sv.median();
}




Real ErrorVector::median() const
{
  ErrorVector ev = (*this);

  return ev.median();
}




Real ErrorVector::variance(const Real mean) const
{
  const unsigned int n   = this->size();
  
  START_LOG ("variance()", "ErrorVector");
  
  Real variance = 0;
  unsigned int nnz = 0;
  
  for (unsigned int i=0; i<n; i++)
    if ((*this)[i] != 0.)
      {
	const Real delta = ( static_cast<Real>((*this)[i]) - mean );
	variance += (delta * delta - variance) / (nnz + 1);

	nnz++;
      }
  
  STOP_LOG ("variance()", "ErrorVector");
  
  return variance;
}




std::vector<unsigned int> ErrorVector::cut_below(Real cut) const
{
  START_LOG ("cut_below()", "ErrorVector");
  
  const unsigned int n = this->size();
  
  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary 
  
  for (unsigned int i=0; i<n; i++)
    if ((*this)[i] != 0.)
      {
	if ((*this)[i] < cut)
	  {
	    cut_indices.push_back(i);
	  }
      }
  
  STOP_LOG ("cut_below()", "ErrorVector");
  
  return cut_indices;
}




std::vector<unsigned int> ErrorVector::cut_above(Real cut) const
{
  START_LOG ("cut_above()", "ErrorVector");
  
  const unsigned int n   = this->size();
  
  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary 
  
  for (unsigned int i=0; i<n; i++)
    if ((*this)[i] != 0.)
      {
	if ((*this)[i] > cut)
	  {
	    cut_indices.push_back(i);
	  }
      }
  
  STOP_LOG ("cut_above()", "ErrorVector");
  
  return cut_indices;
}
