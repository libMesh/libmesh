// $Id: statistics.C,v 1.6 2003-02-03 03:51:50 ddreyer Exp $

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


// C++ includes
#include <algorithm>


// Local includes
#include "statistics.h"


// ------------------------------------------------------------
// StatisticsVector class member functions
template <typename T>
Real StatisticsVector<T>::mean() const
{
  const unsigned int n = this->size();

  Real mean = 0;

  for (unsigned int i=0; i<n; i++)
    {
      mean += ( static_cast<Real>((*this)[i]) - mean ) / (i + 1);
    }
  
  return mean;
}




template <typename T>
Real StatisticsVector<T>::median() 
{
  std::sort(this->begin(), this->end());
  
  const unsigned int n   = this->size();
  const unsigned int lhs = (n-1) / 2;
  const unsigned int rhs = n / 2;
  
  Real median = 0;
  
  if (n == 0)
    return median;
  
  if (lhs == rhs)
    {
      median = static_cast<Real>((*this)[lhs]);
    }
  
  else
    {
      median = ( static_cast<Real>((*this)[lhs]) + 
		 static_cast<Real>((*this)[rhs]) ) / 2.0;
    }
  
  return median;
}




template <typename T>
Real StatisticsVector<T>::median() const
{
  StatisticsVector<T> sv = (*this);

  return sv.median();
}




template <typename T>
Real StatisticsVector<T>::variance() const
{
  const unsigned int n   = this->size();
  
  const Real mean = this->mean();
  
  Real variance = 0;

  for (unsigned int i=0; i<n; i++)
    {
      const Real delta = ( static_cast<Real>((*this)[i]) - mean );
      variance += (delta * delta - variance) / (i + 1);
    }
  
  return variance;
}




template <typename T>
std::vector<unsigned int> StatisticsVector<T>::histogram(unsigned int n_bins)
{
  const unsigned int n   = this->size();
  
  std::sort(this->begin(), this->end());
  
  T max      = *(std::max_element(this->begin(), this->end())); 
  T min      = *(std::min_element(this->begin(), this->end()));
  T bin_size = (max - min) / static_cast<T>(n_bins); 

  std::vector<T> bin_bounds(n_bins+1);
  for (unsigned int i=0; i<n_bins+1; i++)
      bin_bounds[i] = min + i * bin_size;
  
  std::vector<unsigned int> bin_members(n_bins);

  unsigned int data_index = 0;
  for (unsigned int j=1; j<n_bins+1; j++) 
    {
      for (unsigned int i=data_index; i<n; i++) 
	{
	  if ( (*this)[i] > bin_bounds[j] ) 
	    {
	      data_index = i+1; 
	      bin_members[j]++; 
	      break;
	    }
	  
	  bin_members[j-1]++; 
	}
    }
  
  return bin_members; 
}



template <typename T>
std::vector<unsigned int> StatisticsVector<T>::histogram(unsigned int n_bins) const
{
  StatisticsVector<T> sv = (*this);
  
  return sv.histogram(n_bins);
}




template <typename T>
std::vector<unsigned int> StatisticsVector<T>::cut_below(Real cut) const
{
  const unsigned int n   = this->size();
  
  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary 
  
  for (unsigned int i=0; i<n; i++)
    {
      if ((*this)[i] < cut)
	{
	  cut_indices.push_back(i);
	}
    }
  
  return cut_indices;
}




template <typename T>
std::vector<unsigned int> StatisticsVector<T>::cut_above(Real cut) const
{
  const unsigned int n   = this->size();
  
  std::vector<unsigned int> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary 
  
  for (unsigned int i=0; i<n; i++)
    {
      if ((*this)[i] > cut)
	{
	  cut_indices.push_back(i);
	}
    }
  
  return cut_indices;
}




// Explicit Instantions
template class StatisticsVector<Real>;
template class StatisticsVector<int>;
template class StatisticsVector<unsigned int>;
