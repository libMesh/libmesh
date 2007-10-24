// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <algorithm> // for std::min_element, std::max_element


// Local includes
#include "statistics.h"
#include "libmesh_logging.h"



// ------------------------------------------------------------
// StatisticsVector class member functions
template <typename T>
Real StatisticsVector<T>::l2_norm() const
{
  Real normsq = 0.;
  for (unsigned i = 0; i != this->size(); ++i)
    normsq += ((*this)[i] * (*this)[i]);

  return std::sqrt(normsq);
}


template <typename T>
T StatisticsVector<T>::minimum() const
{
  START_LOG ("minimum()", "StatisticsVector");

  const T min = *(std::min_element(this->begin(), this->end()));

  STOP_LOG ("minimum()", "StatisticsVector");

  return min;
}




template <typename T>
T StatisticsVector<T>::maximum() const
{
  START_LOG ("maximum()", "StatisticsVector");

  const T max = *(std::max_element(this->begin(), this->end()));
  
  STOP_LOG ("maximum()", "StatisticsVector");
  
  return max;
}




template <typename T>
Real StatisticsVector<T>::mean() const
{
  START_LOG ("mean()", "StatisticsVector");
  
  const unsigned int n = this->size();

  Real mean = 0;

  for (unsigned int i=0; i<n; i++)
    {
      mean += ( static_cast<Real>((*this)[i]) - mean ) / (i + 1);
    }
  
  STOP_LOG ("mean()", "StatisticsVector");
  
  return mean;
}




template <typename T>
Real StatisticsVector<T>::median() 
{
  const unsigned int n   = this->size();
  
  if (n == 0)
    return 0.;
  
  START_LOG ("median()", "StatisticsVector");
  
  std::sort(this->begin(), this->end());
  
  const unsigned int lhs = (n-1) / 2;
  const unsigned int rhs = n / 2;
  
  Real median = 0;
  
  
  if (lhs == rhs)
    {
      median = static_cast<Real>((*this)[lhs]);
    }
  
  else
    {
      median = ( static_cast<Real>((*this)[lhs]) + 
		 static_cast<Real>((*this)[rhs]) ) / 2.0;
    }
 
  STOP_LOG ("median()", "StatisticsVector");
  
  return median;
}




template <typename T>
Real StatisticsVector<T>::median() const
{
  StatisticsVector<T> sv = (*this);

  return sv.median();
}




template <typename T>
Real StatisticsVector<T>::variance(const Real mean) const
{
  const unsigned int n   = this->size();
  
  START_LOG ("variance()", "StatisticsVector");
  
  Real variance = 0;

  for (unsigned int i=0; i<n; i++)
    {
      const Real delta = ( static_cast<Real>((*this)[i]) - mean );
      variance += (delta * delta - variance) / (i + 1);
    }

  if (n > 1)
    variance *= static_cast<Real>(n) / static_cast<Real>(n - 1);
  
  STOP_LOG ("variance()", "StatisticsVector");
  
  return variance;
}


template <typename T>
    void StatisticsVector<T>::normalize()
{
  const unsigned int n   = this->size();
  const Real max = this->maximum();
  
  for (unsigned int i=0; i<n; i++)
  {
    (*this)[i] = static_cast<T>((*this)[i] / max);
  }
}


template <typename T>
void StatisticsVector<T>::histogram(std::vector<unsigned int>& bin_members,
				    unsigned int n_bins)
{
  const unsigned int n   = this->size();
  
  std::sort(this->begin(), this->end());
  
  T min      = this->minimum();
  T max      = this->maximum();
  T bin_size = (max - min) / static_cast<T>(n_bins); 

  START_LOG ("histogram()", "StatisticsVector");
  
  std::vector<T> bin_bounds(n_bins+1);
  for (unsigned int i=0; i<n_bins+1; i++)
    bin_bounds[i] = min + i * bin_size;
  
  // std::vector<unsigned int> bin_members(n_bins);
  bin_members.resize(n_bins+1);
  
  unsigned int data_index = 0;
  for (unsigned int j=1; j<n_bins+1; j++) 
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

  
  STOP_LOG ("histogram()", "StatisticsVector");
}



template <typename T>
void StatisticsVector<T>::histogram(std::vector<unsigned int>& bin_members,
				    unsigned int n_bins) const
{
  StatisticsVector<T> sv = (*this);
  
  return sv.histogram(bin_members, n_bins);
}




template <typename T>
std::vector<unsigned int> StatisticsVector<T>::cut_below(Real cut) const
{
  START_LOG ("cut_below()", "StatisticsVector");
  
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
  
  STOP_LOG ("cut_below()", "StatisticsVector");
  
  return cut_indices;
}




template <typename T>
std::vector<unsigned int> StatisticsVector<T>::cut_above(Real cut) const
{
  START_LOG ("cut_above()", "StatisticsVector");
  
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
  
  STOP_LOG ("cut_above()", "StatisticsVector");
  
  return cut_indices;
}




//------------------------------------------------------------
// Explicit Instantions
template class StatisticsVector<float>;
template class StatisticsVector<double>;
template class StatisticsVector<int>;
template class StatisticsVector<unsigned int>;
