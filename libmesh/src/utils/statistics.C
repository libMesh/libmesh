// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include <fstream> // std::ofstream
#include <numeric> // std::accumulate

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
  // Must have at least 1 bin
  libmesh_assert (n_bins>0);

  const unsigned int n   = this->size();
  
  std::sort(this->begin(), this->end());

  // The StatisticsVector can hold both integer and float types.
  // We will define all the bins, etc. using Reals.
  Real min      = static_cast<Real>(this->minimum());
  Real max      = static_cast<Real>(this->maximum());
  Real bin_size = (max - min) / static_cast<Real>(n_bins); 

  START_LOG ("histogram()", "StatisticsVector");
  
  std::vector<Real> bin_bounds(n_bins+1);
  for (unsigned int i=0; i<bin_bounds.size(); i++)
    bin_bounds[i] = min + i * bin_size;

  // Give the last bin boundary a little wiggle room: we don't want
  // it to be just barely less than the max, otherwise our bin test below
  // may fail.
  bin_bounds.back() += 1.e-6 * bin_size;
  
  // This vector will store the number of members each bin has.
  bin_members.resize(n_bins);
  
  unsigned int data_index = 0;
  for (unsigned int j=0; j<bin_members.size(); j++) // bin vector indexing
    {
      // libMesh::out << "(debug) Filling bin " << j << std::endl;
      
      for (unsigned int i=data_index; i<n; i++) // data vector indexing
	{
	  //	libMesh::out << "(debug) Processing index=" << i << std::endl;
	  Real current_val = static_cast<Real>( (*this)[i] );
	  
	  // There may be entries in the vector smaller than the value
	  // reported by this->minimum().  (e.g. inactive elements in an
	  // ErrorVector.)  We just skip entries like that.
	  if ( current_val < min )
	    {
	      // 	    libMesh::out << "(debug) Skipping entry v[" << i << "]="
	      // 		      << (*this)[i]
	      // 		      << " which is less than the min value: min="
	      // 		      << min << std::endl;
	      continue;
	    }
	
	  if ( current_val > bin_bounds[j+1] ) // if outside the current bin (bin[j] is bounded
	                                       // by bin_bounds[j] and bin_bounds[j+1])
	    {
	      // libMesh::out.precision(16);
	      // 	    libMesh::out.setf(std::ios_base::fixed);
	      // 	    libMesh::out << "(debug) (*this)[i]= " << (*this)[i]
	      // 		      << " is greater than bin_bounds[j+1]="
	      //		      << bin_bounds[j+1]	 << std::endl;
	      data_index = i; // start searching here for next bin 
	      break; // go to next bin
	    }
	
	  // Otherwise, increment current bin's count
	  bin_members[j]++;
	  // libMesh::out << "(debug) Binned index=" << i << std::endl;
	}
    }
  
#ifdef DEBUG
  // Check the number of binned entries
  const unsigned int n_binned = std::accumulate(bin_members.begin(),
						bin_members.end(),
						static_cast<unsigned int>(0),
						std::plus<unsigned int>());

  if (n != n_binned)
    {
      libMesh::out << "Warning: The number of binned entries, n_binned="
		<< n_binned
		<< ", did not match the total number of entries, n="
		<< n << "." << std::endl;
      //libmesh_error();
    }
#endif

  
  STOP_LOG ("histogram()", "StatisticsVector");
}





template <typename T>
void StatisticsVector<T>::plot_histogram(const std::string& filename,
					 unsigned int n_bins)
{
  // First generate the histogram with the desired number of bins
  std::vector<unsigned int> bin_members;
  this->histogram(bin_members, n_bins);

  // The max, min and bin size are used to generate x-axis values.
  T min      = this->minimum();
  T max      = this->maximum();
  T bin_size = (max - min) / static_cast<T>(n_bins); 
  
  // On processor 0: Write histogram to file
  if (libMesh::processor_id()==0)
    {
      std::ofstream out (filename.c_str());

      out << "clear all\n";
      out << "clf\n";
      //out << "x=linspace(" << min << "," << max << "," << n_bins+1 << ");\n";

      // abscissa values are located at the center of each bin.
      out << "x=[";
      for (unsigned int i=0; i<bin_members.size(); ++i)
	{
	  out << min + (i+0.5)*bin_size << " ";
	}
      out << "];\n";
	
      out << "y=[";
      for (unsigned int i=0; i<bin_members.size(); ++i)
	{
	  out << bin_members[i] << " ";
	}
      out << "];\n";
      out << "bar(x,y);\n";
    }
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
#ifdef TRIPLE_PRECISION
template class StatisticsVector<long double>;
#endif
template class StatisticsVector<int>;
template class StatisticsVector<unsigned int>;
