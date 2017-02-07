// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/statistics.h"
#include "libmesh/libmesh_logging.h"

namespace libMesh
{



// ------------------------------------------------------------
// StatisticsVector class member functions
template <typename T>
Real StatisticsVector<T>::l2_norm() const
{
  Real normsq = 0.;
  const dof_id_type n = cast_int<dof_id_type>(this->size());
  for (dof_id_type i = 0; i != n; ++i)
    normsq += ((*this)[i] * (*this)[i]);

  return std::sqrt(normsq);
}


template <typename T>
T StatisticsVector<T>::minimum() const
{
  LOG_SCOPE ("minimum()", "StatisticsVector");

  const T min = *(std::min_element(this->begin(), this->end()));

  return min;
}




template <typename T>
T StatisticsVector<T>::maximum() const
{
  LOG_SCOPE ("maximum()", "StatisticsVector");

  const T max = *(std::max_element(this->begin(), this->end()));

  return max;
}




template <typename T>
Real StatisticsVector<T>::mean() const
{
  LOG_SCOPE ("mean()", "StatisticsVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  Real the_mean = 0;

  for (dof_id_type i=0; i<n; i++)
    {
      the_mean += ( static_cast<Real>((*this)[i]) - the_mean ) /
        static_cast<Real>(i + 1);
    }

  return the_mean;
}




template <typename T>
Real StatisticsVector<T>::median()
{
  const dof_id_type n = cast_int<dof_id_type>(this->size());

  if (n == 0)
    return 0.;

  LOG_SCOPE ("median()", "StatisticsVector");

  std::sort(this->begin(), this->end());

  const dof_id_type lhs = (n-1) / 2;
  const dof_id_type rhs = n / 2;

  Real the_median = 0;


  if (lhs == rhs)
    {
      the_median = static_cast<Real>((*this)[lhs]);
    }

  else
    {
      the_median = ( static_cast<Real>((*this)[lhs]) +
                     static_cast<Real>((*this)[rhs]) ) / 2.0;
    }

  return the_median;
}




template <typename T>
Real StatisticsVector<T>::median() const
{
  StatisticsVector<T> sv = (*this);

  return sv.median();
}




template <typename T>
Real StatisticsVector<T>::variance(const Real mean_in) const
{
  const dof_id_type n = cast_int<dof_id_type>(this->size());

  LOG_SCOPE ("variance()", "StatisticsVector");

  Real the_variance = 0;

  for (dof_id_type i=0; i<n; i++)
    {
      const Real delta = ( static_cast<Real>((*this)[i]) - mean_in );
      the_variance += (delta * delta - the_variance) /
        static_cast<Real>(i + 1);
    }

  if (n > 1)
    the_variance *= static_cast<Real>(n) / static_cast<Real>(n - 1);

  return the_variance;
}


template <typename T>
void StatisticsVector<T>::normalize()
{
  const dof_id_type n = cast_int<dof_id_type>(this->size());
  const Real max = this->maximum();

  for (dof_id_type i=0; i<n; i++)
    (*this)[i] = static_cast<T>((*this)[i] / max);
}





template <typename T>
void StatisticsVector<T>::histogram(std::vector<dof_id_type> & bin_members,
                                    unsigned int n_bins)
{
  // Must have at least 1 bin
  libmesh_assert (n_bins>0);

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  std::sort(this->begin(), this->end());

  // The StatisticsVector can hold both integer and float types.
  // We will define all the bins, etc. using Reals.
  Real min      = static_cast<Real>(this->minimum());
  Real max      = static_cast<Real>(this->maximum());
  Real bin_size = (max - min) / static_cast<Real>(n_bins);

  LOG_SCOPE ("histogram()", "StatisticsVector");

  std::vector<Real> bin_bounds(n_bins+1);
  for (std::size_t i=0; i<bin_bounds.size(); i++)
    bin_bounds[i] = min + i * bin_size;

  // Give the last bin boundary a little wiggle room: we don't want
  // it to be just barely less than the max, otherwise our bin test below
  // may fail.
  bin_bounds.back() += 1.e-6 * bin_size;

  // This vector will store the number of members each bin has.
  bin_members.resize(n_bins);

  dof_id_type data_index = 0;
  for (std::size_t j=0; j<bin_members.size(); j++) // bin vector indexing
    {
      // libMesh::out << "(debug) Filling bin " << j << std::endl;

      for (dof_id_type i=data_index; i<n; i++) // data vector indexing
        {
          //libMesh::out << "(debug) Processing index=" << i << std::endl;
          Real current_val = static_cast<Real>( (*this)[i] );

          // There may be entries in the vector smaller than the value
          // reported by this->minimum().  (e.g. inactive elements in an
          // ErrorVector.)  We just skip entries like that.
          if ( current_val < min )
            {
              //     libMesh::out << "(debug) Skipping entry v[" << i << "]="
              //       << (*this)[i]
              //       << " which is less than the min value: min="
              //       << min << std::endl;
              continue;
            }

          if ( current_val > bin_bounds[j+1] ) // if outside the current bin (bin[j] is bounded
            // by bin_bounds[j] and bin_bounds[j+1])
            {
              // libMesh::out.precision(16);
              //     libMesh::out.setf(std::ios_base::fixed);
              //     libMesh::out << "(debug) (*this)[i]= " << (*this)[i]
              //       << " is greater than bin_bounds[j+1]="
              //      << bin_bounds[j+1] << std::endl;
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
  const dof_id_type n_binned = std::accumulate(bin_members.begin(),
                                               bin_members.end(),
                                               static_cast<dof_id_type>(0),
                                               std::plus<dof_id_type>());

  if (n != n_binned)
    {
      libMesh::out << "Warning: The number of binned entries, n_binned="
                   << n_binned
                   << ", did not match the total number of entries, n="
                   << n << "." << std::endl;
    }
#endif
}





template <typename T>
void StatisticsVector<T>::plot_histogram(const processor_id_type my_procid,
                                         const std::string & filename,
                                         unsigned int n_bins)
{
  // First generate the histogram with the desired number of bins
  std::vector<dof_id_type> bin_members;
  this->histogram(bin_members, n_bins);

  // The max, min and bin size are used to generate x-axis values.
  T min      = this->minimum();
  T max      = this->maximum();
  T bin_size = (max - min) / static_cast<T>(n_bins);

  // On processor 0: Write histogram to file
  if (my_procid==0)
    {
      std::ofstream out_stream (filename.c_str());

      out_stream << "clear all\n";
      out_stream << "clf\n";
      //out_stream << "x=linspace(" << min << "," << max << "," << n_bins+1 << ");\n";

      // abscissa values are located at the center of each bin.
      out_stream << "x=[";
      for (std::size_t i=0; i<bin_members.size(); ++i)
        {
          out_stream << min + (i+0.5)*bin_size << " ";
        }
      out_stream << "];\n";

      out_stream << "y=[";
      for (std::size_t i=0; i<bin_members.size(); ++i)
        {
          out_stream << bin_members[i] << " ";
        }
      out_stream << "];\n";
      out_stream << "bar(x,y);\n";
    }
}



template <typename T>
void StatisticsVector<T>::histogram(std::vector<dof_id_type> & bin_members,
                                    unsigned int n_bins) const
{
  StatisticsVector<T> sv = (*this);

  return sv.histogram(bin_members, n_bins);
}




template <typename T>
std::vector<dof_id_type> StatisticsVector<T>::cut_below(Real cut) const
{
  LOG_SCOPE ("cut_below()", "StatisticsVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  std::vector<dof_id_type> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (dof_id_type i=0; i<n; i++)
    {
      if ((*this)[i] < cut)
        {
          cut_indices.push_back(i);
        }
    }

  return cut_indices;
}




template <typename T>
std::vector<dof_id_type> StatisticsVector<T>::cut_above(Real cut) const
{
  LOG_SCOPE ("cut_above()", "StatisticsVector");

  const dof_id_type n = cast_int<dof_id_type>(this->size());

  std::vector<dof_id_type> cut_indices;
  cut_indices.reserve(n/2);  // Arbitrary

  for (dof_id_type i=0; i<n; i++)
    if ((*this)[i] > cut)
      cut_indices.push_back(i);

  return cut_indices;
}




//------------------------------------------------------------
// Explicit Instantions
template class StatisticsVector<float>;
template class StatisticsVector<double>;
#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
template class StatisticsVector<long double>;
#endif
template class StatisticsVector<int>;
template class StatisticsVector<unsigned int>;

} // namespace libMesh
