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


#include "libmesh/quadrature_gm.h"
#include "libmesh/quadrature_gauss.h"

namespace libMesh
{

// See also the file:
// quadrature_gm_3D.C
// for additional implementation.




// Constructor
QGrundmann_Moller::QGrundmann_Moller(const unsigned int d,
                                     const Order o) : QBase(d,o)
{
}


// Destructor
QGrundmann_Moller::~QGrundmann_Moller()
{
}



void QGrundmann_Moller::init_1D(const ElemType /*type_in*/,
                                unsigned int p)
{
  QGauss gauss1D(1, static_cast<Order>(_order+2*p));

  // Swap points and weights with the about-to-be destroyed rule.
  _points.swap(gauss1D.get_points());
  _weights.swap(gauss1D.get_weights());
}





void QGrundmann_Moller::gm_rule(unsigned int s, unsigned int dim)
{
  // Only dim==2 or dim==3 is allowed
  libmesh_assert(dim==2 || dim==3);

  // A GM rule of index s can integrate polynomials of degree 2*s+1 exactly
  const unsigned int degree = 2*s+1;

  // The number of points for rule of index s is
  // (dim+1+s)! / (dim+1)! / s!
  // In 3D, this is = 1/24 * Pi_{i=1}^4 (s+i)
  // In 2D, this is =  1/6 * Pi_{i=1}^3 (s+i)
  const unsigned int n_pts = dim==2 ? (s+3)*(s+2)*(s+1) / 6 : (s+4)*(s+3)*(s+2)*(s+1) / 24;
  //libMesh::out << "n_pts=" << n_pts << std::endl;

  // Allocate space for points and weights
  _points.resize(n_pts);
  _weights.resize(n_pts);

  // (-1)^i -> This one flips sign at each iteration of the i-loop below.
  int one_pm=1;

  // Where we store all the integer point compositions/permutations
  std::vector<std::vector<unsigned int> > permutations;

  // Index into the vector where we should start adding the next round of points/weights
  std::size_t offset=0;

  // Implement the GM formula 4.1 on page 286 of the paper
  for (unsigned int i=0; i<=s; ++i)
    {
      // Get all the ordered compositions (and their permutations)
      // of |beta| = s-i into dim+1 parts
      compose_all(s-i, dim+1, permutations);
      //libMesh::out << "n. permutations=" << permutations.size() << std::endl;

      for (std::size_t p=0; p<permutations.size(); ++p)
        {
          // We use the first dim entries of each permutation to
          // construct an integration point.
          for (unsigned int j=0; j<dim; ++j)
            _points[offset+p](j) =
              static_cast<Real>(2.*permutations[p][j] + 1.) /
              static_cast<Real>(  degree + dim - 2.*i     );
        }

      // Compute the weight for this i, being careful to avoid overflow.
      // This technique is borrowed from Burkardt's code as well.
      // Use once for each of the points obtained from the permutations array.
      Real weight = one_pm;

      // This for loop needs to run for dim, degree, or dim+degree-i iterations,
      // whichever is largest.
      const unsigned int weight_loop_index =
        std::max(dim, std::max(degree, degree+dim-i))+1;

      for (unsigned int j=1; j<weight_loop_index; ++j)
        {
          if (j <= degree) // Accumulate (d+n-2i)^d term
            weight *= static_cast<Real>(degree+dim-2*i);

          if (j <= 2*s) // Accumulate 2^{-2s}
            weight *= 0.5;

          if (j <= i) // Accumulate (i!)^{-1}
            weight /= static_cast<Real>(j);

          if (j <= degree+dim-i) // Accumulate ( (d+n-i)! )^{-1}
            weight /= static_cast<Real>(j);
        }

      // This is the weight for each of the points computed previously
      for (std::size_t j=0; j<permutations.size(); ++j)
        _weights[offset+j] = weight;

      // Change sign for next iteration
      one_pm = -one_pm;

      // Update offset for the next set of points
      offset += permutations.size();
    }
}




// This routine for computing compositions and their permutations is
// originally due to:
//
// Albert Nijenhuis, Herbert Wilf,
// Combinatorial Algorithms for Computers and Calculators,
// Second Edition,
// Academic Press, 1978,
// ISBN: 0-12-519260-6,
// LC: QA164.N54.
//
// The routine is deceptively simple: I still don't understand exactly
// why it works, but it does.
void QGrundmann_Moller::compose_all(unsigned int s, // number to be compositioned
                                    unsigned int p, // # of partitions
                                    std::vector<std::vector<unsigned int> > & result)
{
  // Clear out results remaining from previous calls
  result.clear();

  // Allocate storage for a workspace.  The workspace will periodically
  // be copied into the result container.
  std::vector<unsigned int> workspace(p);

  // The first result is always (s,0,...,0)
  workspace[0] = s;
  result.push_back(workspace);

  // the value of the first non-zero entry
  unsigned int head_value=s;

  // When head_index=-1, it refers to "off the front" of the array.  Therefore,
  // this needs to be a regular int rather than unsigned.  I initially tried to
  // do this with head_index unsigned and an else statement below, but then there
  // is the special case: (1,0,...,0) which does not work correctly.
  int head_index = -1;

  // At the end, all the entries will be in the final slot of workspace
  while (workspace.back() != s)
    {
      // If the previous head value is still larger than 1, reset the index
      // to "off the front" of the array
      if (head_value > 1)
        head_index = -1;

      // Either move the index onto the front of the array or on to
      // the next value.
      head_index++;

      // Get current value of the head entry
      head_value = workspace[head_index];

      // Put a zero into the head_index of the array.  If head_index==0,
      // this will be overwritten in the next line with head_value-1.
      workspace[head_index] = 0;

      // The initial entry gets the current head value, minus 1.
      // If head_value > 1, the next loop iteration will start back
      // at workspace[0] again.
      libmesh_assert_greater (head_value, 0);
      workspace[0] = head_value - 1;

      // Increment the head+1 value
      workspace[head_index+1] += 1;

      // Save this composition in the results
      result.push_back(workspace);
    }
}

} // namespace libMesh
