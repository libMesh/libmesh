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
#include <iomanip>

// Local includes
#include "libmesh/radial_basis_interpolation.h"
#include "libmesh/radial_basis_functions.h"
#include "libmesh/mesh_tools.h" // BoundingBox
#include "libmesh/libmesh_logging.h"
#include "libmesh/eigen_core_support.h"



namespace libMesh
{
//--------------------------------------------------------------------------------
// RadialBasisInterpolation methods
template <unsigned int KDDim, class RBF>
void RadialBasisInterpolation<KDDim,RBF>::clear()
{
  // Call base class clear method
  InverseDistanceInterpolation<KDDim>::clear();
}



template <unsigned int KDDim, class RBF>
void RadialBasisInterpolation<KDDim,RBF>::prepare_for_use()
{
  // Call base class methods for prep
  InverseDistanceInterpolation<KDDim>::prepare_for_use();
  InverseDistanceInterpolation<KDDim>::construct_kd_tree();

#ifndef LIBMESH_HAVE_EIGEN

  libmesh_error_msg("ERROR: this functionality presently requires Eigen!");

#else
  LOG_SCOPE ("prepare_for_use()", "RadialBasisInterpolation<>");

  // Construct a bounding box for our source points
  _src_bbox.invalidate();

  const std::size_t  n_src_pts = this->_src_pts.size();
  const unsigned int n_vars    = this->n_field_variables();
  libmesh_assert_equal_to (this->_src_vals.size(), n_src_pts*this->n_field_variables());

  {
    Point
      &p_min(_src_bbox.min()),
      &p_max(_src_bbox.max());

    for (std::size_t p=0; p<n_src_pts; p++)
      {
        const Point & p_src(_src_pts[p]);

        for (unsigned int d=0; d<LIBMESH_DIM; d++)
          {
            p_min(d) = std::min(p_min(d), p_src(d));
            p_max(d) = std::max(p_max(d), p_src(d));
          }
      }
  }

  libMesh::out << "bounding box is \n"
               << _src_bbox.min() << '\n'
               << _src_bbox.max() << std::endl;


  // Construct the Radial Basis Function, giving it the size of the domain
  if(_r_override < 0)
    _r_bbox = (_src_bbox.max() - _src_bbox.min()).norm();
  else
    _r_bbox = _r_override;

  RBF rbf(_r_bbox);

  libMesh::out << "bounding box is \n"
               << _src_bbox.min() << '\n'
               << _src_bbox.max() << '\n'
               << "r_bbox = " << _r_bbox << '\n'
               << "rbf(r_bbox/2) = " << rbf(_r_bbox/2) << std::endl;


  // Construct the projection Matrix
  typedef Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> DynamicMatrix;
  //typedef Eigen::Matrix<Number, Eigen::Dynamic,              1, Eigen::ColMajor> DynamicVector;

  DynamicMatrix A(n_src_pts, n_src_pts), x(n_src_pts,n_vars), b(n_src_pts,n_vars);

  for (std::size_t i=0; i<n_src_pts; i++)
    {
      const Point & x_i (_src_pts[i]);

      // Diagonal
      A(i,i) = rbf(0.);

      for (std::size_t j=i+1; j<n_src_pts; j++)
        {
          const Point & x_j (_src_pts[j]);

          const Real r_ij = (x_j - x_i).norm();

          A(i,j) = A(j,i) = rbf(r_ij);
        }

      // set source data
      for (unsigned int var=0; var<n_vars; var++)
        b(i,var) = _src_vals[i*n_vars + var];
    }


  // Solve the linear system
  x = A.ldlt().solve(b);
  //x = A.fullPivLu().solve(b);

  // save  the weights for each variable
  _weights.resize (this->_src_vals.size());

  for (std::size_t i=0; i<n_src_pts; i++)
    for (unsigned int var=0; var<n_vars; var++)
      _weights[i*n_vars + var] = x(i,var);

#endif

}



template <unsigned int KDDim, class RBF>
void RadialBasisInterpolation<KDDim,RBF>::interpolate_field_data (const std::vector<std::string> & field_names,
                                                                  const std::vector<Point> & tgt_pts,
                                                                  std::vector<Number> & tgt_vals) const
{
  LOG_SCOPE ("interpolate_field_data()", "RadialBasisInterpolation<>");

  libmesh_experimental();

  const unsigned int
    n_vars    = this->n_field_variables();

  const std::size_t
    n_src_pts = this->_src_pts.size(),
    n_tgt_pts = tgt_pts.size();

  libmesh_assert_equal_to (_weights.size(),    this->_src_vals.size());
  libmesh_assert_equal_to (field_names.size(), this->n_field_variables());

  // If we already have field variables, we assume we are appending.
  // that means the names and ordering better be identical!
  if (this->_names.size() != field_names.size())
    libmesh_error_msg("ERROR:  when adding field data to an existing list the \nvariable list must be the same!");

  for (std::size_t v=0; v<this->_names.size(); v++)
    if (_names[v] != field_names[v])
      libmesh_error_msg("ERROR:  when adding field data to an existing list the \nvariable list must be the same!");


  RBF rbf(_r_bbox);

  tgt_vals.resize (n_tgt_pts*n_vars); /**/ std::fill (tgt_vals.begin(), tgt_vals.end(), Number(0.));

  for (std::size_t tgt=0; tgt<n_tgt_pts; tgt++)
    {
      const Point & p (tgt_pts[tgt]);

      for (std::size_t i=0; i<n_src_pts; i++)
        {
          const Point & x_i(_src_pts[i]);
          const Real
            r_i   = (p - x_i).norm(),
            phi_i = rbf(r_i);

          for (unsigned int var=0; var<n_vars; var++)
            tgt_vals[tgt*n_vars + var] += _weights[i*n_vars + var]*phi_i;
        }
    }
}



// ------------------------------------------------------------
// Explicit Instantiations
template class RadialBasisInterpolation<3, WendlandRBF<3,0> >;
template class RadialBasisInterpolation<3, WendlandRBF<3,2> >;
template class RadialBasisInterpolation<3, WendlandRBF<3,4> >;
template class RadialBasisInterpolation<3, WendlandRBF<3,8> >;

} // namespace libMesh
