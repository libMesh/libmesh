// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/point.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/parallel.h"



namespace libMesh
{
  
  //--------------------------------------------------------------------------------
  // MeshfreeInterpolation methods
  void MeshfreeInterpolation::print_info (std::ostream& os) const
  {
    os << "MeshfreeInterpolation"
       << "\n n_field_variables()=" << this->n_field_variables() 
       << "\n n_source_points()=" << _src_pts.size()
       << std::endl;
  }


  
  std::ostream& operator << (std::ostream& os, const MeshfreeInterpolation& mfi)
  {
    mfi.print_info(os);
    return os;
  }



  void MeshfreeInterpolation::clear ()
  {
    _is_prepared = false;
    _names.clear();
    _src_pts.clear();
    _src_vals.clear();
  }



  void MeshfreeInterpolation::add_field_data (const std::vector<std::string> &field_names,
					      const std::vector<Point>  &pts,
					      const std::vector<Number> &vals)
  {
    _is_prepared = false;
    
    libmesh_assert_equal_to (field_names.size()*pts.size(), vals.size());

    // If we already have field variables, we assume we are appending.
    // that means the names and ordering better be identical!
    if (!_names.empty())
      {
	if (_names.size() != field_names.size())
	  {
	    libMesh::err << "ERROR:  when adding field data to an existing list the\n"
			 << "varaible list must be the same!\n";
	    libmesh_error();	      	      
	  }
	for (unsigned int v=0; v<_names.size(); v++)
	  if (_names[v] != field_names[v])
	    {
	      libMesh::err << "ERROR:  when adding field data to an existing list the\n"
			   << "varaible list must be the same!\n";
	      libmesh_error();	      	      
	    }	    
      }
    
    // otherwise copy the names
    else
      _names = field_names;

    // append the data
    _src_pts.insert (_src_pts.end(),
		     pts.begin(), 
		     pts.end());

    _src_vals.insert (_src_vals.end(),
		      vals.begin(),
		      vals.end());

    libmesh_assert_equal_to (_src_vals.size(), 
			     _src_pts.size()*this->n_field_variables());
  }


  
  void MeshfreeInterpolation::prepare_for_use ()
  {
    switch (_parallelization_strategy)
      {
      case SYNC_SOURCES:
	this->gather_remote_data();
	break;

      case INVALID_STRATEGY:
	libmesh_error();
	break;

      default:
	// no preparation required for other strategies
	break;
      }
    _is_prepared = true;
  }



  void MeshfreeInterpolation::gather_remote_data ()
  {
#ifndef LIBMESH_HAVE_MPI
    
    // no MPI -- no-op
    return;

#else

    // This function must be run on all processors at once
    parallel_only();

    START_LOG ("gather_remote_data()", "MeshfreeInterpolation");

    std::vector<Real> send_buf, recv_buf;
    
    libmesh_assert_equal_to (_src_vals.size(),
			     _src_pts.size()*this->n_field_variables());
    
    send_buf.reserve (_src_pts.size()*(3 + this->n_field_variables()));

    // Everyone packs their data at the same time
    for (unsigned int p_idx=0, v_idx=0; p_idx<_src_pts.size(); p_idx++)
      {
	const Point &pt(_src_pts[p_idx]);
	
	send_buf.push_back(pt(0));
	send_buf.push_back(pt(1));
	send_buf.push_back(pt(2));

	for (unsigned int var=0; var<this->n_field_variables(); var++)
	  {
	    libmesh_assert_less (v_idx, _src_vals.size());
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
	    send_buf.push_back (_src_vals[v_idx].real()); 
	    send_buf.push_back (_src_vals[v_idx].imag()); 
	    v_idx++;
	    
#else
	    send_buf.push_back (_src_vals[v_idx++]); 
#endif
  	  }
      }

    // Send our data to everyone else.  Note that MPI-1 said you could not 
    // use the same buffer in nonblocking sends, but that restriction has 
    // recently been removed.
    std::vector<Parallel::Request> send_request(libMesh::n_processors()-1);
    
    // Use a tag for best practices.  In debug mode parallel_only() blocks above
    // so we can be sure there is no other shenanigarry going on, but in optimized
    // mode there is no such guarantee - other prcoessors could be somewhere else
    // completing some other communication, and we don't want to intercept that.
    Parallel::MessageTag tag = Parallel::Communicator_World.get_unique_tag ( 6000 );

    for (unsigned int proc=0, cnt=0; proc<libMesh::n_processors(); proc++)
      if (proc != libMesh::processor_id())
	CommWorld.send (proc, send_buf, send_request[cnt++], tag);
	  
    // All data has been sent.  Receive remote data in any order
    for (unsigned int comm_step=0; comm_step<(libMesh::n_processors()-1); comm_step++)
      {
	// blocking receive
	CommWorld.receive (Parallel::any_source, recv_buf, tag);

	// Add their data to our list
	Point  pt;
	Number val;
	std::vector<Real>::const_iterator it=recv_buf.begin();
	while (it != recv_buf.end())
	  {
	    pt(0) = *it, ++it;
	    pt(1) = *it, ++it;
	    pt(2) = *it, ++it;

	    _src_pts.push_back(pt);
	    
	    for (unsigned int var=0; var<this->n_field_variables(); var++)
	      {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
		val.real() = *it, ++it;
		val.imag() = *it, ++it;
#else
		val = *it, ++it;
#endif
		_src_vals.push_back(val);
	      }
	  }
      }

    Parallel::wait (send_request);

    STOP_LOG  ("gather_remote_data()", "MeshfreeInterpolation");

#endif // LIBMESH_HAVE_MPI
  }



  //--------------------------------------------------------------------------------
  // InverseDistanceInterpolation methods
  template <unsigned int KDDim>
  void InverseDistanceInterpolation<KDDim>::construct_kd_tree ()
  {
#ifdef LIBMESH_HAVE_NANOFLANN

    START_LOG ("construct_kd_tree()", "InverseDistanceInterpolation<>");

    // Initialize underlying KD tree
    if (_kd_tree.get() == NULL)
      _kd_tree.reset (new kd_tree_t (KDDim,
				     _point_list_adaptor,
				     nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)));
    
    libmesh_assert (_kd_tree.get() != NULL);
    
    _kd_tree->buildIndex();

    STOP_LOG ("construct_kd_tree()", "InverseDistanceInterpolation<>");
#endif
  }



  template <unsigned int KDDim>
  void InverseDistanceInterpolation<KDDim>::interpolate_field_data (const std::vector<std::string> &field_names,
								    const std::vector<Point>  &tgt_pts,
								    std::vector<Number> &tgt_vals) const
  {
    // forcibly initialize, if needed
#if LIBMESH_HAVE_NANOFLANN
    if (_kd_tree.get() == NULL)
      const_cast<InverseDistanceInterpolation<KDDim>*>(this)->construct_kd_tree();
#endif

    START_LOG ("interpolate_field_data()", "InverseDistanceInterpolation<>");

    libmesh_assert_equal_to (field_names.size(), this->n_field_variables());

    // If we already have field variables, we assume we are appending.
    // that means the names and ordering better be identical!
    if (_names.size() != field_names.size())
      {
	libMesh::err << "ERROR:  when adding field data to an existing list the\n"
		     << "varaible list must be the same!\n";
	libmesh_error();	      	      
      }
    for (unsigned int v=0; v<_names.size(); v++)
      if (_names[v] != field_names[v])
	{
	  libMesh::err << "ERROR:  when adding field data to an existing list the\n"
		       << "varaible list must be the same!\n";
	  libmesh_error();	      	      
	}	    
    
    tgt_vals.resize (tgt_pts.size()*this->n_field_variables());
    
    std::vector<Number>::iterator out_it = tgt_vals.begin();

#ifdef LIBMESH_HAVE_NANOFLANN
    {
      const size_t num_results = std::min((size_t) _n_interp_pts, _src_pts.size());

      std::vector<size_t> ret_index(num_results);
      std::vector<Real>   ret_dist_sqr(num_results);

      for (std::vector<Point>::const_iterator tgt_it=tgt_pts.begin(); 
	   tgt_it != tgt_pts.end(); ++tgt_it)
	{
	  const Point &tgt(*tgt_it);
	  const Real query_pt[] = { tgt(0), tgt(1), tgt(2) };

	  _kd_tree->knnSearch(&query_pt[0], num_results, &ret_index[0], &ret_dist_sqr[0]);

	  this->interpolate (tgt, ret_index, ret_dist_sqr, out_it);

	  // std::cout << "knnSearch(): num_results=" << num_results << "\n";
	  // for (size_t i=0;i<num_results;i++)
	  //   std::cout << "idx[" << i << "]=" 
	  // 	      << std::setw(6) << ret_index[i] 
	  // 	      << "\t dist["<< i << "]=" << ret_dist_sqr[i] 
	  // 	      << "\t val[" << std::setw(6) << ret_index[i] << "]=" << _src_vals[ret_index[i]]
	  // 	      << std::endl;
	  // std::cout << "\n";

	  // std::cout << "ival=" << _vals[0] << '\n';
	}
    }    
#else
    
    libMesh::err << "ERROR:  This functionality requires the library to be configured\n"
		 << "with nanoflann KD-Tree approximate nearest neighbor support!\n"
		 << std::endl;
    libmesh_error();

#endif

    STOP_LOG ("interpolate_field_data()", "InverseDistanceInterpolation<>");
  }

  template <unsigned int KDDim>
  void InverseDistanceInterpolation<KDDim>::interpolate (const Point               & /* pt */,
							 const std::vector<size_t> &src_indices,
							 const std::vector<Real>   &src_dist_sqr,
							 std::vector<Number>::iterator &out_it) const
  {
    // We explicitly assume that the input source points are sorted from closest to
    // farthests.  assert that assumption in DEBUG mode.
#ifdef DEBUG
    if (!src_dist_sqr.empty())
      {
	Real min_dist = src_dist_sqr.front();
	std::vector<Real>::const_iterator it = src_dist_sqr.begin();

	for (++it; it!= src_dist_sqr.end(); ++it)
	  {
	    if (*it < min_dist) libmesh_error();
	    min_dist = *it;
	  }
      }    
#endif


    libmesh_assert_equal_to (src_dist_sqr.size(), src_indices.size());


    // Compute the interpolation weights & interpolated value
    const unsigned int n_fv = this->n_field_variables();
    _vals.resize(n_fv); /**/ std::fill (_vals.begin(), _vals.end(), Number(0.));
    
    Real tot_weight = 0.;

    std::vector<Real>::const_iterator src_dist_sqr_it=src_dist_sqr.begin();
    std::vector<size_t>::const_iterator src_idx_it=src_indices.begin();

    // Loop over source points
    while ((src_dist_sqr_it != src_dist_sqr.end()) &&
	   (src_idx_it      != src_indices.end()))
      {
	libmesh_assert_greater_equal (*src_dist_sqr_it, 0.);
	
	const Real
	  dist_sq = std::max(*src_dist_sqr_it, std::numeric_limits<Real>::epsilon()),
	  weight = 1./std::pow(dist_sq, _half_power);

	tot_weight += weight;

	const unsigned int src_idx = *src_idx_it;

	// loop over field variables
	for (unsigned int v=0; v<n_fv; v++)
	  {
	    libmesh_assert_less (src_idx*n_fv+v, _src_vals.size());
	    _vals[v] += _src_vals[src_idx*n_fv+v]*weight;
	  }

	++src_dist_sqr_it;
	++src_idx_it;
      }

    // don't forget normalizing term & set the output buffer!
    for (unsigned int v=0; v<n_fv; v++, ++out_it)
      {	
	_vals[v] /= tot_weight;
	
	*out_it = _vals[v];
      }	
  }


  
// ------------------------------------------------------------
// Explicit Instantiations
  template class InverseDistanceInterpolation<1>;
  template class InverseDistanceInterpolation<2>;
  template class InverseDistanceInterpolation<3>;

} // namespace libMesh

