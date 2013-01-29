#include "libmesh/meshfree_solution_transfer.h"

#include "libmesh/mesh.h"
#include "libmesh/system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/threads.h"
#include "libmesh/meshfree_interpolation.h"
#include "libmesh/function_base.h"

// C++ includes
#include <cstddef>

namespace libMesh {

// Forward Declarations
template <typename T>
class DenseVector;

// Helper function for doing the projection
class MeshlessInterpolationFunction : public FunctionBase<Number>
{  
public:
  MeshlessInterpolationFunction (const MeshfreeInterpolation &mfi,
				 Threads::spin_mutex &mutex) :
      _mfi(mfi),
      _mutex(mutex)
  {}

  void init () {}
  void clear () {}

  virtual AutoPtr<FunctionBase<Number> > clone () const
  { 
    return AutoPtr<FunctionBase<Number> > (new MeshlessInterpolationFunction (_mfi, _mutex) );
  }

  Number operator() (const Point& p,
		     const Real /*time*/)
  {
    _pts.clear();
    _pts.push_back(p);
    _vals.resize(1);

    Threads::spin_mutex::scoped_lock lock(_mutex);
    
    _mfi.interpolate_field_data(_mfi.field_variables(), _pts, _vals);

    return _vals.front();
  }


  void operator() (const Point& p,
		   const Real time,
		   DenseVector<Number>& output)
  {
    output.resize(1);
    output(0) = (*this)(p,time);
    return;
  }

private:
  const MeshfreeInterpolation &_mfi;
  mutable std::vector<Point> _pts;
  mutable std::vector<Number> _vals;
  Threads::spin_mutex &_mutex;
};

void
MeshfreeSolutionTransfer::transfer(const Variable & from_var, const Variable & to_var)
{  
  System * from_sys = from_var.system();
  System * to_sys = to_var.system();

  EquationSystems & from_es = from_sys->get_equation_systems();

  MeshBase & from_mesh = from_es.get_mesh();
   
  InverseDistanceInterpolation<LIBMESH_DIM> idi (4, 2);

  std::vector<Point>  &src_pts  (idi.get_source_points());
  std::vector<Number> &src_vals (idi.get_source_vals());
  
  std::vector<std::string> field_vars;      
  field_vars.push_back(from_var.name());
  idi.set_field_variables(field_vars);

  // We now will loop over every node in the source mesh
  // and add it to a source point list, along with the solution
  {
    MeshBase::const_node_iterator nd  = from_mesh.local_nodes_begin();
    MeshBase::const_node_iterator end = from_mesh.local_nodes_end();

    for (; nd!=end; ++nd)
    {
      const Node *node(*nd);
      src_pts.push_back(*node);
      src_vals.push_back((*from_sys->solution)(node->dof_number(from_sys->number(),from_var.number(),0)));
    }	  	
  }

  // We have only set local values - prepare for use by gathering remote gata
  idi.prepare_for_use();

  // Create a MeshlessInterpolationFunction that uses our InverseDistanceInterpolation
  // object.  Since each MeshlessInterpolationFunction shares the same InverseDistanceInterpolation
  // object in a threaded environment we must also provide a locking mechanism.
  Threads::spin_mutex mutex;
  MeshlessInterpolationFunction mif(idi, mutex);

  // project the solution
  to_sys->project_solution(&mif);      
}

} // namespace libMesh
