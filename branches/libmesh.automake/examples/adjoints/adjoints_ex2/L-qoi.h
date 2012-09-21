
#ifndef L_QOI_H
#define L_QOI_H

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/diff_qoi.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class LaplaceQoI : public DifferentiableQoI
{
public:
  LaplaceQoI(){}
  virtual ~LaplaceQoI(){} 

  virtual void init_qoi( std::vector<Number>& sys_qoi );

  virtual void postprocess( ){} 
  
  virtual void element_qoi_derivative(DiffContext &context, const QoISet & qois);  

  virtual void element_qoi (DiffContext &context, const QoISet & qois); 

  virtual AutoPtr<DifferentiableQoI> clone( ) {
    AutoPtr<DifferentiableQoI> my_clone( new LaplaceQoI );
    *my_clone = *this;
    return my_clone;
  }

};
#endif // L_QOI_H
