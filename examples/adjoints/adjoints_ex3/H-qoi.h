#ifndef H_QOI_H
#define H_QOI_H

#include "libmesh/libmesh_common.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/point.h"
#include "libmesh/quadrature.h"
#include "libmesh/diff_qoi.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class CoupledSystemQoI : public DifferentiableQoI
{
public:
  CoupledSystemQoI(){}
  virtual ~CoupledSystemQoI(){}

  virtual void init_qoi( std::vector<Number>& sys_qoi);
  virtual void postprocess( ){}

  virtual void side_qoi_derivative(DiffContext &context, const QoISet & qois);

  virtual void side_qoi(DiffContext &context, const QoISet & qois);

  virtual UniquePtr<DifferentiableQoI> clone()
  {
    DifferentiableQoI* my_clone = new CoupledSystemQoI;
    *my_clone = *this;
    return UniquePtr<DifferentiableQoI>(my_clone);
  }

};
#endif // H_QOI_H
