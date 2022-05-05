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
  CoupledSystemQoI() = default;
  virtual ~CoupledSystemQoI() = default;

  virtual void init_qoi_count(System & sys);
  virtual void init_context(DiffContext &);

  virtual void postprocess(){}

  virtual void side_qoi_derivative(DiffContext & context,
                                   const QoISet & qois);

  virtual void side_qoi(DiffContext & context,
                        const QoISet & qois);

  virtual std::unique_ptr<DifferentiableQoI> clone()
  {
    DifferentiableQoI * my_clone = new CoupledSystemQoI;
    *my_clone = *this;
    return std::unique_ptr<DifferentiableQoI>(my_clone);
  }

protected:
  unsigned int u_var, p_var;
};

#endif // H_QOI_H
