// $Id: inf_fe.C,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



// Local includes
#include "inf_fe.h"



#ifdef ENABLE_INFINITE_ELEMENTS




// ------------------------------------------------------------
// InfFE class members
unsigned int InfFE::n_shape_functions() const
{
  return n_shape_functions(elem_type,
			   get_order());
};



unsigned int InfFE::n_shape_functions(const ElemType t, const Order o)
{
  return n_dofs(t, o);
};




unsigned int InfFE::n_dofs(const ElemType t, const Order o)
{
  // this one needs different switch statements, since different
  // Order descriptors are used.
  // if we do _not_ have this InfFE::n_dofs(), we would _REALLY_ have
  // to mess up FEBase::n_dofs().
  error();
  return 1;
};




unsigned int InfFE::n_dofs_at_node(const ElemType t,
				   const Order o,
				   const unsigned int n)
{
  // see n_dofs, why we need this.
  error();
  return 1;
};





unsigned int InfFE::n_dofs_per_elem(const ElemType t,
				    const Order o)
{
  // see n_dofs, why we need this.
  error();
  return 1;
};






void InfFE::reinit(const Elem* elem)
{
  // if we don't have this, we would call the wrong init_shape_functions,
  // compute_map() etc.
  error();
};



void InfFE::reinit(QBase* qside,
		   const Elem* elem,
		   const unsigned int side)
{
  // if we don't have this, we would call the wrong init_shape_functions,
  // compute_map() etc.
  error();
};






void InfFE::init_shape_functions(const QBase* qrule,
				 const Elem* elem)
{
  // if we don't re-define this here, we would call the
  // FEBase::shape(....), not the InfFE::shape(...).
  error();
};





void InfFE::init_shape_functions(const QBase* qrule,
				 const Elem* elem,
				 const unsigned int s)
{
  // see above, why we need this.
  error();
};





// Note that there is _no_ need to re-define
// FEBase::compute_shape_functions(). We can directly re-use
// this from FEBase.







//------------------------------------------------------------
// Mapping issues


bool InfFE::on_reference_element(const Point& p, 
				 const ElemType t, 
				 const real eps)
{
  // due to infinite extend, we must only in 2D (note that the INFHEX
  // only extend to r=2*a!)
  error();
  return false;
};




void InfFE::nodal_soln(const MeshBase& mesh,
		       const Elem* elem,
		       const Order order,
		       const std::vector<number>& elem_soln,
		       std::vector<number>&       nodal_soln)
{
  // since we use different polynomials, the nodal_soln is also
  // quite different from the FEBase version.
  error();
  return;
};




void InfFE::compute_map(const QBase* q,
			const Elem* e)
{
  // for conventional mapping this would be sufficient,
  // and since compute_map() only uses _true_ FEBase members, which
  // we do not need to overload here in InfFE, we would be save
  // to use it. Still, i need other coordinate mappings...
};

  
void InfFE::compute_map(const QBase* q,
			const Elem* e,
			const unsigned int s)
{
  // see above
};



//------------------------------------------------------------
// shape functions




real InfFE::shape(const unsigned int d,
		  const ElemType t,
		  const Order o,
		  const unsigned int i,
		  const Point& p)
{
  // obviously, these are different from FEBase
  std::cout << "This is InfFE::shape()!" << std::endl;
  return 0.;
};





real InfFE::shape(const unsigned int d,
		  const Elem* elem,
		  const Order o,
		  const unsigned int i,
		  const Point& p)
{
  // obviously, these are different from FEBase
  std::cout << "This is InfFE::shape()!" << std::endl;
  return 0.;
};





real InfFE::shape_deriv(const unsigned int d,
			const ElemType t,
			const Order o,
			const unsigned int i,
			const unsigned int j,
			const Point& p)
{
  // obviously, these are different from FEBase
  std::cout << "This is InfFE::shape_deriv()!" << std::endl;
  return 0.;
};





real InfFE::shape_deriv(const unsigned int d,
			const Elem* elem,
			const Order o,
			const unsigned int i,
			const unsigned int j,
			const Point& p)
{
  // obviously, these are different from FEBase
  std::cout << "This is InfFE::shape_deriv()!" << std::endl;
  return 0.;
};



#endif
