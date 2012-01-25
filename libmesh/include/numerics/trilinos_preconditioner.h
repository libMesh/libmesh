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



#ifndef __trilinos_preconditioner_h__
#define __trilinos_preconditioner_h__

#include "libmesh_config.h"

#ifdef LIBMESH_HAVE_TRILINOS

// C++ includes

// Local includes
#include "preconditioner.h"
#include "libmesh_common.h"
#include "enum_solver_package.h"
#include "enum_preconditioner_type.h"
#include "reference_counted_object.h"
#include "libmesh.h"

#include "Epetra_Operator.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_ParameterList.hpp"

namespace libMesh
{

// forward declarations
template <typename T> class AutoPtr;
template <typename T> class SparseMatrix;
template <typename T> class NumericVector;
template <typename T> class ShellMatrix;

/**
 * This class provides an interface to the suite of preconditioners available
 * from Trilinos.
 *
 * @author David Andrs, 2011
 */

template <typename T>
class TrilinosPreconditioner :
  public Preconditioner<T>,
  public Epetra_Operator
{
public:

  /**
   *  Constructor. Initializes PetscPreconditioner data structures
   */
  TrilinosPreconditioner ();

  /**
   * Destructor.
   */
  virtual ~TrilinosPreconditioner ();

  /**
   * Computes the preconditioned vector "y" based on input "x".
   * Usually by solving Py=x to get the action of P^-1 x.
   */
  virtual void apply(const NumericVector<T> & x, NumericVector<T> & y);

  /**
   * Release all memory and clear data structures.
   */
  virtual void clear () {}

  /**
   * Initialize data structures if not done so already.
   */
  virtual void init ();

  void set_params(Teuchos::ParameterList & list);

  /**
   * Returns the actual Trilinos preconditioner object.
   */
  Epetra_FECrsMatrix * mat() { return _mat; }

  /**
   */
  void set_preconditioner_type (const PreconditionerType & preconditioner_type);

  /**
   * Compute the preconditioner. In Trilinos, we need to call this explicitly.
   */
  void compute();

protected:

  /**
   * Trilinos preconditioner
   */
  Epetra_Operator * _prec;

  /**
   * Trilinos matrix that's been pulled out of the _matrix object.
   */
  Epetra_FECrsMatrix * _mat;

  /**
   * Parameter list to be used for building the preconditioner
   */
  Teuchos::ParameterList _param_list;

  // Epetra_Operator interface
  virtual int SetUseTranspose(bool UseTranspose);
  virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
  virtual int ApplyInverse(const Epetra_MultiVector &r, Epetra_MultiVector &z) const;
  virtual double NormInf() const;
  virtual const char *Label() const;
  virtual bool UseTranspose() const;
  virtual bool HasNormInf() const;
  virtual const Epetra_Comm &Comm() const;
  virtual const Epetra_Map &OperatorDomainMap() const;
  virtual const Epetra_Map &OperatorRangeMap() const;
};




/*----------------------- inline functions ----------------------------------*/
template <typename T>
inline
TrilinosPreconditioner<T>::TrilinosPreconditioner () :
  Preconditioner<T>(),
  _prec(NULL),
  _mat(NULL)
{
}



template <typename T>
inline
TrilinosPreconditioner<T>::~TrilinosPreconditioner ()
{
  this->clear ();
}

} // namespace libMesh

#endif // #ifdef LIBMESH_HAVE_TRILINOS
#endif // #ifdef __trilinos_preconditioner_h__
