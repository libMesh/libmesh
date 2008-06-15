// $Id: epetra_vector.C 2606 2008-01-23 20:21:47Z roystgnr $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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

// Local Includes
#include "trilinos_epetra_vector.h"

#ifdef HAVE_TRILINOS

#include "parallel.h"
#include "utility.h"
#include "dense_vector.h"
#include "parallel.h"

// Trilinos Includes
#include <Epetra_Import.h>



template <typename T>
T EpetraVector<T>::sum () const
{
  libmesh_assert(this->closed());
  
  const unsigned int nl = _vec->MyLength();

  T sum=0.0;

  T * values = _vec->Values();
  
  for (unsigned int i=0; i<nl; i++)
    sum += values[i];

  Parallel::sum<T>(sum);
  
  return sum;
}

template <typename T>
Real EpetraVector<T>::l1_norm () const
{
  libmesh_assert(this->closed());

  Real value;

  _vec->Norm1(&value);
  
  return value;
}

template <typename T>
Real EpetraVector<T>::l2_norm () const
{
  libmesh_assert(this->closed());
  
  Real value;

  _vec->Norm2(&value);
  
  return value;
}

template <typename T>
Real EpetraVector<T>::linfty_norm () const
{
  libmesh_assert(this->closed());
  
  Real value;

  _vec->NormInf(&value);
  
  return value;
}

template <typename T>
NumericVector<T>&
EpetraVector<T>::operator += (const NumericVector<T>& v)
{
  libmesh_assert(this->closed());
  
  this->add(1., v);
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator -= (const NumericVector<T>& v)
{
  libmesh_assert(this->closed());
  
  this->add(-1., v);
  
  return *this;
}



template <typename T>
void EpetraVector<T>::set (const unsigned int i_in, const T value_in)
{
  int i = static_cast<int> (i_in);
  T value = value_in;

  libmesh_assert(i_in<this->size());

  _vec->ReplaceGlobalValues(1,&value, &i);

  this->_is_closed = false;
}



template <typename T>
void EpetraVector<T>::add (const unsigned int i_in, const T value_in)
{
  int i = static_cast<int> (i_in);
  T value = value_in;

  libmesh_assert(i_in<this->size());
  
  _vec->SumIntoGlobalValues(1,&value,&i);

  this->_is_closed = false;
}



template <typename T>
void EpetraVector<T>::add_vector (const std::vector<T>& v,
				  const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (v.size() == dof_indices.size());

  _vec->SumIntoGlobalValues(v.size(),
			    const_cast<T*>(&v[0]),
			    (int*) &dof_indices[0]);
}



template <typename T>
void EpetraVector<T>::add_vector (const NumericVector<T>& V,
				 const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  for (unsigned int i=0; i<V.size(); i++)
    this->add (dof_indices[i], V(i));
}



// TODO: fill this in after creating an EpetraMatrix
template <typename T>
void EpetraVector<T>::add_vector (const NumericVector<T>& /* V_in */,
				  const SparseMatrix<T>& /* A_in */)
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   const EpetraVector<T>* V = dynamic_cast<const EpetraVector<T>*>(&V_in);
//   const EpetraMatrix<T>* A = dynamic_cast<const EpetraMatrix<T>*>(&A_in);

//   libmesh_assert (V != NULL);
//   libmesh_assert (A != NULL);
  
//   int ierr=0;

//   A->close();

//   // The const_cast<> is not elegant, but it is required since Epetra
//   // is not const-correct.  
//   ierr = MatMultAdd(const_cast<EpetraMatrix<T>*>(A)->mat(), V->_vec, _vec, _vec);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr); 
}



template <typename T>
void EpetraVector<T>::add_vector (const DenseVector<T>& /* V_in */,
				  const std::vector<unsigned int>& /* dof_indices */)
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   libmesh_assert (V_in.size() == dof_indices.size());

//   const EpetraVector<T>* V = dynamic_cast<const EpetraVector<T>*>(&V_in);

//   this->add_vector(V->_vec->get_values[0], dof_indices);
}



template <typename T>
void EpetraVector<T>::add (const T v_in)
{
  const unsigned int nl = _vec->MyLength();

  T * values = _vec->Values();
  
  for (unsigned int i=0; i<nl; i++)
    values[i]+=v_in;

  this->_is_closed = false;
}


template <typename T>
void EpetraVector<T>::add (const NumericVector<T>& v)
{
  this->add (1., v);
}


template <typename T>
void EpetraVector<T>::add (const T a_in, const NumericVector<T>& v_in)
{
  const EpetraVector<T>* v = dynamic_cast<const EpetraVector<T>*>(&v_in);

  libmesh_assert (v != NULL);
  libmesh_assert(this->size() == v->size());

  _vec->Update(a_in,*v->_vec, 1.);
}



template <typename T>
void EpetraVector<T>::insert (const std::vector<T>& v,
			      const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (v.size() == dof_indices.size());

  _vec->ReplaceGlobalValues(v.size(),
			    const_cast<T*>(&v[0]),
			    (int*) &dof_indices[0]);
}



template <typename T>
void EpetraVector<T>::insert (const NumericVector<T>& V,
			     const std::vector<unsigned int>& dof_indices)
{
  libmesh_assert (V.size() == dof_indices.size());

  // TODO: If V is an EpetraVector this can be optimized
  for (unsigned int i=0; i<V.size(); i++)
    this->set (dof_indices[i], V(i));
}



template <typename T>
void EpetraVector<T>::insert (const DenseVector<T>& /* V_in */,
			      const std::vector<unsigned int>& /* dof_indices */)
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   libmesh_assert (V_in.size() == dof_indices.size());

//   const EpetraVector<T>* V = dynamic_cast<const EpetraVector<T>*>(&V_in);

//   this->insert(V->_vec->getValues(), dof_indices);
}



template <typename T>
void EpetraVector<T>::scale (const T factor_in)
{
  _vec->Scale(factor_in);
}



template <typename T>
T EpetraVector<T>::dot (const NumericVector<T>& V_in) const
{
  const EpetraVector<T>* V = dynamic_cast<const EpetraVector<T>*>(&V_in);

  libmesh_assert(V);
  
  T result=0.0;

  _vec->Dot(*V->_vec, &result);

  return result;
}


template <typename T>
NumericVector<T>& 
EpetraVector<T>::operator = (const T s_in)
{
  _vec->PutScalar(s_in);
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator = (const NumericVector<T>& v_in)
{
  const EpetraVector<T>* v = dynamic_cast<const EpetraVector<T>*>(&v_in);

  libmesh_assert (v != NULL);
  
  *this = *v;
  
  return *this;
}



template <typename T>
EpetraVector<T>&
EpetraVector<T>::operator = (const EpetraVector<T>& v)
{
  (*_vec) = *v._vec;
  
  return *this;
}



template <typename T>
NumericVector<T>&
EpetraVector<T>::operator = (const std::vector<T>& v)
{
  T * values = _vec->Values();
  
  /**
   * Case 1:  The vector is the same size of
   * The global vector.  Only add the local components.
   */
  if(this->size() == v.size())
  {
    const unsigned int nl=this->local_size();
    const unsigned int fli=this->first_local_index();
    
    for(unsigned int i=0;i<nl;i++)
      values[i]=v[fli+i];
  }

  /**
   * Case 2: The vector is the same size as our local
   * piece.  Insert directly to the local piece.
   */
  else
  {
    libmesh_assert(v.size()==this->local_size());
    
    const unsigned int nl=this->local_size();
    
    for(unsigned int i=0;i<nl;i++)
      values[i]=v[i];
  }

  return *this;
}



template <typename T>
void EpetraVector<T>::localize (NumericVector<T>& v_local_in) const
{
  EpetraVector<T>* v_local = dynamic_cast<EpetraVector<T>*>(&v_local_in);

  libmesh_assert(v_local);

  Epetra_Map rootMap = Epetra_Util::Create_Root_Map( *_map, libMesh::processor_id() );
  Epetra_Import importer(rootMap, *_map);
  v_local->_vec->ReplaceMap(rootMap);
  
  v_local->_vec->Import(*_vec, importer, Insert);
}



template <typename T>
void EpetraVector<T>::localize (NumericVector<T>& /* v_local_in */,
				const std::vector<unsigned int>& /* send_list */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   EpetraVector<T>* v_local =
//   dynamic_cast<EpetraVector<T>*>(&v_local_in);

//   libmesh_assert (v_local != NULL);
//   libmesh_assert (v_local->local_size() == this->size());
//   libmesh_assert (send_list.size()     <= v_local->size());
  
//   int ierr=0;
//   const int n_sl = send_list.size();

//   IS is;
//   VecScatter scatter;

//   std::vector<int> idx(n_sl);
  
//   for (unsigned int i=0; i<n_sl; i++)
//     idx[i] = static_cast<int>(send_list[i]);
  
//   // Create the index set & scatter object
//   if (idx.empty())
//     ierr = ISCreateGeneral(libMesh::COMM_WORLD, n_sl, EPETRA_NULL, &is);
//   else
//     ierr = ISCreateGeneral(libMesh::COMM_WORLD, n_sl, &idx[0], &is);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   ierr = VecScatterCreate(_vec,          is,
// 			  v_local->_vec, is,
// 			  &scatter);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

  
//   // Perform the scatter
// #if EPETRA_VERSION_LESS_THAN(2,3,3)
	 
//   ierr = VecScatterBegin(_vec, v_local->_vec, INSERT_VALUES,
// 			 SCATTER_FORWARD, scatter);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//   ierr = VecScatterEnd  (_vec, v_local->_vec, INSERT_VALUES,
// 			 SCATTER_FORWARD, scatter);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

// #else
	 
//   // API argument order change in Epetra 2.3.3
//   ierr = VecScatterBegin(scatter, _vec, v_local->_vec,
//                          INSERT_VALUES, SCATTER_FORWARD);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//   ierr = VecScatterEnd  (scatter, _vec, v_local->_vec,
//                          INSERT_VALUES, SCATTER_FORWARD);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

// #endif
	 

//   // Clean up
//   ierr = ISDestroy (is);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//   ierr = VecScatterDestroy(scatter);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
}


template <typename T>
void EpetraVector<T>::localize (const unsigned int /* first_local_idx */,
				const unsigned int /* last_local_idx */,
				const std::vector<unsigned int>& /* send_list */)
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   // Only good for serial vectors.
//   libmesh_assert (this->size() == this->local_size());
//   libmesh_assert (last_local_idx > first_local_idx);
//   libmesh_assert (send_list.size() <= this->size());
//   libmesh_assert (last_local_idx < this->size());
  
//   const unsigned int size       = this->size();
//   const unsigned int local_size = (last_local_idx - first_local_idx + 1);
//   int ierr=0;  
  
//   // Don't bother for serial cases
//   if ((first_local_idx == 0) &&
//       (local_size == size))
//     return;
  
  
//   // Build a parallel vector, initialize it with the local
//   // parts of (*this)
//   EpetraVector<T> parallel_vec;

//   parallel_vec.init (size, local_size);


//   // Copy part of *this into the parallel_vec
//   {
//     IS is;
//     VecScatter scatter;

//     // Create idx, idx[i] = i+first_local_idx;
//     std::vector<int> idx(local_size);
//     Utility::iota (idx.begin(), idx.end(), first_local_idx);

//     // Create the index set & scatter object
//     ierr = ISCreateGeneral(libMesh::COMM_WORLD, local_size, &idx[0], &is); 
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);

//     ierr = VecScatterCreate(_vec,              is,
// 			    parallel_vec._vec, is,
// 			    &scatter);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);

//     // Perform the scatter
// #if EPETRA_VERSION_LESS_THAN(2,3,3)

//     ierr = VecScatterBegin(_vec, parallel_vec._vec, INSERT_VALUES,
// 			   SCATTER_FORWARD, scatter);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//     ierr = VecScatterEnd  (_vec, parallel_vec._vec, INSERT_VALUES,
// 			   SCATTER_FORWARD, scatter);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);

// #else
	   
//       // API argument order change in Epetra 2.3.3
//     ierr = VecScatterBegin(scatter, _vec, parallel_vec._vec,
// 			   INSERT_VALUES, SCATTER_FORWARD);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//     ierr = VecScatterEnd  (scatter, _vec, parallel_vec._vec,
// 			   INSERT_VALUES, SCATTER_FORWARD);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);
	   
// #endif

//     // Clean up
//     ierr = ISDestroy (is);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//     ierr = VecScatterDestroy(scatter);
//            CHKERRABORT(libMesh::COMM_WORLD,ierr);
//   }

//   // localize like normal
//   parallel_vec.close();
//   parallel_vec.localize (*this, send_list);
//   this->close();
}



template <typename T>
void EpetraVector<T>::localize (std::vector<T>& /* v_local */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   // This function must be run on all processors at once
//   parallel_only();

//   int ierr=0;
//   const int n = this->size();
//   const int nl = this->local_size();
//   EpetraScalar *values;

//   v_local.clear();
//   v_local.resize(n, 0.);

//   ierr = VecGetArray (_vec, &values);
// 	 CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   unsigned int ioff = first_local_index();

//   for (unsigned int i=0; i<nl; i++)
//     v_local[i+ioff] = static_cast<T>(values[i]);

//   ierr = VecRestoreArray (_vec, &values);
// 	 CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   Parallel::sum(v_local);
}



// Full specialization for Real datatypes
#ifdef USE_REAL_NUMBERS

template <>
void EpetraVector<Real>::localize_to_one (std::vector<Real>& /* v_local */,
					  const unsigned int /* pid */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   int ierr=0;
//   const int n  = size();
//   const int nl = local_size();
//   EpetraScalar *values;

  
//   v_local.resize(n);

  
//   // only one processor
//   if (n == nl)
//     {      
//       ierr = VecGetArray (_vec, &values);
// 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

//       for (unsigned int i=0; i<n; i++)
// 	v_local[i] = static_cast<Real>(values[i]);

//       ierr = VecRestoreArray (_vec, &values);
// 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   // otherwise multiple processors
//   else
//     {
//       unsigned int ioff = this->first_local_index();
//       std::vector<Real> local_values (n, 0.);
      
//       {
// 	ierr = VecGetArray (_vec, &values);
// 	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
// 	for (unsigned int i=0; i<nl; i++)
// 	  local_values[i+ioff] = static_cast<Real>(values[i]);
	
// 	ierr = VecRestoreArray (_vec, &values);
// 	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
//       }
      

//       MPI_Reduce (&local_values[0], &v_local[0], n, MPI_REAL, MPI_SUM,
// 		  pid, libMesh::COMM_WORLD);
//     }
}

#endif


// Full specialization for Complex datatypes
#ifdef USE_COMPLEX_NUMBERS

template <>
void EpetraVector<Complex>::localize_to_one (std::vector<Complex>& /* v_local */,
					     const unsigned int /* pid */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   int ierr=0;
//   const int n  = size();
//   const int nl = local_size();
//   EpetraScalar *values;

  
//   v_local.resize(n);

  
//   for (unsigned int i=0; i<n; i++)
//     v_local[i] = 0.;
  
//   // only one processor
//   if (n == nl)
//     {      
//       ierr = VecGetArray (_vec, &values);
// 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);

//       for (int i=0; i<n; i++)
// 	v_local[i] = static_cast<Complex>(values[i]);

//       ierr = VecRestoreArray (_vec, &values);
// 	     CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   // otherwise multiple processors
//   else
//     {
//       unsigned int ioff = this->first_local_index();

//       /* in here the local values are stored, acting as send buffer for MPI
//        * initialize to zero, since we collect using MPI_SUM
//        */
//       std::vector<Real> real_local_values(n, 0.);
//       std::vector<Real> imag_local_values(n, 0.);

//       {
// 	ierr = VecGetArray (_vec, &values);
// 	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
// 	// provide my local share to the real and imag buffers
// 	for (int i=0; i<nl; i++)
// 	  {
// 	    real_local_values[i+ioff] = static_cast<Complex>(values[i]).real();
// 	    imag_local_values[i+ioff] = static_cast<Complex>(values[i]).imag();
// 	  }

// 	ierr = VecRestoreArray (_vec, &values);
// 	       CHKERRABORT(libMesh::COMM_WORLD,ierr);
//       }
   
//       /* have buffers of the real and imaginary part of v_local.
//        * Once MPI_Reduce() collected all the real and imaginary
//        * parts in these std::vector<double>, the values can be 
//        * copied to v_local
//        */
//       std::vector<Real> real_v_local(n);
//       std::vector<Real> imag_v_local(n);

//       // collect entries from other proc's in real_v_local, imag_v_local
//       MPI_Reduce (&real_local_values[0], &real_v_local[0], n, 
// 		  MPI_DOUBLE, MPI_SUM,
// 		  pid, libMesh::COMM_WORLD);	

//       MPI_Reduce (&imag_local_values[0], &imag_v_local[0], n, 
// 		  MPI_DOUBLE, MPI_SUM,
// 		  pid, libMesh::COMM_WORLD);	

//       // copy real_v_local and imag_v_local to v_local
//       for (int i=0; i<n; i++)
// 	v_local[i] = Complex(real_v_local[i], imag_v_local[i]);
//     }  
}

#endif



template <typename T>
void EpetraVector<T>::print_matlab (const std::string /* name */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   libmesh_assert (this->initialized());
//   libmesh_assert (this->closed());
  
//   int ierr=0; 
//   EpetraViewer epetra_viewer;


//   ierr = EpetraViewerCreate (libMesh::COMM_WORLD,
// 			    &epetra_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   /**
//    * Create an ASCII file containing the matrix
//    * if a filename was provided.  
//    */
//   if (name != "NULL")
//     {
//       ierr = EpetraViewerASCIIOpen( libMesh::COMM_WORLD,
// 				   name.c_str(),
// 				   &epetra_viewer);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
      
//       ierr = EpetraViewerSetFormat (epetra_viewer,
// 				   EPETRA_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = VecView (_vec, epetra_viewer);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }

//   /**
//    * Otherwise the matrix will be dumped to the screen.
//    */
//   else
//     {
//       ierr = EpetraViewerSetFormat (EPETRA_VIEWER_STDOUT_WORLD,
// 				   EPETRA_VIEWER_ASCII_MATLAB);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
  
//       ierr = VecView (_vec, EPETRA_VIEWER_STDOUT_WORLD);
//              CHKERRABORT(libMesh::COMM_WORLD,ierr);
//     }


//   /**
//    * Destroy the viewer.
//    */
//   ierr = EpetraViewerDestroy (epetra_viewer);
//          CHKERRABORT(libMesh::COMM_WORLD,ierr);
}





template <typename T>
void EpetraVector<T>::create_subvector(NumericVector<T>& /* subvector */,
				       const std::vector<unsigned int>& /* rows */) const
{
  LIBMESH_THROW(libMesh::NotImplemented());

//   // Epetra data structures
//   IS parent_is, subvector_is;
//   VecScatter scatter;
//   int ierr = 0;
  
//   // Make sure the passed int subvector is really a EpetraVector
//   EpetraVector<T>* epetra_subvector = dynamic_cast<EpetraVector<T>*>(&subvector);
//   libmesh_assert(epetra_subvector != NULL);
  
//   // If the epetra_subvector is already initialized, we assume that the
//   // user has already allocated the *correct* amount of space for it.
//   // If not, we use the appropriate Epetra routines to initialize it.
//   if (!epetra_subvector->initialized())
//     {
//       // Initialize the epetra_subvector to have enough space to hold
//       // the entries which will be scattered into it.  Note: such an
//       // init() function (where we let Epetra decide the number of local
//       // entries) is not currently offered by the EpetraVector
//       // class.  Should we differentiate here between sequential and
//       // parallel vector creation based on libMesh::n_processors() ?
//       ierr = VecCreateMPI(libMesh::COMM_WORLD,
// 			  EPETRA_DECIDE,          // n_local
// 			  rows.size(),           // n_global
// 			  &(epetra_subvector->_vec)); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//       ierr = VecSetFromOptions (epetra_subvector->_vec); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//       // Mark the subvector as initialized
//       epetra_subvector->_is_initialized = true;
//     }
  
//   // Use iota to fill an array with entries [0,1,2,3,4,...rows.size()]
//   std::vector<int> idx(rows.size());
//   Utility::iota (idx.begin(), idx.end(), 0);

//   // Construct index sets
//   ierr = ISCreateGeneral(libMesh::COMM_WORLD,
// 			 rows.size(),
// 			 (int*) &rows[0],
// 			 &parent_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   ierr = ISCreateGeneral(libMesh::COMM_WORLD,
// 			 rows.size(),
// 			 (int*) &idx[0],
// 			 &subvector_is); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   // Construct the scatter object
//   ierr = VecScatterCreate(this->_vec,
// 			  parent_is,
// 			  epetra_subvector->_vec,
// 			  subvector_is,
// 			  &scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   // Actually perform the scatter
// #if EPETRA_VERSION_LESS_THAN(2,3,3)
//   ierr = VecScatterBegin(this->_vec,
// 			 epetra_subvector->_vec,
// 			 INSERT_VALUES,
// 			 SCATTER_FORWARD,
// 			 scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   ierr = VecScatterEnd(this->_vec,
// 		       epetra_subvector->_vec,
// 		       INSERT_VALUES,
// 		       SCATTER_FORWARD,
// 		       scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr);
// #else
//   // API argument order change in Epetra 2.3.3
//   ierr = VecScatterBegin(scatter,
// 			 this->_vec,
// 			 epetra_subvector->_vec,
// 			 INSERT_VALUES,
// 			 SCATTER_FORWARD); CHKERRABORT(libMesh::COMM_WORLD,ierr);

//   ierr = VecScatterEnd(scatter,
// 		       this->_vec,
// 		       epetra_subvector->_vec,
// 		       INSERT_VALUES,
// 		       SCATTER_FORWARD); CHKERRABORT(libMesh::COMM_WORLD,ierr);
// #endif
  
//   // Clean up 
//   ierr = ISDestroy(parent_is);       CHKERRABORT(libMesh::COMM_WORLD,ierr);
//   ierr = ISDestroy(subvector_is);    CHKERRABORT(libMesh::COMM_WORLD,ierr);
//   ierr = VecScatterDestroy(scatter); CHKERRABORT(libMesh::COMM_WORLD,ierr); 
}




//------------------------------------------------------------------
// Explicit instantiations
template class EpetraVector<Number>;

#endif // #ifdef HAVE_EPETRA
