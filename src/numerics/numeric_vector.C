// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local Includes
#include "libmesh/numeric_vector.h"
#include "libmesh/distributed_vector.h"
#include "libmesh/laspack_vector.h"
#include "libmesh/eigen_sparse_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/trilinos_epetra_vector.h"
#include "libmesh/shell_matrix.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/int_range.h"

// gzstream for reading/writing compressed files as a stream
#ifdef LIBMESH_HAVE_GZSTREAM
# include "libmesh/ignore_warnings.h" // shadowing in gzstream.h
# include "gzstream.h"
# include "libmesh/restore_warnings.h"
#endif

// C++ includes
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath> // for std::abs
#include <sstream>
#include <limits>
#include <memory>

#ifdef LIBMESH_HAVE_CXX11_REGEX
#include <regex>
#endif


namespace libMesh
{



//------------------------------------------------------------------
// NumericVector methods

// Full specialization for Real datatypes
template <typename T>
std::unique_ptr<NumericVector<T>>
NumericVector<T>::build(const Parallel::Communicator & comm,
                        SolverPackage solver_package,
                        ParallelType parallel_type)
{
  // Build the appropriate vector
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      return std::make_unique<LaspackVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      return std::make_unique<PetscVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_TRILINOS_HAVE_EPETRA
    case TRILINOS_SOLVERS:
      return std::make_unique<EpetraVector<T>>(comm, parallel_type);
#endif

#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      return std::make_unique<EigenSparseVector<T>>(comm, parallel_type);
#endif

    default:
      return std::make_unique<DistributedVector<T>>(comm, parallel_type);
    }
}



template <typename T>
void NumericVector<T>::set_type(ParallelType t)
{
  // Check for no-op
  if (_type == t)
    return;

  // If the NumericVector is not yet initialized, then it is generally
  // safe to change the ParallelType, with minor restrictions.
  if (!this->initialized())
    {
      // If ghosted vectors are not enabled and the user requested a
      // GHOSTED vector, fall back on SERIAL.
#ifndef LIBMESH_ENABLE_GHOSTED
      if (t == GHOSTED)
        {
          _type = SERIAL;
          return;
        }
#endif

      _type = t;
      return;
    }

  // If we made it here, then the NumericVector was already
  // initialized and we don't currently allow the ParallelType to be
  // changed, although this could potentially be added later.
  libmesh_not_implemented();
}

template <typename T>
void NumericVector<T>::insert (const T * v,
                               const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert (v);

  for (auto i : index_range(dof_indices))
    this->set (dof_indices[i], v[i]);
}



template <typename T>
void NumericVector<T>::insert (const NumericVector<T> & V,
                               const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert_equal_to (V.size(), dof_indices.size());
  libmesh_assert (V.readable());

  for (auto i : index_range(dof_indices))
    this->set (dof_indices[i], V(i));
}



template <typename T>
int NumericVector<T>::compare (const NumericVector<T> & other_vector,
                               const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index())
  {
    if (std::abs((*this)(i) - other_vector(i)) > threshold)
      first_different_i = i;
    else
      i++;
  }

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::local_relative_compare (const NumericVector<T> & other_vector,
                                              const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  do
    {
      if (std::abs((*this)(i) - other_vector(i)) > threshold *
          std::max(std::abs((*this)(i)), std::abs(other_vector(i))))
        first_different_i = i;
      else
        i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}


template <typename T>
int NumericVector<T>::global_relative_compare (const NumericVector<T> & other_vector,
                                               const Real threshold) const
{
  libmesh_assert(this->compatible(other_vector));

  int first_different_i = std::numeric_limits<int>::max();
  numeric_index_type i = first_local_index();

  const Real my_norm = this->linfty_norm();
  const Real other_norm = other_vector.linfty_norm();
  const Real abs_threshold = std::max(my_norm, other_norm) * threshold;

  do
    {
      if (std::abs((*this)(i) - other_vector(i) ) > abs_threshold)
        first_different_i = i;
      else
        i++;
    }
  while (first_different_i==std::numeric_limits<int>::max()
         && i<last_local_index());

  // Find the correct first differing index in parallel
  this->comm().min(first_different_i);

  if (first_different_i == std::numeric_limits<int>::max())
    return -1;

  return first_different_i;
}



template <typename T>
void NumericVector<T>::print_matlab(const std::string & filename) const
{
  parallel_object_only();

  libmesh_assert (this->initialized());

  const numeric_index_type first_dof = this->first_local_index(),
                           end_dof   = this->last_local_index();

  // We'll print the vector from processor 0 to make sure
  // it's serialized properly
  if (this->processor_id() == 0)
    {
      std::unique_ptr<std::ofstream> file;

      if (filename != "")
        file = std::make_unique<std::ofstream>(filename.c_str());

      std::ostream & os = (filename == "") ? libMesh::out : *file;

      // Print a header similar to PETSC_VIEWER_ASCII_MATLAB
      os << "%Vec Object: () " << this->n_processors() << " MPI processes\n"
         << "%  type: " << (this->n_processors() > 1 ? "mpi" : "seq") << "\n"
         << "Vec = [\n";

      for (numeric_index_type i : make_range(end_dof))
        os << (*this)(i) << '\n';

      std::vector<T> cbuf;
      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          this->comm().receive(p, cbuf);

          for (auto c : cbuf)
            os << c << '\n';
        }

      os << "];\n" << std::endl;
    }
  else
    {
      std::vector<T> cbuf(end_dof - first_dof);

      for (numeric_index_type i : make_range(end_dof - first_dof))
        cbuf[i] = (*this)(first_dof+i);
      this->comm().send(0,cbuf);
    }
}



template <typename T>
void NumericVector<T>::read_matlab(const std::string & filename)
{
  LOG_SCOPE("read_matlab()", "NumericVector");

#ifndef LIBMESH_HAVE_CXX11_REGEX
  libmesh_not_implemented();  // What is your compiler?!?  Email us!
  libmesh_ignore(filename);
#else
  parallel_object_only();

  const bool gzipped_file = (filename.rfind(".gz") == filename.size() - 3);

  // If we don't already have this size, we'll need to reinit, and
  // determine which entries each processor is in charge of.
  std::vector<numeric_index_type> first_entries, end_entries;

  numeric_index_type first_entry = 0,
                     end_entry = 0,
                     n = 0;

  // We'll use an istream here; it might be an ifstream if we're
  // opening a raw ASCII file or a gzstream if we're opening a
  // compressed one.
  std::unique_ptr<std::istream> file;

  // First read through the file, saving size and entry data
  std::vector<T> entries;

  {
  // We'll read the vector on processor 0 rather than try to juggle
  // parallel I/O.
  if (this->processor_id() == 0)
    {
      // We'll be using regular expressions to make ourselves slightly
      // more robust to formatting.
      const std::regex start_regex // assignment like "Vec_0x84000002_1 = ["
        ("\\s*\\w+\\s*=\\s*\\[");
      const std::regex end_regex // end of assignment
        ("^[^%]*\\]");

      if (gzipped_file)
        {
#ifdef LIBMESH_HAVE_GZSTREAM
          auto inf = std::make_unique<igzstream>();
          libmesh_assert(inf);
          inf->open(filename.c_str(), std::ios::in);
          file = std::move(inf);
#else
          libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
        }
      else
        {
          auto inf = std::make_unique<std::ifstream>();
          libmesh_assert(inf);

          std::string new_name = Utility::unzip_file(filename);

          inf->open(new_name.c_str(), std::ios::in);
          file = std::move(inf);
        }

      const std::string whitespace = " \t";

      bool have_started = false;
      bool have_ended = false;

      for (std::string line; std::getline(*file, line);)
        {
          std::smatch sm;

          // First, try to match an entry.  This is the most common
          // case so we won't rely on slow std::regex for it.
          // stringstream is at least an improvement over that.
          std::istringstream l(line);
          T value;
          l >> value;

          if (!l.fail())
            {
              libmesh_error_msg_if
                (!have_started, "Confused by premature entries in vector file " << filename);

              entries.push_back(value);
              ++n;
            }

          else if (std::regex_search(line, start_regex))
            have_started = true;

          else if (std::regex_search(line, end_regex))
            {
              have_ended = true;
              break;
            }
        }

      libmesh_error_msg_if
        (!have_started, "Confused by missing assignment-beginning in vector file " << filename);

      libmesh_error_msg_if
        (!have_ended, "Confused by missing assignment-ending in vector file " << filename);
    }

  this->comm().broadcast(n);

  if (this->initialized() &&
      n == this->size())
    {
      first_entry = this->first_local_index(),
      end_entry = this->last_local_index();
    }
  else
    {
      // Determine which rows/columns each processor will be in charge of
      first_entry =     this->processor_id() * n / this->n_processors(),
      end_entry   = (this->processor_id()+1) * n / this->n_processors();
    }

  this->comm().gather(0, first_entry, first_entries);
  this->comm().gather(0, end_entry, end_entries);

  } // Done reading entry data and broadcasting vector size

  // If we're not already initialized compatibly with the file then
  // we'll initialize here
  bool need_init = !this->initialized() ||
                   (this->size() != n) ||
                   (this->local_size() != end_entry - first_entry);

  this->comm().max(need_init);
  
  if (need_init)
    this->init(n, end_entry - first_entry);

  // Set the vector values last.  The iota call here is inefficient
  // but it's probably better than calling a single virtual function
  // per index.
  if (this->processor_id() == 0)
    for (auto p : make_range(this->n_processors()))
      {
        const numeric_index_type first_entry_p = first_entries[p];
        const numeric_index_type n_local = end_entries[p] - first_entry_p;
        std::vector<numeric_index_type> indices(n_local);
        std::iota(indices.begin(), indices.end(), first_entry_p);
        this->insert(entries.data() + first_entries[p],
                     indices);
      }

  this->close();
#endif
}


/*
// Full specialization for float datatypes (DistributedVector wants this)

template <>
int NumericVector<float>::compare (const NumericVector<float> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}

// Full specialization for double datatypes
template <>
int NumericVector<double>::compare (const NumericVector<double> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}

#ifdef LIBMESH_DEFAULT_TRIPLE_PRECISION
// Full specialization for long double datatypes
template <>
int NumericVector<long double>::compare (const NumericVector<long double> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if (std::abs((*this)(i) - other_vector(i) ) > threshold)
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<last_local_index());

return rvalue;
}
#endif


// Full specialization for Complex datatypes
template <>
int NumericVector<Complex>::compare (const NumericVector<Complex> & other_vector,
const Real threshold) const
{
libmesh_assert (this->initialized());
libmesh_assert (other_vector.initialized());
libmesh_assert_equal_to (this->first_local_index(), other_vector.first_local_index());
libmesh_assert_equal_to (this->last_local_index(), other_vector.last_local_index());

int rvalue     = -1;
numeric_index_type i = first_local_index();

do
{
if ((std::abs((*this)(i).real() - other_vector(i).real()) > threshold) || (std::abs((*this)(i).imag() - other_vector(i).imag()) > threshold))
rvalue = i;
else
i++;
}
while (rvalue==-1 && i<this->last_local_index());

return rvalue;
}
*/


template <class T>
Real NumericVector<T>::subset_l1_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    norm += std::abs(v(index));

  this->comm().sum(norm);

  return norm;
}

template <class T>
Real NumericVector<T>::subset_l2_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    norm += TensorTools::norm_sq(v(index));

  this->comm().sum(norm);

  return std::sqrt(norm);
}

template <class T>
Real NumericVector<T>::subset_linfty_norm (const std::set<numeric_index_type> & indices) const
{
  libmesh_assert (this->readable());

  const NumericVector<T> & v = *this;

  Real norm = 0;

  for (const auto & index : indices)
    {
      Real value = std::abs(v(index));
      if (value > norm)
        norm = value;
    }

  this->comm().max(norm);

  return norm;
}



template <class T>
Real NumericVector<T>::l2_norm_diff (const NumericVector<T> & v) const
{
  libmesh_assert(this->compatible(v));

  Real norm = 0;
  for (const auto i : make_range(this->first_local_index(), this->last_local_index()))
    norm += TensorTools::norm_sq((*this)(i) - v(i));

  this->comm().sum(norm);

  return std::sqrt(norm);
}



template <class T>
Real NumericVector<T>::l1_norm_diff (const NumericVector<T> & v) const
{
  libmesh_assert(this->compatible(v));

  Real norm = 0;
  for (const auto i : make_range(this->first_local_index(), this->last_local_index()))
    norm += libMesh::l1_norm_diff((*this)(i), v(i));

  this->comm().sum(norm);

  return norm;
}



template <typename T>
void NumericVector<T>::add_vector (const T * v,
                                   const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v);

  for (auto i : index_range(dof_indices))
    this->add (dof_indices[i], v[i]);
}



template <typename T>
void NumericVector<T>::add_vector (const NumericVector<T> & v,
                                   const std::vector<numeric_index_type> & dof_indices)
{
  libmesh_assert(v.readable());

  const std::size_t n = dof_indices.size();
  libmesh_assert_equal_to(v.size(), n);
  for (numeric_index_type i=0; i != n; i++)
    this->add (dof_indices[i], v(i));
}



template <typename T>
void NumericVector<T>::add_vector (const NumericVector<T> & v,
                                   const ShellMatrix<T> & a)
{
  libmesh_assert(this->compatible(v));

  a.vector_mult_add(*this,v);
}



template <typename T>
bool NumericVector<T>::readable () const
{
  return this->initialized() && this->closed();
}


template <typename T>
bool NumericVector<T>::compatible (const NumericVector<T> & v) const
{
  return this->readable() && v.readable() &&
         this->size() == v.size() &&
         this->local_size() == v.local_size() &&
         this->first_local_index() == v.first_local_index() &&
         this->last_local_index() == v.last_local_index();
}


//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT NumericVector<Number>;

} // namespace libMesh
