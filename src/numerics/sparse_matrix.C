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
#include "libmesh/sparse_matrix.h"

// libMesh includes
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/diagonal_matrix.h"
#include "libmesh/laspack_matrix.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/eigen_sparse_matrix.h"
#include "libmesh/parallel.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/trilinos_epetra_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_solver_package.h"


// gzstream for reading compressed files as a stream
#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h"
#endif

#ifdef LIBMESH_HAVE_HDF5
// Sticking to the C API to be safest on portability
# include "hdf5.h"
#endif

// C++ includes
#include <memory>
#include <fstream>

#ifdef LIBMESH_HAVE_CXX11_REGEX
#include <regex>
#endif

// Anonymous namespace for helper functions
namespace
{
#ifdef LIBMESH_HAVE_HDF5
  void check_open(const std::string & filename,
                  hid_t id,
                  const std::string & objname)
  {
    if (id == H5I_INVALID_HID)
      libmesh_error_msg("Couldn't open " + objname + " in " + filename);

  }

  template <typename T>
  void check_hdf5(const std::string & filename, T hdf5val, const std::string & objname)
  {
    if (hdf5val < 0)
      libmesh_error_msg("HDF5 error from " + objname + " in " + filename);
  }
#endif
}



namespace libMesh
{


//------------------------------------------------------------------
// SparseMatrix Methods


// Constructor
template <typename T>
SparseMatrix<T>::SparseMatrix (const Parallel::Communicator & comm_in) :
  ParallelObject(comm_in),
  _dof_map(nullptr),
  _sp(nullptr),
  _is_initialized(false),
  _use_hash_table(false)
{}



template <typename T>
void SparseMatrix<T>::attach_dof_map (const DofMap & dof_map)
{
  _dof_map = &dof_map;
  if (!_sp)
    _sp = dof_map.get_sparsity_pattern();
}



template <typename T>
void SparseMatrix<T>::attach_sparsity_pattern (const SparsityPattern::Build & sp)
{
  _sp = &sp;
}



// default implementation is to fall back to non-blocked method
template <typename T>
void SparseMatrix<T>::add_block_matrix (const DenseMatrix<T> & dm,
                                        const std::vector<numeric_index_type> & brows,
                                        const std::vector<numeric_index_type> & bcols)
{
  libmesh_assert_equal_to (dm.m() / brows.size(), dm.n() / bcols.size());

  const numeric_index_type blocksize = cast_int<numeric_index_type>
    (dm.m() / brows.size());

  libmesh_assert_equal_to (dm.m()%blocksize, 0);
  libmesh_assert_equal_to (dm.n()%blocksize, 0);

  std::vector<numeric_index_type> rows, cols;

  rows.reserve(blocksize*brows.size());
  cols.reserve(blocksize*bcols.size());

  for (auto & row : brows)
    {
      numeric_index_type i = row * blocksize;

      for (unsigned int v=0; v<blocksize; v++)
        rows.push_back(i++);
    }

  for (auto & col : bcols)
    {
      numeric_index_type j = col * blocksize;

      for (unsigned int v=0; v<blocksize; v++)
        cols.push_back(j++);
    }

  this->add_matrix (dm, rows, cols);
}



// Full specialization of print method for Complex datatypes
template <>
void SparseMatrix<Complex>::print(std::ostream & os, const bool sparse) const
{
  // std::complex<>::operator<<() is defined, but use this form

  if (sparse)
    {
      libmesh_not_implemented();
    }

  os << "Real part:" << std::endl;
  for (auto i : make_range(this->m()))
    {
      for (auto j : make_range(this->n()))
        os << std::setw(8) << (*this)(i,j).real() << " ";
      os << std::endl;
    }

  os << std::endl << "Imaginary part:" << std::endl;
  for (auto i : make_range(this->m()))
    {
      for (auto j : make_range(this->n()))
        os << std::setw(8) << (*this)(i,j).imag() << " ";
      os << std::endl;
    }
}






// Full specialization for Real datatypes
template <typename T>
std::unique_ptr<SparseMatrix<T>>
SparseMatrix<T>::build(const Parallel::Communicator & comm,
                       const SolverPackage solver_package,
                       const MatrixBuildType matrix_build_type /* = AUTOMATIC */)
{
  // Avoid unused parameter warnings when no solver packages are enabled.
  libmesh_ignore(comm);

  if (matrix_build_type == MatrixBuildType::DIAGONAL)
    return std::make_unique<DiagonalMatrix<T>>(comm);

  // Build the appropriate vector
  switch (solver_package)
    {

#ifdef LIBMESH_HAVE_LASPACK
    case LASPACK_SOLVERS:
      return std::make_unique<LaspackMatrix<T>>(comm);
#endif


#ifdef LIBMESH_HAVE_PETSC
    case PETSC_SOLVERS:
      return std::make_unique<PetscMatrix<T>>(comm);
#endif


#ifdef LIBMESH_TRILINOS_HAVE_EPETRA
    case TRILINOS_SOLVERS:
      return std::make_unique<EpetraMatrix<T>>(comm);
#endif


#ifdef LIBMESH_HAVE_EIGEN
    case EIGEN_SOLVERS:
      return std::make_unique<EigenSparseMatrix<T>>(comm);
#endif

    default:
      libmesh_error_msg("ERROR:  Unrecognized solver package: " << solver_package);
    }
}


template <typename T>
void SparseMatrix<T>::vector_mult (NumericVector<T> & dest,
                                   const NumericVector<T> & arg) const
{
  dest.zero();
  this->vector_mult_add(dest,arg);
}



template <typename T>
void SparseMatrix<T>::vector_mult_add (NumericVector<T> & dest,
                                       const NumericVector<T> & arg) const
{
  /* This functionality is actually implemented in the \p
     NumericVector class.  */
  dest.add_vector(arg,*this);
}



template <typename T>
void SparseMatrix<T>::zero_rows (std::vector<numeric_index_type> &, T)
{
  /* This functionality isn't implemented or stubbed in every subclass yet */
  libmesh_not_implemented();
}


template <typename T>
std::size_t SparseMatrix<T>::n_nonzeros() const
{
  if (!_sp)
    return 0;
  return _sp->n_nonzeros();
}


template <typename T>
void SparseMatrix<T>::print(std::ostream & os, const bool sparse) const
{
  parallel_object_only();

  libmesh_assert (this->initialized());

  const numeric_index_type first_dof = this->row_start(),
                           end_dof   = this->row_stop();

  // We'll print the matrix from processor 0 to make sure
  // it's serialized properly
  if (this->processor_id() == 0)
    {
      libmesh_assert_equal_to (first_dof, 0);
      for (numeric_index_type i : make_range(end_dof))
        {
          if (sparse)
            {
              for (auto j : make_range(this->n()))
                {
                  T c = (*this)(i,j);
                  if (c != static_cast<T>(0.0))
                    {
                      os << i << " " << j << " " << c << std::endl;
                    }
                }
            }
          else
            {
              for (auto j : make_range(this->n()))
                os << (*this)(i,j) << " ";
              os << std::endl;
            }
        }

      std::vector<numeric_index_type> ibuf, jbuf;
      std::vector<T> cbuf;
      numeric_index_type currenti = end_dof;
      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          this->comm().receive(p, ibuf);
          this->comm().receive(p, jbuf);
          this->comm().receive(p, cbuf);
          libmesh_assert_equal_to (ibuf.size(), jbuf.size());
          libmesh_assert_equal_to (ibuf.size(), cbuf.size());

          if (ibuf.empty())
            continue;
          libmesh_assert_greater_equal (ibuf.front(), currenti);
          libmesh_assert_greater_equal (ibuf.back(), ibuf.front());

          std::size_t currentb = 0;
          for (;currenti <= ibuf.back(); ++currenti)
            {
              if (sparse)
                {
                  for (numeric_index_type j=0; j<this->n(); j++)
                    {
                      if (currentb < ibuf.size() &&
                          ibuf[currentb] == currenti &&
                          jbuf[currentb] == j)
                        {
                          os << currenti << " " << j << " " << cbuf[currentb] << std::endl;
                          currentb++;
                        }
                    }
                }
              else
                {
                  for (auto j : make_range(this->n()))
                    {
                      if (currentb < ibuf.size() &&
                          ibuf[currentb] == currenti &&
                          jbuf[currentb] == j)
                        {
                          os << cbuf[currentb] << " ";
                          currentb++;
                        }
                      else
                        os << static_cast<T>(0.0) << " ";
                    }
                  os << std::endl;
                }
            }
        }
      if (!sparse)
        {
          for (; currenti != this->m(); ++currenti)
            {
              for (numeric_index_type j=0; j<this->n(); j++)
                os << static_cast<T>(0.0) << " ";
              os << std::endl;
            }
        }
    }
  else
    {
      std::vector<numeric_index_type> ibuf, jbuf;
      std::vector<T> cbuf;

      // We'll assume each processor has access to entire
      // matrix rows, so (*this)(i,j) is valid if i is a local index.
      for (numeric_index_type i : make_range(first_dof, end_dof))
        {
          for (auto j : make_range(this->n()))
            {
              T c = (*this)(i,j);
              if (c != static_cast<T>(0.0))
                {
                  ibuf.push_back(i);
                  jbuf.push_back(j);
                  cbuf.push_back(c);
                }
            }
        }
      this->comm().send(0,ibuf);
      this->comm().send(0,jbuf);
      this->comm().send(0,cbuf);
    }
}


template <typename T>
void SparseMatrix<T>::print_matlab(const std::string & name) const
{
  parallel_object_only();

  libmesh_assert (this->initialized());

  const numeric_index_type first_dof = this->row_start(),
                           end_dof   = this->row_stop();

  // We'll print the matrix from processor 0 to make sure
  // it's serialized properly
  if (this->processor_id() == 0)
    {
      std::unique_ptr<std::ofstream> file;

      if (name != "")
        file = std::make_unique<std::ofstream>(name.c_str());

      std::ostream & os = (name == "") ? libMesh::out : *file;

      std::size_t sparsity_nonzeros = this->n_nonzeros();

      std::size_t real_nonzeros = 0;

      libmesh_assert_equal_to(first_dof, 0);
      for (numeric_index_type i : make_range(end_dof))
        {
          for (auto j : make_range(this->n()))
            {
              T c = (*this)(i,j);
              if (c != static_cast<T>(0.0))
                ++real_nonzeros;
            }
        }


      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          std::size_t nonzeros_on_p = 0;
          this->comm().receive(p, nonzeros_on_p);
          real_nonzeros += nonzeros_on_p;
        }

      if (sparsity_nonzeros &&
          sparsity_nonzeros != real_nonzeros)
        libmesh_warning(sparsity_nonzeros <<
                        " nonzeros allocated, but " <<
                        real_nonzeros << " used.");

      // We probably want to be more consistent than that, if our
      // sparsity is overallocated.

      // Print a header similar to PETSc's mat_view ascii_matlab
      os << "%Mat Object: () " << this->n_processors() << " MPI processes\n"
         << "%  type: " << (this->n_processors() > 1 ? "mpi" : "seq") << "aij\n"
         << "% Size = " << this->m() << ' ' << this->n() << '\n'
         << "% Nonzeros = " << real_nonzeros << '\n'
         << "zzz = zeros(" << real_nonzeros << ",3);\n"
         << "zzz = [\n";

      for (numeric_index_type i : make_range(end_dof))
        {
          // FIXME - we need a base class way to iterate over a
          // SparseMatrix row.
          for (auto j : make_range(this->n()))
            {
              T c = (*this)(i,j);
              if (c != static_cast<T>(0.0))
                {
                  // Convert from 0-based to 1-based indexing
                  os << (i+1) << ' ' << (j+1) << "  " << c << '\n';
                }
            }
        }

      std::vector<numeric_index_type> ibuf, jbuf;
      std::vector<T> cbuf;
      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          this->comm().receive(p, ibuf);
          this->comm().receive(p, jbuf);
          this->comm().receive(p, cbuf);
          libmesh_assert_equal_to (ibuf.size(), jbuf.size());
          libmesh_assert_equal_to (ibuf.size(), cbuf.size());

          for (auto n : index_range(ibuf))
            os << ibuf[n] << ' ' << jbuf[n] << "  " << cbuf[n] << '\n';
        }

      os << "];\n" << "Mat_sparse = spconvert(zzz);" << std::endl;
    }
  else
    {
      std::vector<numeric_index_type> ibuf, jbuf;
      std::vector<T> cbuf;
      std::size_t my_nonzeros = 0;

      // We'll assume each processor has access to entire
      // matrix rows, so (*this)(i,j) is valid if i is a local index.
      for (numeric_index_type i : make_range(first_dof, end_dof))
        {
          for (auto j : make_range(this->n()))
            {
              T c = (*this)(i,j);
              if (c != static_cast<T>(0.0))
                {
                  ibuf.push_back(i);
                  jbuf.push_back(j);
                  cbuf.push_back(c);
                  ++my_nonzeros;
                }
            }
        }
      this->comm().send(0,my_nonzeros);
      this->comm().send(0,ibuf);
      this->comm().send(0,jbuf);
      this->comm().send(0,cbuf);
    }
}



template <typename T>
void SparseMatrix<T>::print_coreform_hdf5(const std::string & filename,
                                          const std::string & groupname) const
{
#ifndef LIBMESH_HAVE_HDF5
  libmesh_ignore(filename, groupname);
  libmesh_error_msg("ERROR: need HDF5 support to handle .h5 files!!!");
#else
  LOG_SCOPE("print_coreform_hdf5()", "SparseMatrix");

  const numeric_index_type first_dof = this->row_start(),
                           end_dof   = this->row_stop();

  std::vector<std::size_t> cols, row_offsets;
  std::vector<double> vals;

  if (this->processor_id() == 0)
    row_offsets.push_back(0);

  for (numeric_index_type i : make_range(first_dof, end_dof))
    {
      for (auto j : make_range(this->n()))
        {
          T c = (*this)(i,j);
          if (c != static_cast<T>(0.0))
            {
              cols.push_back(j);
              vals.push_back(c);
            }
        }
      // This is a *local* row offset; proc 0 may need to adjust later
      row_offsets.push_back(cols.size());
    }

  if (this->processor_id() == 0)
    {
      const hid_t file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC,
                                   H5P_DEFAULT, H5P_DEFAULT);
      if (file == H5I_INVALID_HID)
        libmesh_file_error(filename);

      const hid_t group =
        H5Gcreate2(file, groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT,
                   H5P_DEFAULT);
      check_open(filename, group, groupname);

      auto write_size_attribute = [&filename, &group]
        (const std::string & attribute_name, unsigned long long writeval)
      {
        const hid_t fspace = H5Screate(H5S_SCALAR);
        check_hdf5(filename, fspace, attribute_name + " fspace");

        const hid_t attr = H5Acreate2(group, attribute_name.c_str(),
                                      H5T_STD_I64LE, fspace,
                                      H5P_DEFAULT, H5P_DEFAULT);
        check_hdf5(filename, attr, attribute_name);

        // HDF5 is supposed to handle both upscaling and endianness
        // conversions here
        const herr_t errval = H5Awrite(attr, H5T_NATIVE_ULLONG, &writeval);
        check_hdf5(filename, errval, attribute_name + " write");

        H5Aclose(attr);
        H5Sclose(fspace);
      };

      write_size_attribute("num_cols", this->n());
      write_size_attribute("num_rows", this->m());

      auto write_vector = [&filename, &group]
        (const std::string & dataname, auto hdf5_file_type,
         auto hdf5_native_type, auto & datavec)
      {
        const hsize_t len[1] = {cast_int<hsize_t>(datavec.size())};

        const hid_t space = H5Screate_simple(1, len, nullptr);
        check_hdf5(filename, space, dataname + " space");

        const hid_t data =
          H5Dcreate2(group, dataname.c_str(), hdf5_file_type, space,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check_hdf5(filename, data, dataname + " data");

        const hid_t errval =
          H5Dwrite(data, hdf5_native_type, H5S_ALL, H5S_ALL,
                   H5P_DEFAULT, datavec.data());
        check_hdf5(filename, errval, dataname + " write");

        H5Dclose(data);
        H5Sclose(space);
      };

      std::vector<std::size_t> vals_sizes, first_dofs;
      this->comm().gather(0, vals.size(), vals_sizes);
      this->comm().gather(0, cast_int<std::size_t>(first_dof), first_dofs);
      first_dofs.push_back(this->m());

      this->comm().allgather(cols);
      this->comm().allgather(vals);
      this->comm().allgather(row_offsets);

      libmesh_assert_equal_to(vals.size(),
                              cols.size());
      libmesh_assert_equal_to(vals.size(),
                              std::accumulate(vals_sizes.begin(),
                                              vals_sizes.end(),
                                              std::size_t(0)));

      std::size_t extra_offset = 0;
      for (auto p : make_range(processor_id_type(1), this->n_processors()))
        {
          extra_offset += vals_sizes[p-1];
          for (auto i : make_range(first_dofs[p]+1, first_dofs[p+1]+1))
            row_offsets[i] += extra_offset;
        }

      write_vector("cols", H5T_STD_U64LE, H5T_NATIVE_ULLONG, cols);
      write_vector("row_offsets", H5T_STD_U64LE, H5T_NATIVE_ULLONG, row_offsets);
      write_vector("vals", H5T_IEEE_F64LE, H5T_NATIVE_DOUBLE, vals);

      H5Gclose(group);
      H5Fclose(file);
    }
  else
    {
      std::vector<std::size_t> dummy;
      this->comm().gather(0, vals.size(), dummy);
      this->comm().gather(0, cast_int<std::size_t>(first_dof), dummy);
      this->comm().allgather(cols);
      this->comm().allgather(vals);
      this->comm().allgather(row_offsets);
    }
#endif // LIBMESH_HAVE_HDF5
}



template <typename T>
void SparseMatrix<T>::print_petsc_binary(const std::string &) const
{
  libmesh_not_implemented_msg
    ("libMesh cannot write PETSc binary-format files from non-PETSc matrices");
}



template <typename T>
void SparseMatrix<T>::print_petsc_hdf5(const std::string &) const
{
  libmesh_not_implemented_msg
    ("libMesh cannot write PETSc HDF5-format files from non-PETSc matrices");
}



template <typename T>
void SparseMatrix<T>::read(const std::string & filename)
{
  {
    std::ifstream in (filename.c_str());
    libmesh_error_msg_if
      (!in.good(), "ERROR: cannot read file:\n\t" <<
       filename);
  }

  std::string_view basename = Utility::basename_of(filename);

  const bool gzipped_file = Utility::ends_with(filename, ".gz");

  if (gzipped_file)
    basename.remove_suffix(3);

  if (Utility::ends_with(basename, ".matlab") ||
      Utility::ends_with(basename, ".m"))
    this->read_matlab(filename);
  else if (Utility::ends_with(basename, ".petsc64"))
    {
#ifndef LIBMESH_HAVE_PETSC
      libmesh_error_msg("Cannot load PETSc matrix file " <<
                        filename << " without PETSc-enabled libMesh.");
#elif LIBMESH_DOF_ID_BYTES != 8
      libmesh_error_msg("Cannot load 64-bit PETSc matrix file " <<
                        filename << " with non-64-bit libMesh.");
#endif
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped PETSc matrices are not currently supported");
      this->read_petsc_binary(filename);
    }
  else if (Utility::ends_with(basename, ".petsc32"))
    {
#ifndef LIBMESH_HAVE_PETSC
      libmesh_error_msg("Cannot load PETSc matrix file " <<
                        filename << " without PETSc-enabled libMesh.");
#elif LIBMESH_DOF_ID_BYTES != 4
      libmesh_error_msg("Cannot load 32-bit PETSc matrix file " <<
                        filename << " with non-32-bit libMesh.");
#endif
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped PETSc matrices are not currently supported");
      this->read_petsc_binary(filename);
    }
  else if (Utility::ends_with(basename, ".h5"))
    {
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped HDF5 matrices are not currently supported");
      this->read_coreform_hdf5(filename);
    }
  else
    libmesh_error_msg(" ERROR: Unrecognized matrix file extension on: "
                      << basename
                      << "\n   I understand the following:\n\n"
                      << "     *.h5        -- CoreForm HDF5 sparse matrix format\n"
                      << "     *.matlab    -- Matlab sparse matrix format\n"
                      << "     *.matlab.gz -- Matlab sparse matrix format, gzipped\n"
                      << "     *.m         -- Matlab sparse matrix format\n"
                      << "     *.m.gz      -- Matlab sparse matrix format, gzipped\n"
                      << "     *.petsc32   -- PETSc binary format, 32-bit\n"
                      << "     *.petsc64   -- PETSc binary format, 64-bit\n"
                     );
}



template <typename T>
void SparseMatrix<T>::print(const std::string & filename) const
{
  {
    std::ofstream outstr (filename.c_str());
    libmesh_error_msg_if
      (!outstr.good(), "ERROR: cannot write to file:\n\t" <<
       filename);
  }

  std::string_view basename = Utility::basename_of(filename);

  const bool gzipped_file = Utility::ends_with(filename, ".gz");

  if (gzipped_file)
    basename.remove_suffix(3);

  if (Utility::ends_with(basename, ".matlab") ||
      Utility::ends_with(basename, ".m"))
    this->print_matlab(filename);
  else if (Utility::ends_with(basename, ".petsc64"))
    {
#ifndef LIBMESH_HAVE_PETSC
      libmesh_error_msg("Cannot load PETSc matrix file " <<
                        filename << " without PETSc-enabled libMesh.");
#elif LIBMESH_DOF_ID_BYTES != 8
      libmesh_error_msg("Cannot load 64-bit PETSc matrix file " <<
                        filename << " with non-64-bit libMesh.");
#endif
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped PETSc matrices are not currently supported");
      this->print_petsc_binary(filename);
    }
  else if (Utility::ends_with(basename, ".petsc32"))
    {
#ifndef LIBMESH_HAVE_PETSC
      libmesh_error_msg("Cannot load PETSc matrix file " <<
                        filename << " without PETSc-enabled libMesh.");
#elif LIBMESH_DOF_ID_BYTES != 4
      libmesh_error_msg("Cannot load 32-bit PETSc matrix file " <<
                        filename << " with non-32-bit libMesh.");
#endif
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped PETSc matrices are not currently supported");
      this->print_petsc_binary(filename);
    }
  else if (Utility::ends_with(basename, ".h5"))
    {
      if (gzipped_file)
        libmesh_not_implemented_msg("Gzipped HDF5 matrices are not currently supported");
      this->print_coreform_hdf5(filename);
    }
  else
    libmesh_error_msg(" ERROR: Unrecognized matrix file extension on: "
                      << basename
                      << "\n   I understand the following:\n\n"
                      << "     *.h5        -- CoreForm HDF5 sparse matrix format\n"
                      << "     *.matlab    -- Matlab sparse matrix format\n"
                      << "     *.matlab.gz -- Matlab sparse matrix format, gzipped\n"
                      << "     *.m         -- Matlab sparse matrix format\n"
                      << "     *.m.gz      -- Matlab sparse matrix format, gzipped\n"
                      << "     *.petsc32   -- PETSc binary format, 32-bit\n"
                      << "     *.petsc64   -- PETSc binary format, 64-bit\n"
                     );
}



template <typename T>
void SparseMatrix<T>::read_coreform_hdf5(const std::string & filename,
                                         const std::string & groupname)
{
#ifndef LIBMESH_HAVE_HDF5
  libmesh_ignore(filename, groupname);
  libmesh_error_msg("ERROR: need HDF5 support to handle .h5 files!!!");
#else
  LOG_SCOPE("read_coreform_hdf5()", "SparseMatrix");

  std::size_t num_rows, num_cols;

  // These are only used on pid 0, but avoid "uninitialized" warnings
  hid_t group = H5I_INVALID_HID;
  hid_t file = H5I_INVALID_HID;

  if (this->processor_id() == 0)
    {
      file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      if (file == H5I_INVALID_HID)
        libmesh_file_error(filename);

      group = H5Gopen(file, groupname.c_str(), H5P_DEFAULT);
      check_open(filename, group, groupname);

      auto read_size_attribute = [&filename, &group]
        (const std::string & attribute_name)
      {
        unsigned long long returnval = 0;

        const hid_t attr = H5Aopen(group, attribute_name.c_str(), H5P_DEFAULT);
        check_open(filename, attr, attribute_name);

        const hid_t attr_type = H5Aget_type(attr);
        check_hdf5(filename, attr_type, attribute_name + " type");

        // HDF5 will convert between the file's integer type and ours, but
        // we do expect an integer type.
        if (H5Tget_class(attr_type) != H5T_INTEGER)
          libmesh_error_msg("Non-integer type for " + attribute_name + " in " + filename);

        H5Tclose(attr_type);

        // HDF5 is supposed to handle both upscaling and endianness
        // conversions here
        const herr_t errval = H5Aread(attr, H5T_NATIVE_ULLONG, &returnval);
        check_hdf5(filename, errval, attribute_name + " read");

        H5Aclose(attr);

        return returnval;
      };

      num_cols = read_size_attribute("num_cols");
      num_rows = read_size_attribute("num_rows");

      this->comm().broadcast(num_cols);
      this->comm().broadcast(num_rows);
    }
  else
    {
      this->comm().broadcast(num_cols);
      this->comm().broadcast(num_rows);
    }

  numeric_index_type new_row_start, new_row_stop,
                     new_col_start, new_col_stop;

  // If we need to reinit, we need to determine which rows+columns
  // each processor is in charge of.
  std::vector<numeric_index_type> new_row_starts, new_row_stops,
                                  new_col_starts, new_col_stops;

  if (this->initialized() &&
      num_cols == this->n() &&
      num_rows == this->m())
    {
      new_row_start = this->row_start(),
      new_row_stop  = this->row_stop();

      new_col_start = this->col_start(),
      new_col_stop  = this->col_stop();
    }
  else
    {
      // Determine which rows/columns each processor will be in charge of
      new_row_start =     this->processor_id() * num_rows / this->n_processors(),
      new_row_stop  = (this->processor_id()+1) * num_rows / this->n_processors();

      new_col_start =     this->processor_id() * num_cols / this->n_processors(),
      new_col_stop  = (this->processor_id()+1) * num_cols / this->n_processors();
    }

  this->comm().gather(0, new_row_start, new_row_starts);
  this->comm().gather(0, new_row_stop, new_row_stops);
  this->comm().gather(0, new_col_start, new_col_starts);
  this->comm().gather(0, new_col_stop, new_col_stops);

  numeric_index_type on_diagonal_nonzeros = 0,
                     off_diagonal_nonzeros = 0;

  std::vector<std::size_t> cols, row_offsets;
  std::vector<double> vals;

  if (this->processor_id() == 0)
    {
      auto read_vector = [&filename, &group]
        (const std::string & dataname, auto hdf5_class,
         auto hdf5_type, auto & datavec)
      {
        const hid_t data = H5Dopen1(group, dataname.c_str());
        check_open(filename, data, dataname.c_str());

        const hid_t data_type = H5Dget_type(data);
        check_hdf5(filename, data_type, dataname + " type");

        // HDF5 will convert between the file's integer type and ours, but
        // we do expect an integer type.
        if (H5Tget_class(data_type) != hdf5_class)
          libmesh_error_msg("Unexpected type for " + dataname + " in " + filename);

        H5Tclose(data_type);

        const hid_t dataspace = H5Dget_space(data);
        check_hdf5(filename, dataspace, dataname + " space");

        int ndims = H5Sget_simple_extent_ndims(dataspace);
        if (ndims != 1)
          libmesh_error_msg("Non-vector space for " + dataname + " in " + filename);

        hsize_t len, maxlen;
        herr_t errval = H5Sget_simple_extent_dims(dataspace, &len, &maxlen);
        check_hdf5(filename, errval, dataname + " dims");

        datavec.resize(len);

        errval = H5Dread(data, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, datavec.data());
        check_hdf5(filename, errval, dataname + " read");

        H5Dclose(data);
      };

      read_vector("cols", H5T_INTEGER, H5T_NATIVE_ULLONG, cols);
      read_vector("row_offsets", H5T_INTEGER, H5T_NATIVE_ULLONG, row_offsets);
      read_vector("vals", H5T_FLOAT, H5T_NATIVE_DOUBLE, vals);

      if (cols.size() != vals.size())
        libmesh_error_msg("Inconsistent cols/vals sizes in " + filename);

      if (row_offsets.size() != num_rows + 1)
        libmesh_error_msg("Inconsistent row_offsets size in " + filename);

      // Data for the row we're working on
      numeric_index_type current_row = 0;
      processor_id_type current_proc = 0;
      numeric_index_type current_on_diagonal_nonzeros = 0;
      numeric_index_type current_off_diagonal_nonzeros = 0;
      if (row_offsets[0] != 0)
        libmesh_error_msg("Unexpected row_offsets[0] in " + filename);

      for (auto i : index_range(cols))
        {
          while (row_offsets[current_row+1] <= i)
            {
              ++current_row;
              if (row_offsets[current_row] < row_offsets[current_row-1])
                libmesh_error_msg("Non-monotonic row_offsets in " + filename);
              current_on_diagonal_nonzeros = 0;
              current_off_diagonal_nonzeros = 0;
            }

          while (current_row >= new_row_stops[current_proc])
            ++current_proc;

          // 0-based indexing in file
          if (cols[i] >= new_col_starts[current_proc] &&
              cols[i] < new_col_stops[current_proc])
            {
              ++current_on_diagonal_nonzeros;
              on_diagonal_nonzeros =
                std::max(on_diagonal_nonzeros,
                         current_on_diagonal_nonzeros);
            }
          else
            {
              ++current_off_diagonal_nonzeros;
              off_diagonal_nonzeros =
                std::max(off_diagonal_nonzeros,
                         current_off_diagonal_nonzeros);
            }
        }
    }

  this->comm().broadcast(on_diagonal_nonzeros);
  this->comm().broadcast(off_diagonal_nonzeros);

  this->init(num_rows, num_cols,
             new_row_stop-new_row_start,
             new_col_stop-new_col_start,
             on_diagonal_nonzeros,
             off_diagonal_nonzeros);

  // Set the matrix values last.
  if (this->processor_id() == 0)
    {
      numeric_index_type current_row = 0;
      for (auto i : index_range(cols))
        {
          while (row_offsets[current_row+1] <= i)
            {
              ++current_row;
              libmesh_assert_greater_equal (row_offsets[current_row],
                                            row_offsets[current_row-1]);
            }
          this->set(current_row, cols[i], vals[i]);
        }

      H5Gclose(group);
      H5Fclose(file);
    }

  this->close();
#endif // LIBMESH_HAVE_HDF5
}



template <typename T>
void SparseMatrix<T>::read_matlab(const std::string & filename)
{
  LOG_SCOPE("read_matlab()", "SparseMatrix");

#ifndef LIBMESH_HAVE_CXX11_REGEX
  libmesh_not_implemented();  // What is your compiler?!?  Email us!
  libmesh_ignore(filename);
#else
  parallel_object_only();

  const bool gzipped_file = Utility::ends_with(filename, ".gz");

  // The sizes we get from the file
  std::size_t m = 0,
              n = 0;

  // If we don't already have this size, we'll need to reinit, and
  // determine which rows+columns each processor is in charge of.
  std::vector<numeric_index_type> new_row_starts, new_row_stops,
                                  new_col_starts, new_col_stops;

  numeric_index_type new_row_start, new_row_stop,
                     new_col_start, new_col_stop;

  // We'll read through the file three times: once to get a reliable
  // value for the matrix size (so we can divvy it up among
  // processors), then again to get the sparsity to send to each
  // processor, then a final time to get the entries to send to each
  // processor.
  //
  // We'll use an istream here; it might be an ifstream if we're
  // opening a raw ASCII file or a gzstream if we're opening a
  // compressed one.
  std::unique_ptr<std::istream> file;

  // We'll need a temporary structure to cache matrix entries, because
  // we need to read through the whole file before we know the size
  // and sparsity structure with which we can init().
  //
  // Reading through the file three times via `seekg` doesn't work
  // with our gzstream wrapper, and seems to take three times as long
  // even with a plain ifstream.  What happened to disk caching!?
  std::vector<std::tuple<numeric_index_type, numeric_index_type, T>> entries;

  // First read through the file, saving size and entry data
  {
  // We'll read the matrix on processor 0 rather than try to juggle
  // parallel I/O.
  if (this->processor_id() == 0)
    {
      // We'll be using regular expressions to make ourselves slightly
      // more robust to formatting.
      const std::regex start_regex // assignment like "zzz = ["
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

      // If we have a matrix with all-zero trailing rows, the only
      // way to get the size is if it ended up in a comment
      const std::regex size_regex // comment like "% size = 8 8"
        ("%\\s*[Ss][Ii][Zz][Ee]\\s*=\\s*(\\d+)\\s+(\\d+)");
      const std::string whitespace = " \t";

      bool have_started = false;
      bool have_ended = false;
      std::size_t largest_i_seen = 0, largest_j_seen = 0;

      // Data for the row we're working on
      // Use 1-based indexing for current_row, as in the file
      std::size_t current_row = 1;

      for (std::string line; std::getline(*file, line);)
        {
          std::smatch sm;

          // First, try to match an entry.  This is the most common
          // case so we won't rely on slow std::regex for it.
          // stringstream is at least an improvement over that.

          // Look for row/col/val like "1 1 -2.0e-4"

          std::istringstream l(line);

          std::size_t i, j;
          T value;

          l >> i >> j >> value;

          if (!l.fail())
            {
              libmesh_error_msg_if
                (!have_started, "Confused by premature entries in matrix file " << filename);

              entries.emplace_back(cast_int<numeric_index_type>(i),
                                   cast_int<numeric_index_type>(j),
                                   value);

              libmesh_error_msg_if
                (!i || !j, "Expected 1-based indexing in matrix file "
                 << filename);

              current_row = std::max(current_row, i);

              libmesh_error_msg_if
                (i < current_row,
                 "Can't handle out-of-order entries in matrix file "
                 << filename);

              largest_i_seen = std::max(i, largest_i_seen);
              largest_j_seen = std::max(j, largest_j_seen);
            }

          else if (std::regex_search(line, sm, size_regex))
            {
              const std::string msize = sm[1];
              const std::string nsize = sm[2];
              m = std::stoull(msize);
              n = std::stoull(nsize);
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
        (!have_started, "Confused by missing assignment beginning in matrix file " << filename);

      libmesh_error_msg_if
        (!have_ended, "Confused by missing assignment ending in matrix file " << filename);

      libmesh_error_msg_if
        (m > largest_i_seen, "Confused by missing final row(s) in matrix file " << filename);

      libmesh_error_msg_if
        (m > 0 && m < largest_i_seen, "Confused by extra final row(s) in matrix file " << filename);

      if (!m)
        m = largest_i_seen;

      libmesh_error_msg_if
        (n > largest_j_seen, "Confused by missing final column(s) in matrix file " << filename);

      libmesh_error_msg_if
        (n > 0 && n < largest_j_seen, "Confused by extra final column(s) in matrix file " << filename);

      if (!n)
        n = largest_j_seen;

      this->comm().broadcast(m);
      this->comm().broadcast(n);
    }
  else
    {
      this->comm().broadcast(m);
      this->comm().broadcast(n);
    }

  if (this->initialized() &&
      m == this->m() &&
      n == this->n())
    {
      new_row_start = this->row_start(),
      new_row_stop  = this->row_stop();

      new_col_start = this->col_start(),
      new_col_stop  = this->col_stop();
    }
  else
    {
      // Determine which rows/columns each processor will be in charge of
      new_row_start =     this->processor_id() * m / this->n_processors(),
      new_row_stop  = (this->processor_id()+1) * m / this->n_processors();

      new_col_start =     this->processor_id() * n / this->n_processors(),
      new_col_stop  = (this->processor_id()+1) * n / this->n_processors();
    }

  this->comm().gather(0, new_row_start, new_row_starts);
  this->comm().gather(0, new_row_stop, new_row_stops);
  this->comm().gather(0, new_col_start, new_col_starts);
  this->comm().gather(0, new_col_stop, new_col_stops);

  } // Done reading entry data and broadcasting matrix size

  // Calculate the matrix sparsity and initialize it second
  {
  // Deduce the sparsity pattern, or at least the maximum number of
  // on- and off- diagonal non-zeros per row.
  numeric_index_type  on_diagonal_nonzeros =0,
                     off_diagonal_nonzeros =0;

  if (this->processor_id() == 0)
    {
      // Data for the row we're working on
      // Use 1-based indexing for current_row, as in the file
      numeric_index_type current_row = 1;
      processor_id_type current_proc = 0;
      numeric_index_type current_on_diagonal_nonzeros = 0;
      numeric_index_type current_off_diagonal_nonzeros = 0;

      for (auto [i, j, value] : entries)
        {
          if (i > current_row)
            {
              current_row = i;
              // +1 for 1-based indexing in file
              while (current_row >= (new_row_stops[current_proc]+1))
                ++current_proc;
              current_on_diagonal_nonzeros = 0;
              current_off_diagonal_nonzeros = 0;
            }

          // +1 for 1-based indexing in file
          if (j >= (new_col_starts[current_proc]+1) &&
              j < (new_col_stops[current_proc]+1))
            {
              ++current_on_diagonal_nonzeros;
              on_diagonal_nonzeros =
                std::max(on_diagonal_nonzeros,
                         current_on_diagonal_nonzeros);
            }
          else
            {
              ++current_off_diagonal_nonzeros;
              off_diagonal_nonzeros =
                std::max(off_diagonal_nonzeros,
                         current_off_diagonal_nonzeros);
            }
        }
    }

  this->comm().broadcast(on_diagonal_nonzeros);
  this->comm().broadcast(off_diagonal_nonzeros);

  this->init(m, n,
             new_row_stop-new_row_start,
             new_col_stop-new_col_start,
             on_diagonal_nonzeros,
             off_diagonal_nonzeros);
  }

  // Set the matrix values last.
  // Convert from 1-based to 0-based indexing
  if (this->processor_id() == 0)
    for (auto [i, j, value] : entries)
      this->set(i-1, j-1, value);

  this->close();
#endif
}



template <typename T>
void SparseMatrix<T>::read_petsc_binary(const std::string &)
{
  libmesh_not_implemented_msg
    ("libMesh cannot read PETSc binary-format files into non-PETSc matrices");
}



template <typename T>
void SparseMatrix<T>::read_petsc_hdf5(const std::string &)
{
  libmesh_not_implemented_msg
    ("libMesh cannot read PETSc HDF5-format files into non-PETSc matrices");
}



template <typename T>
void SparseMatrix<T>::scale(const T scale)
{
  libmesh_assert(this->closed());

  for (const auto i : make_range(this->row_start(), this->row_stop()))
    for (const auto j : make_range(this->col_start(), this->col_stop()))
      this->set(i, j, (*this)(i, j) * scale);
}



template <typename T>
Real SparseMatrix<T>::l1_norm_diff(const SparseMatrix<T> & other_mat) const
{
  auto diff_mat = this->clone();
  diff_mat->add(-1.0, other_mat);
  return diff_mat->l1_norm();
}

//------------------------------------------------------------------
// Explicit instantiations
template class LIBMESH_EXPORT SparseMatrix<Number>;

} // namespace libMesh
