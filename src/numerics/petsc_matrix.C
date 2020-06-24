// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <unistd.h> // mkstemp
#include <fstream>

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_PETSC

// Local includes
#include "libmesh/petsc_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/parallel.h"


#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE

namespace
{
using namespace libMesh;

// historic libMesh n_nz & n_oz arrays are set up for PETSc's AIJ format.
// however, when the blocksize is >1, we need to transform these into
// their BAIJ counterparts.
inline
void transform_preallocation_arrays (const PetscInt blocksize,
                                     const std::vector<numeric_index_type> & n_nz,
                                     const std::vector<numeric_index_type> & n_oz,
                                     std::vector<numeric_index_type>       & b_n_nz,
                                     std::vector<numeric_index_type>       & b_n_oz)
{
  libmesh_assert_equal_to (n_nz.size(), n_oz.size());
  libmesh_assert_equal_to (n_nz.size()%blocksize, 0);

  b_n_nz.clear(); /**/ b_n_nz.reserve(n_nz.size()/blocksize);
  b_n_oz.clear(); /**/ b_n_oz.reserve(n_oz.size()/blocksize);

  for (std::size_t nn=0, nnzs=n_nz.size(); nn<nnzs; nn += blocksize)
    {
      b_n_nz.push_back (n_nz[nn]/blocksize);
      b_n_oz.push_back (n_oz[nn]/blocksize);
    }
}
}

#endif



namespace libMesh
{


//-----------------------------------------------------------------------
// PetscMatrix members


// Constructor
template <typename T>
PetscMatrix<T>::PetscMatrix(const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(true),
  _mat_type(AIJ)
{}



// Constructor taking an existing Mat but not the responsibility
// for destroying it
template <typename T>
PetscMatrix<T>::PetscMatrix(Mat mat_in,
                            const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(false),
  _mat_type(AIJ)
{
  this->_mat = mat_in;
  this->_is_initialized = true;
}



// Destructor
template <typename T>
PetscMatrix<T>::~PetscMatrix()
{
  this->clear();
}

template <typename T>
void PetscMatrix<T>::set_matrix_type(PetscMatrixType mat_type)
{
  _mat_type = mat_type;
}

template <typename T>
void PetscMatrix<T>::init (const numeric_index_type m_in,
                           const numeric_index_type n_in,
                           const numeric_index_type m_l,
                           const numeric_index_type n_l,
                           const numeric_index_type nnz,
                           const numeric_index_type noz,
                           const numeric_index_type blocksize_in)
{
  // So compilers don't warn when !LIBMESH_ENABLE_BLOCKED_STORAGE
  libmesh_ignore(blocksize_in);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;


  PetscErrorCode ierr = 0;
  PetscInt m_global   = static_cast<PetscInt>(m_in);
  PetscInt n_global   = static_cast<PetscInt>(n_in);
  PetscInt m_local    = static_cast<PetscInt>(m_l);
  PetscInt n_local    = static_cast<PetscInt>(n_l);
  PetscInt n_nz       = static_cast<PetscInt>(nnz);
  PetscInt n_oz       = static_cast<PetscInt>(noz);

  ierr = MatCreate(this->comm().get(), &_mat);
  LIBMESH_CHKERR(ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
  LIBMESH_CHKERR(ierr);
  PetscInt blocksize  = static_cast<PetscInt>(blocksize_in);
  ierr = MatSetBlockSize(_mat,blocksize);
  LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  if (blocksize > 1)
    {
      // specified blocksize, bs>1.
      // double check sizes.
      libmesh_assert_equal_to (m_local  % blocksize, 0);
      libmesh_assert_equal_to (n_local  % blocksize, 0);
      libmesh_assert_equal_to (m_global % blocksize, 0);
      libmesh_assert_equal_to (n_global % blocksize, 0);
      libmesh_assert_equal_to (n_nz     % blocksize, 0);
      libmesh_assert_equal_to (n_oz     % blocksize, 0);

      ierr = MatSetType(_mat, MATBAIJ); // Automatically chooses seqbaij or mpibaij
      LIBMESH_CHKERR(ierr);
      ierr = MatSeqBAIJSetPreallocation(_mat, blocksize, n_nz/blocksize, NULL);
      LIBMESH_CHKERR(ierr);
      ierr = MatMPIBAIJSetPreallocation(_mat, blocksize,
                                        n_nz/blocksize, NULL,
                                        n_oz/blocksize, NULL);
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {
      switch (_mat_type) {
        case AIJ:
          ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
          LIBMESH_CHKERR(ierr);
          ierr = MatSeqAIJSetPreallocation(_mat, n_nz, NULL);
          LIBMESH_CHKERR(ierr);
          ierr = MatMPIAIJSetPreallocation(_mat, n_nz, NULL, n_oz, NULL);
          LIBMESH_CHKERR(ierr);
          break;

        case HYPRE:
#if !PETSC_VERSION_LESS_THAN(3,9,4) && LIBMESH_HAVE_PETSC_HYPRE
          ierr = MatSetType(_mat, MATHYPRE);
          LIBMESH_CHKERR(ierr);
          ierr = MatHYPRESetPreallocation(_mat, n_nz, NULL, n_oz, NULL);
          LIBMESH_CHKERR(ierr);
#else
          libmesh_error_msg("PETSc 3.9.4 or higher with hypre is required for MatHypre");
#endif
          break;

        default: libmesh_error_msg("Unsupported petsc matrix type");
      }
    }

  // Make it an error for PETSc to allocate new nonzero entries during assembly
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  LIBMESH_CHKERR(ierr);

  // Is prefix information available somewhere? Perhaps pass in the system name?
  ierr = MatSetOptionsPrefix(_mat, "");
  LIBMESH_CHKERR(ierr);
  ierr = MatSetFromOptions(_mat);
  LIBMESH_CHKERR(ierr);

  this->zero ();
}


template <typename T>
void PetscMatrix<T>::init (const numeric_index_type m_in,
                           const numeric_index_type n_in,
                           const numeric_index_type m_l,
                           const numeric_index_type n_l,
                           const std::vector<numeric_index_type> & n_nz,
                           const std::vector<numeric_index_type> & n_oz,
                           const numeric_index_type blocksize_in)
{
  // So compilers don't warn when !LIBMESH_ENABLE_BLOCKED_STORAGE
  libmesh_ignore(blocksize_in);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;

  // Make sure the sparsity pattern isn't empty unless the matrix is 0x0
  libmesh_assert_equal_to (n_nz.size(), m_l);
  libmesh_assert_equal_to (n_oz.size(), m_l);

  PetscErrorCode ierr = 0;
  PetscInt m_global   = static_cast<PetscInt>(m_in);
  PetscInt n_global   = static_cast<PetscInt>(n_in);
  PetscInt m_local    = static_cast<PetscInt>(m_l);
  PetscInt n_local    = static_cast<PetscInt>(n_l);

  ierr = MatCreate(this->comm().get(), &_mat);
  LIBMESH_CHKERR(ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
  LIBMESH_CHKERR(ierr);
  PetscInt blocksize  = static_cast<PetscInt>(blocksize_in);
  ierr = MatSetBlockSize(_mat,blocksize);
  LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  if (blocksize > 1)
    {
      // specified blocksize, bs>1.
      // double check sizes.
      libmesh_assert_equal_to (m_local  % blocksize, 0);
      libmesh_assert_equal_to (n_local  % blocksize, 0);
      libmesh_assert_equal_to (m_global % blocksize, 0);
      libmesh_assert_equal_to (n_global % blocksize, 0);

      ierr = MatSetType(_mat, MATBAIJ); // Automatically chooses seqbaij or mpibaij
      LIBMESH_CHKERR(ierr);

      // transform the per-entry n_nz and n_oz arrays into their block counterparts.
      std::vector<numeric_index_type> b_n_nz, b_n_oz;

      transform_preallocation_arrays (blocksize,
                                      n_nz, n_oz,
                                      b_n_nz, b_n_oz);

      ierr = MatSeqBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? nullptr : b_n_nz.data()));
      LIBMESH_CHKERR(ierr);

      ierr = MatMPIBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? nullptr : b_n_nz.data()),
                                         0,
                                         numeric_petsc_cast(b_n_oz.empty() ? nullptr : b_n_oz.data()));
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {
      switch (_mat_type) {
        case AIJ:
          ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
          LIBMESH_CHKERR(ierr);
          ierr = MatSeqAIJSetPreallocation (_mat,
                                            0,
                                            numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()));
          LIBMESH_CHKERR(ierr);
          ierr = MatMPIAIJSetPreallocation (_mat,
                                            0,
                                            numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()),
                                            0,
                                            numeric_petsc_cast(n_oz.empty() ? nullptr : n_oz.data()));
          break;

        case HYPRE:
#if !PETSC_VERSION_LESS_THAN(3,9,4) && LIBMESH_HAVE_PETSC_HYPRE
          ierr = MatSetType(_mat, MATHYPRE);
          LIBMESH_CHKERR(ierr);
          ierr = MatHYPRESetPreallocation (_mat,
                                           0,
                                           numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()),
                                           0,
                                           numeric_petsc_cast(n_oz.empty() ? nullptr : n_oz.data()));
          LIBMESH_CHKERR(ierr);
#else
          libmesh_error_msg("PETSc 3.9.4 or higher with hypre is required for MatHypre");
#endif
          break;

        default: libmesh_error_msg("Unsupported petsc matrix type");
      }

    }

  // Make it an error for PETSc to allocate new nonzero entries during assembly
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  LIBMESH_CHKERR(ierr);

  // Is prefix information available somewhere? Perhaps pass in the system name?
  ierr = MatSetOptionsPrefix(_mat, "");
  LIBMESH_CHKERR(ierr);
  ierr = MatSetFromOptions(_mat);
  LIBMESH_CHKERR(ierr);


  this->zero();
}


template <typename T>
void PetscMatrix<T>::init (const ParallelType)
{
  libmesh_assert(this->_dof_map);

  // Clear initialized matrices
  if (this->initialized())
    this->clear();

  this->_is_initialized = true;


  const numeric_index_type my_m = this->_dof_map->n_dofs();
  const numeric_index_type my_n = my_m;
  const numeric_index_type n_l  = this->_dof_map->n_dofs_on_processor(this->processor_id());
  const numeric_index_type m_l  = n_l;


  const std::vector<numeric_index_type> & n_nz = this->_dof_map->get_n_nz();
  const std::vector<numeric_index_type> & n_oz = this->_dof_map->get_n_oz();

  // Make sure the sparsity pattern isn't empty unless the matrix is 0x0
  libmesh_assert_equal_to (n_nz.size(), m_l);
  libmesh_assert_equal_to (n_oz.size(), m_l);

  PetscErrorCode ierr = 0;
  PetscInt m_global   = static_cast<PetscInt>(my_m);
  PetscInt n_global   = static_cast<PetscInt>(my_n);
  PetscInt m_local    = static_cast<PetscInt>(m_l);
  PetscInt n_local    = static_cast<PetscInt>(n_l);

  ierr = MatCreate(this->comm().get(), &_mat);
  LIBMESH_CHKERR(ierr);
  ierr = MatSetSizes(_mat, m_local, n_local, m_global, n_global);
  LIBMESH_CHKERR(ierr);
  PetscInt blocksize  = static_cast<PetscInt>(this->_dof_map->block_size());
  ierr = MatSetBlockSize(_mat,blocksize);
  LIBMESH_CHKERR(ierr);

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  if (blocksize > 1)
    {
      // specified blocksize, bs>1.
      // double check sizes.
      libmesh_assert_equal_to (m_local  % blocksize, 0);
      libmesh_assert_equal_to (n_local  % blocksize, 0);
      libmesh_assert_equal_to (m_global % blocksize, 0);
      libmesh_assert_equal_to (n_global % blocksize, 0);

      ierr = MatSetType(_mat, MATBAIJ); // Automatically chooses seqbaij or mpibaij
      LIBMESH_CHKERR(ierr);

      // transform the per-entry n_nz and n_oz arrays into their block counterparts.
      std::vector<numeric_index_type> b_n_nz, b_n_oz;

      transform_preallocation_arrays (blocksize,
                                      n_nz, n_oz,
                                      b_n_nz, b_n_oz);

      ierr = MatSeqBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? nullptr : b_n_nz.data()));
      LIBMESH_CHKERR(ierr);

      ierr = MatMPIBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? nullptr : b_n_nz.data()),
                                         0,
                                         numeric_petsc_cast(b_n_oz.empty() ? nullptr : b_n_oz.data()));
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {
      switch (_mat_type) {
        case AIJ:
          ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
          LIBMESH_CHKERR(ierr);
          ierr = MatSeqAIJSetPreallocation (_mat,
                                            0,
                                            numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()));
          LIBMESH_CHKERR(ierr);
          ierr = MatMPIAIJSetPreallocation (_mat,
                                            0,
                                            numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()),
                                            0,
                                            numeric_petsc_cast(n_oz.empty() ? nullptr : n_oz.data()));
          break;

        case HYPRE:
#if !PETSC_VERSION_LESS_THAN(3,9,4) && LIBMESH_HAVE_PETSC_HYPRE
          ierr = MatSetType(_mat, MATHYPRE);
          LIBMESH_CHKERR(ierr);
          ierr = MatHYPRESetPreallocation (_mat,
                                           0,
                                           numeric_petsc_cast(n_nz.empty() ? nullptr : n_nz.data()),
                                           0,
                                           numeric_petsc_cast(n_oz.empty() ? nullptr : n_oz.data()));
          LIBMESH_CHKERR(ierr);
#else
          libmesh_error_msg("PETSc 3.9.4 or higher with hypre is required for MatHypre");
#endif
          break;

        default: libmesh_error_msg("Unsupported petsc matrix type");
      }
    }

  // Make it an error for PETSc to allocate new nonzero entries during assembly
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  LIBMESH_CHKERR(ierr);

  // Is prefix information available somewhere? Perhaps pass in the system name?
  ierr = MatSetOptionsPrefix(_mat, "");
  LIBMESH_CHKERR(ierr);
  ierr = MatSetFromOptions(_mat);
  LIBMESH_CHKERR(ierr);

  this->zero();
}


template <typename T>
void PetscMatrix<T>::update_preallocation_and_zero ()
{
  libmesh_not_implemented();
}

template <typename T>
void PetscMatrix<T>::reset_preallocation()
{
#if !PETSC_VERSION_LESS_THAN(3,9,0)
  libmesh_assert (this->initialized());

  auto ierr = MatResetPreallocation(_mat);
  LIBMESH_CHKERR(ierr);
#else
  libmesh_warning("Your version of PETSc doesn't support resetting of "
                  "preallocation, so we will use your most recent sparsity "
                  "pattern. This may result in a degradation of performance\n");
#endif
}

template <typename T>
void PetscMatrix<T>::zero ()
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  PetscErrorCode ierr=0;

  PetscInt m_l, n_l;

  ierr = MatGetLocalSize(_mat,&m_l,&n_l);
  LIBMESH_CHKERR(ierr);

  if (n_l)
    {
      ierr = MatZeroEntries(_mat);
      LIBMESH_CHKERR(ierr);
    }
}

template <typename T>
void PetscMatrix<T>::zero_rows (std::vector<numeric_index_type> & rows, T diag_value)
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  PetscErrorCode ierr=0;

  // As of petsc-dev at the time of 3.1.0, MatZeroRows now takes two additional
  // optional arguments.  The optional arguments (x,b) can be used to specify the
  // solutions for the zeroed rows (x) and right hand side (b) to update.
  // Could be useful for setting boundary conditions...
  if (!rows.empty())
    ierr = MatZeroRows(_mat, cast_int<PetscInt>(rows.size()),
                       numeric_petsc_cast(rows.data()), PS(diag_value),
                       NULL, NULL);
  else
    ierr = MatZeroRows(_mat, 0, NULL, PS(diag_value), NULL, NULL);

  LIBMESH_CHKERR(ierr);
}

template <typename T>
std::unique_ptr<SparseMatrix<T>> PetscMatrix<T>::zero_clone () const
{
  if (!this->closed())
    libmesh_error_msg("Matrix must be closed before it can be cloned!");

  // Copy the nonzero pattern only
  Mat copy;
  PetscErrorCode ierr = MatDuplicate(_mat, MAT_DO_NOT_COPY_VALUES, &copy);
  LIBMESH_CHKERR(ierr);

  // Call wrapping PetscMatrix constructor, have it take over
  // ownership.
  auto ret = libmesh_make_unique<PetscMatrix<T>>(copy, this->comm());
  ret->set_destroy_mat_on_exit(true);

  // Work around an issue on older compilers.  We are able to simply
  // "return ret;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(ret.release());
}



template <typename T>
std::unique_ptr<SparseMatrix<T>> PetscMatrix<T>::clone () const
{
  if (!this->closed())
    libmesh_error_msg("Matrix must be closed before it can be cloned!");

  // Copy the nonzero pattern and numerical values
  Mat copy;
  PetscErrorCode ierr = MatDuplicate(_mat, MAT_COPY_VALUES, &copy);
  LIBMESH_CHKERR(ierr);

  // Call wrapping PetscMatrix constructor, have it take over
  // ownership.
  auto ret = libmesh_make_unique<PetscMatrix<T>>(copy, this->comm());
  ret->set_destroy_mat_on_exit(true);

  // Work around an issue on older compilers.  We are able to simply
  // "return ret;" on newer compilers
  return std::unique_ptr<SparseMatrix<T>>(ret.release());
}

template <typename T>
void PetscMatrix<T>::clear ()
{
  PetscErrorCode ierr=0;

  if ((this->initialized()) && (this->_destroy_mat_on_exit))
    {
      semiparallel_only();

      ierr = MatDestroy (&_mat);
      LIBMESH_CHKERR(ierr);

      this->_is_initialized = false;
    }
}



template <typename T>
Real PetscMatrix<T>::l1_norm () const
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  PetscErrorCode ierr=0;
  PetscReal petsc_value;
  Real value;

  libmesh_assert (this->closed());

  ierr = MatNorm(_mat, NORM_1, &petsc_value);
  LIBMESH_CHKERR(ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
Real PetscMatrix<T>::linfty_norm () const
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  PetscErrorCode ierr=0;
  PetscReal petsc_value;
  Real value;

  libmesh_assert (this->closed());

  ierr = MatNorm(_mat, NORM_INFINITY, &petsc_value);
  LIBMESH_CHKERR(ierr);

  value = static_cast<Real>(petsc_value);

  return value;
}



template <typename T>
void PetscMatrix<T>::print_matlab (const std::string & name) const
{
  libmesh_assert (this->initialized());

  semiparallel_only();

  if (!this->closed())
    {
      libmesh_deprecated();
      libmesh_warning("The matrix must be assembled before calling PetscMatrix::print_matlab().\n"
                      "Please update your code, as this warning will become an error in a future release.");
      const_cast<PetscMatrix<T> *>(this)->close();
    }

  PetscErrorCode ierr=0;
  PetscViewer petsc_viewer;


  ierr = PetscViewerCreate (this->comm().get(),
                            &petsc_viewer);
  LIBMESH_CHKERR(ierr);

  /**
   * Create an ASCII file containing the matrix
   * if a filename was provided.
   */
  if (name != "")
    {
      ierr = PetscViewerASCIIOpen( this->comm().get(),
                                   name.c_str(),
                                   &petsc_viewer);
      LIBMESH_CHKERR(ierr);

#if PETSC_VERSION_LESS_THAN(3,7,0)
      ierr = PetscViewerSetFormat (petsc_viewer,
                                   PETSC_VIEWER_ASCII_MATLAB);
#else
      ierr = PetscViewerPushFormat (petsc_viewer,
                                    PETSC_VIEWER_ASCII_MATLAB);
#endif

      LIBMESH_CHKERR(ierr);

      ierr = MatView (_mat, petsc_viewer);
      LIBMESH_CHKERR(ierr);
    }

  /**
   * Otherwise the matrix will be dumped to the screen.
   */
  else
    {
#if PETSC_VERSION_LESS_THAN(3,7,0)
      ierr = PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD,
                                   PETSC_VIEWER_ASCII_MATLAB);
#else
      ierr = PetscViewerPushFormat (PETSC_VIEWER_STDOUT_WORLD,
                                    PETSC_VIEWER_ASCII_MATLAB);
#endif

      LIBMESH_CHKERR(ierr);

      ierr = MatView (_mat, PETSC_VIEWER_STDOUT_WORLD);
      LIBMESH_CHKERR(ierr);
    }


  /**
   * Destroy the viewer.
   */
  ierr = PetscViewerDestroy (&petsc_viewer);
  LIBMESH_CHKERR(ierr);
}





template <typename T>
void PetscMatrix<T>::print_personal(std::ostream & os) const
{
  libmesh_assert (this->initialized());

  // Routine must be called in parallel on parallel matrices
  // and serial on serial matrices.
  semiparallel_only();

  // #ifndef NDEBUG
  //   if (os != std::cout)
  //     libMesh::err << "Warning! PETSc can only print to std::cout!" << std::endl;
  // #endif

  // Matrix must be in an assembled state to be printed
  if (!this->closed())
    {
      libmesh_deprecated();
      libmesh_warning("The matrix must be assembled before calling PetscMatrix::print_personal().\n"
                      "Please update your code, as this warning will become an error in a future release.");
      const_cast<PetscMatrix<T> *>(this)->close();
    }

  PetscErrorCode ierr=0;

  // Print to screen if ostream is stdout
  if (os.rdbuf() == std::cout.rdbuf())
    {
      ierr = MatView(_mat, PETSC_VIEWER_STDOUT_SELF);
      LIBMESH_CHKERR(ierr);
    }

  // Otherwise, print to the requested file, in a roundabout way...
  else
    {
      // We will create a temporary filename, and file, for PETSc to
      // write to.
      std::string temp_filename;

      {
        // Template for temporary filename
        char c[] = "temp_petsc_matrix.XXXXXX";

        // Generate temporary, unique filename only on processor 0.  We will
        // use this filename for PetscViewerASCIIOpen, before copying it into
        // the user's stream
        if (this->processor_id() == 0)
          {
            int fd = mkstemp(c);

            // Check to see that mkstemp did not fail.
            if (fd == -1)
              libmesh_error_msg("mkstemp failed in PetscMatrix::print_personal()");

            // mkstemp returns a file descriptor for an open file,
            // so let's close it before we hand it to PETSc!
            ::close (fd);
          }

        // Store temporary filename as string, makes it easier to broadcast
        temp_filename = c;
      }

      // Now broadcast the filename from processor 0 to all processors.
      this->comm().broadcast(temp_filename);

      // PetscViewer object for passing to MatView
      PetscViewer petsc_viewer;

      // This PETSc function only takes a string and handles the opening/closing
      // of the file internally.  Since print_personal() takes a reference to
      // an ostream, we have to do an extra step...  print_personal() should probably
      // have a version that takes a string to get rid of this problem.
      ierr = PetscViewerASCIIOpen( this->comm().get(),
                                   temp_filename.c_str(),
                                   &petsc_viewer);
      LIBMESH_CHKERR(ierr);

      // Probably don't need to set the format if it's default...
      //      ierr = PetscViewerSetFormat (petsc_viewer,
      //   PETSC_VIEWER_DEFAULT);
      //      LIBMESH_CHKERR(ierr);

      // Finally print the matrix using the viewer
      ierr = MatView (_mat, petsc_viewer);
      LIBMESH_CHKERR(ierr);

      if (this->processor_id() == 0)
        {
          // Now the inefficient bit: open temp_filename as an ostream and copy the contents
          // into the user's desired ostream.  We can't just do a direct file copy, we don't even have the filename!
          std::ifstream input_stream(temp_filename.c_str());
          os << input_stream.rdbuf();  // The "most elegant" way to copy one stream into another.
          // os.close(); // close not defined in ostream

          // Now remove the temporary file
          input_stream.close();
          std::remove(temp_filename.c_str());
        }
    }
}






template <typename T>
void PetscMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                const std::vector<numeric_index_type> & rows,
                                const std::vector<numeric_index_type> & cols)
{
  libmesh_assert (this->initialized());

  const numeric_index_type n_rows = dm.m();
  const numeric_index_type n_cols = dm.n();

  libmesh_assert_equal_to (rows.size(), n_rows);
  libmesh_assert_equal_to (cols.size(), n_cols);

  PetscErrorCode ierr=0;
  ierr = MatSetValues(_mat,
                      n_rows, numeric_petsc_cast(rows.data()),
                      n_cols, numeric_petsc_cast(cols.data()),
                      pPS(const_cast<T*>(dm.get_values().data())),
                      ADD_VALUES);
  LIBMESH_CHKERR(ierr);
}






template <typename T>
void PetscMatrix<T>::add_block_matrix(const DenseMatrix<T> & dm,
                                      const std::vector<numeric_index_type> & brows,
                                      const std::vector<numeric_index_type> & bcols)
{
  libmesh_assert (this->initialized());

  const numeric_index_type n_brows =
    cast_int<numeric_index_type>(brows.size());
  const numeric_index_type n_bcols =
    cast_int<numeric_index_type>(bcols.size());

  PetscErrorCode ierr=0;

#ifndef NDEBUG
  const numeric_index_type n_rows  =
    cast_int<numeric_index_type>(dm.m());
  const numeric_index_type n_cols  =
    cast_int<numeric_index_type>(dm.n());
  const numeric_index_type blocksize = n_rows / n_brows;

  libmesh_assert_equal_to (n_cols / n_bcols, blocksize);
  libmesh_assert_equal_to (blocksize*n_brows, n_rows);
  libmesh_assert_equal_to (blocksize*n_bcols, n_cols);

  PetscInt petsc_blocksize;
  ierr = MatGetBlockSize(_mat, &petsc_blocksize);
  LIBMESH_CHKERR(ierr);
  libmesh_assert_equal_to (blocksize, static_cast<numeric_index_type>(petsc_blocksize));
#endif

  // These casts are required for PETSc <= 2.1.5
  ierr = MatSetValuesBlocked(_mat,
                             n_brows, numeric_petsc_cast(brows.data()),
                             n_bcols, numeric_petsc_cast(bcols.data()),
                             pPS(const_cast<T*>(dm.get_values().data())),
                             ADD_VALUES);
  LIBMESH_CHKERR(ierr);
}





template <typename T>
void PetscMatrix<T>::_get_submatrix(SparseMatrix<T> & submatrix,
                                    const std::vector<numeric_index_type> & rows,
                                    const std::vector<numeric_index_type> & cols,
                                    const bool reuse_submatrix) const
{
  if (!this->closed())
    {
      libmesh_deprecated();
      libmesh_warning("The matrix must be assembled before calling PetscMatrix::create_submatrix().\n"
                      "Please update your code, as this warning will become an error in a future release.");
      const_cast<PetscMatrix<T> *>(this)->close();
    }

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T> * petsc_submatrix = cast_ptr<PetscMatrix<T> *>(&submatrix);

  // If we're not reusing submatrix and submatrix is already initialized
  // then we need to clear it, otherwise we get a memory leak.
  if (!reuse_submatrix && submatrix.initialized())
    submatrix.clear();

  // Construct row and column index sets.
  PetscErrorCode ierr=0;
  IS isrow, iscol;

  ierr = ISCreateGeneral(this->comm().get(),
                         cast_int<PetscInt>(rows.size()),
                         numeric_petsc_cast(rows.data()),
                         PETSC_USE_POINTER,
                         &isrow); LIBMESH_CHKERR(ierr);

  ierr = ISCreateGeneral(this->comm().get(),
                         cast_int<PetscInt>(cols.size()),
                         numeric_petsc_cast(cols.data()),
                         PETSC_USE_POINTER,
                         &iscol); LIBMESH_CHKERR(ierr);

  // Extract submatrix
  ierr = LibMeshCreateSubMatrix(_mat,
                                isrow,
                                iscol,
                                (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                                &(petsc_submatrix->_mat));  LIBMESH_CHKERR(ierr);

  // Specify that the new submatrix is initialized and close it.
  petsc_submatrix->_is_initialized = true;
  petsc_submatrix->close();

  // Clean up PETSc data structures
  ierr = ISDestroy(&isrow); LIBMESH_CHKERR(ierr);
  ierr = ISDestroy(&iscol); LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscMatrix<T>::get_diagonal (NumericVector<T> & dest) const
{
  // Make sure the NumericVector passed in is really a PetscVector
  PetscVector<T> & petsc_dest = cast_ref<PetscVector<T> &>(dest);

  // Needs a const_cast since PETSc does not work with const.
  PetscErrorCode ierr =
    MatGetDiagonal(const_cast<PetscMatrix<T> *>(this)->mat(),petsc_dest.vec()); LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscMatrix<T>::get_transpose (SparseMatrix<T> & dest) const
{
  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T> & petsc_dest = cast_ref<PetscMatrix<T> &>(dest);

  // If we aren't reusing the matrix then need to clear dest,
  // otherwise we get a memory leak
  if (&petsc_dest != this)
    dest.clear();

  PetscErrorCode ierr;
  if (&petsc_dest == this)
    // The MAT_REUSE_MATRIX flag was replaced by MAT_INPLACE_MATRIX
    // in PETSc 3.7.0
#if PETSC_VERSION_LESS_THAN(3,7,0)
    ierr = MatTranspose(_mat,MAT_REUSE_MATRIX,&petsc_dest._mat);
#else
    ierr = MatTranspose(_mat, MAT_INPLACE_MATRIX, &petsc_dest._mat);
#endif
  else
    ierr = MatTranspose(_mat,MAT_INITIAL_MATRIX,&petsc_dest._mat);
  LIBMESH_CHKERR(ierr);

  // Specify that the transposed matrix is initialized and close it.
  petsc_dest._is_initialized = true;
  petsc_dest.close();
}





template <typename T>
void PetscMatrix<T>::close ()
{
  semiparallel_only();

  // BSK - 1/19/2004
  // strictly this check should be OK, but it seems to
  // fail on matrix-free matrices.  Do they falsely
  // state they are assembled?  Check with the developers...
  //   if (this->closed())
  //     return;

  PetscErrorCode ierr=0;

  ierr = MatAssemblyBegin (_mat, MAT_FINAL_ASSEMBLY);
  LIBMESH_CHKERR(ierr);
  ierr = MatAssemblyEnd   (_mat, MAT_FINAL_ASSEMBLY);
  LIBMESH_CHKERR(ierr);
}

template <typename T>
void PetscMatrix<T>::flush ()
{
  semiparallel_only();

  PetscErrorCode ierr=0;

  ierr = MatAssemblyBegin (_mat, MAT_FLUSH_ASSEMBLY);
  LIBMESH_CHKERR(ierr);
  ierr = MatAssemblyEnd   (_mat, MAT_FLUSH_ASSEMBLY);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
numeric_index_type PetscMatrix<T>::m () const
{
  libmesh_assert (this->initialized());

  PetscInt petsc_m=0, petsc_n=0;
  PetscErrorCode ierr=0;

  ierr = MatGetSize (_mat, &petsc_m, &petsc_n);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(petsc_m);
}

template <typename T>
numeric_index_type PetscMatrix<T>::local_m () const
{
  libmesh_assert (this->initialized());

  PetscInt m = 0;

  auto ierr = MatGetLocalSize (_mat, &m, NULL);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(m);
}

template <typename T>
numeric_index_type PetscMatrix<T>::n () const
{
  libmesh_assert (this->initialized());

  PetscInt petsc_m=0, petsc_n=0;
  PetscErrorCode ierr=0;

  ierr = MatGetSize (_mat, &petsc_m, &petsc_n);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(petsc_n);
}

template <typename T>
numeric_index_type PetscMatrix<T>::local_n () const
{
  libmesh_assert (this->initialized());

  PetscInt n = 0;

  auto ierr = MatGetLocalSize (_mat, NULL, &n);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(n);
}

template <typename T>
void PetscMatrix<T>::get_local_size (numeric_index_type & m,
                                     numeric_index_type & n) const
{
  libmesh_assert (this->initialized());

  PetscInt petsc_m = 0, petsc_n = 0;

  auto ierr = MatGetLocalSize (_mat, &petsc_m, &petsc_n);
  LIBMESH_CHKERR(ierr);

  m = static_cast<numeric_index_type>(petsc_m);
  n = static_cast<numeric_index_type>(petsc_n);
}

template <typename T>
numeric_index_type PetscMatrix<T>::row_start () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;
  PetscErrorCode ierr=0;

  ierr = MatGetOwnershipRange(_mat, &start, &stop);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(start);
}



template <typename T>
numeric_index_type PetscMatrix<T>::row_stop () const
{
  libmesh_assert (this->initialized());

  PetscInt start=0, stop=0;
  PetscErrorCode ierr=0;

  ierr = MatGetOwnershipRange(_mat, &start, &stop);
  LIBMESH_CHKERR(ierr);

  return static_cast<numeric_index_type>(stop);
}



template <typename T>
void PetscMatrix<T>::set (const numeric_index_type i,
                          const numeric_index_type j,
                          const T value)
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
                      &petsc_value, INSERT_VALUES);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscMatrix<T>::add (const numeric_index_type i,
                          const numeric_index_type j,
                          const T value)
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscInt i_val=i, j_val=j;

  PetscScalar petsc_value = static_cast<PetscScalar>(value);
  ierr = MatSetValues(_mat, 1, &i_val, 1, &j_val,
                      &petsc_value, ADD_VALUES);
  LIBMESH_CHKERR(ierr);
}



template <typename T>
void PetscMatrix<T>::add_matrix(const DenseMatrix<T> & dm,
                                const std::vector<numeric_index_type> & dof_indices)
{
  this->add_matrix (dm, dof_indices, dof_indices);
}







template <typename T>
void PetscMatrix<T>::add (const T a_in, const SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());

  // sanity check. but this cannot avoid
  // crash due to incompatible sparsity structure...
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  PetscScalar a = static_cast<PetscScalar>      (a_in);
  const PetscMatrix<T> * X = cast_ptr<const PetscMatrix<T> *> (&X_in);

  libmesh_assert (X);

  PetscErrorCode ierr=0;

  // the matrix from which we copy the values has to be assembled/closed
  libmesh_assert(X->closed());

  semiparallel_only();

  ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
  LIBMESH_CHKERR(ierr);
}


template <typename T>
void PetscMatrix<T>::mat_mult (SparseMatrix<T> & X_in, SparseMatrix<T> & Y_out)
{
  libmesh_assert (this->initialized());

   // sanity check.
   libmesh_assert_equal_to (this->n(), X_in.m());
   libmesh_assert_equal_to (this->m(), Y_out.m());
   libmesh_assert_equal_to (X_in.n(), Y_out.n());

   const PetscMatrix<T> * X = cast_ptr<const PetscMatrix<T> *> (&X_in);
   PetscMatrix<T> * Y = cast_ptr<PetscMatrix<T> *> (&Y_out);

   libmesh_assert (X);
   libmesh_assert (Y);

   PetscErrorCode ierr=0;

   // the matrix from which we copy the values has to be assembled/closed
   libmesh_assert(X->closed());

   semiparallel_only();

   ierr = MatMatMult(_mat, X->_mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y->_mat);
   LIBMESH_CHKERR(ierr);

}

template <typename T>
void PetscMatrix<T>::add_sparse_matrix (const SparseMatrix<T> & spm,
                        const std::vector<numeric_index_type> & rows,
                        const std::vector<numeric_index_type> & cols,
                        const T scalar)
{
  libmesh_assert_equal_to(spm.m(), rows.size());
  libmesh_assert_equal_to(spm.n(), cols.size());

  for(unsigned int i = 0; i<rows.size(); ++i)
  {
    std::vector<numeric_index_type> indices;
    std::vector<T> values;
    spm.get_row(i, indices, values);
    for (unsigned int j=0; j<indices.size(); ++j)
      this->add (rows[i], cols[indices[j]], scalar*values[j]);
  }
}



template <typename T>
T PetscMatrix<T>::operator () (const numeric_index_type i_in,
                               const numeric_index_type j_in) const
{
  libmesh_assert (this->initialized());

  // PETSc 2.2.1 & newer
  const PetscScalar * petsc_row;
  const PetscInt    * petsc_cols;

  // If the entry is not in the sparse matrix, it is 0.
  T value=0.;

  PetscErrorCode
    ierr=0;
  PetscInt
    ncols=0,
    i_val=static_cast<PetscInt>(i_in),
    j_val=static_cast<PetscInt>(j_in);


  // the matrix needs to be closed for this to work
  // this->close();
  // but closing it is a semiparallel operation; we want operator()
  // to run on one processor.
  libmesh_assert(this->closed());

  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
  LIBMESH_CHKERR(ierr);

  // Perform a binary search to find the contiguous index in
  // petsc_cols (resp. petsc_row) corresponding to global index j_val
  std::pair<const PetscInt *, const PetscInt *> p =
    std::equal_range (petsc_cols, petsc_cols + ncols, j_val);

  // Found an entry for j_val
  if (p.first != p.second)
    {
      // The entry in the contiguous row corresponding
      // to the j_val column of interest
      const std::size_t j =
        std::distance (const_cast<PetscInt *>(petsc_cols),
                       const_cast<PetscInt *>(p.first));

      libmesh_assert_less (j, ncols);
      libmesh_assert_equal_to (petsc_cols[j], j_val);

      value = static_cast<T> (petsc_row[j]);
    }

  ierr  = MatRestoreRow(_mat, i_val,
                        &ncols, &petsc_cols, &petsc_row);
  LIBMESH_CHKERR(ierr);

  return value;
}

template <typename T>
void PetscMatrix<T>::get_row (numeric_index_type i_in,
                              std::vector<numeric_index_type> & indices,
                              std::vector<T> & values) const
{
  libmesh_assert (this->initialized());

  const PetscScalar * petsc_row;
  const PetscInt    * petsc_cols;

  PetscErrorCode ierr=0;
  PetscInt
    ncols=0,
    i_val = static_cast<PetscInt>(i_in);

  // the matrix needs to be closed for this to work
  // this->close();
  // but closing it is a semiparallel operation; we want operator()
  // to run on one processor.
  libmesh_assert(this->closed());

  // PETSc makes no effort at being thread safe. Helgrind complains about
  // possible data races even just in PetscFunctionBegin (due to things
  // like stack counter incrementing). Perhaps we could ignore
  // this, but there are legitimate data races for Mat data members like
  // mat->getrowactive between MatGetRow and MatRestoreRow. Moreover,
  // there could be a write into mat->rowvalues during MatGetRow from
  // one thread while we are attempting to read from mat->rowvalues
  // (through petsc_cols) during data copy in another thread. So
  // the safe thing to do is to lock the whole method

#ifdef LIBMESH_HAVE_CXX11_THREAD
  std::lock_guard<std::mutex>
#else
  Threads::spin_mutex::scoped_lock
#endif
    lock(_petsc_matrix_mutex);

  ierr = MatGetRow(_mat, i_val, &ncols, &petsc_cols, &petsc_row);
  LIBMESH_CHKERR(ierr);

  // Copy the data
  indices.resize(static_cast<std::size_t>(ncols));
  values.resize(static_cast<std::size_t>(ncols));

  for (auto i : index_range(indices))
  {
    indices[i] = static_cast<numeric_index_type>(petsc_cols[i]);
    values[i] = static_cast<T>(petsc_row[i]);
  }

  ierr  = MatRestoreRow(_mat, i_val,
                        &ncols, &petsc_cols, &petsc_row);
  LIBMESH_CHKERR(ierr);
}


template <typename T>
bool PetscMatrix<T>::closed() const
{
  libmesh_assert (this->initialized());

  PetscErrorCode ierr=0;
  PetscBool assembled;

  ierr = MatAssembled(_mat, &assembled);
  LIBMESH_CHKERR(ierr);

  return (assembled == PETSC_TRUE);
}

template <typename T>
void PetscMatrix<T>::set_destroy_mat_on_exit(bool destroy)
{
  this->_destroy_mat_on_exit = destroy;
}


template <typename T>
void PetscMatrix<T>::swap(PetscMatrix<T> & m_in)
{
  std::swap(_mat, m_in._mat);
  std::swap(_destroy_mat_on_exit, m_in._destroy_mat_on_exit);
}



//------------------------------------------------------------------
// Explicit instantiations
template class PetscMatrix<Number>;

} // namespace libMesh


#endif // #ifdef LIBMESH_HAVE_PETSC
