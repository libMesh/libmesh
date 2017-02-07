// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#if !PETSC_RELEASE_LESS_THAN(3,6,0)
# include "libmesh/ignore_warnings.h"
# include "petsc/private/matimpl.h"
# include "libmesh/restore_warnings.h"
#elif !PETSC_VERSION_LESS_THAN(3,5,0)
# include "libmesh/ignore_warnings.h"
# include "petsc-private/matimpl.h"
# include "libmesh/restore_warnings.h"
#endif

// For some reason, the blocked matrix API calls below seem to break with PETSC 3.2 & presumably earier.
// For example:
// [0]PETSC ERROR: --------------------- Error Message ------------------------------------
// [0]PETSC ERROR: Nonconforming object sizes!
// [0]PETSC ERROR: Attempt to set block size 3 with BAIJ 1!
// [0]PETSC ERROR: ------------------------------------------------------------------------
// so as a cowardly workaround disable the functionality in this translation unit for older PETSc's
#if PETSC_VERSION_LESS_THAN(3,3,0)
#  undef LIBMESH_ENABLE_BLOCKED_STORAGE
#endif



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

  for (std::size_t nn=0; nn<n_nz.size(); nn += blocksize)
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
  _destroy_mat_on_exit(true)
{}



// Constructor taking an existing Mat but not the responsibility
// for destroying it
template <typename T>
PetscMatrix<T>::PetscMatrix(Mat mat_in,
                            const Parallel::Communicator & comm_in) :
  SparseMatrix<T>(comm_in),
  _destroy_mat_on_exit(false)
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

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  PetscInt blocksize  = static_cast<PetscInt>(blocksize_in);
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
      ierr = MatSetBlockSize(_mat, blocksize);
      LIBMESH_CHKERR(ierr);
      ierr = MatSeqBAIJSetPreallocation(_mat, blocksize, n_nz/blocksize, PETSC_NULL);
      LIBMESH_CHKERR(ierr);
      ierr = MatMPIBAIJSetPreallocation(_mat, blocksize,
                                        n_nz/blocksize, PETSC_NULL,
                                        n_oz/blocksize, PETSC_NULL);
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {
      ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
      LIBMESH_CHKERR(ierr);
      ierr = MatSeqAIJSetPreallocation(_mat, n_nz, PETSC_NULL);
      LIBMESH_CHKERR(ierr);
      ierr = MatMPIAIJSetPreallocation(_mat, n_nz, PETSC_NULL, n_oz, PETSC_NULL);
      LIBMESH_CHKERR(ierr);
    }

  // Make it an error for PETSc to allocate new nonzero entries during assembly
#if PETSC_VERSION_LESS_THAN(3,0,0)
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR);
#else
  ierr = MatSetOption(_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
#endif
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

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  PetscInt blocksize  = static_cast<PetscInt>(blocksize_in);
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
      ierr = MatSetBlockSize(_mat, blocksize);
      LIBMESH_CHKERR(ierr);

      // transform the per-entry n_nz and n_oz arrays into their block counterparts.
      std::vector<numeric_index_type> b_n_nz, b_n_oz;

      transform_preallocation_arrays (blocksize,
                                      n_nz, n_oz,
                                      b_n_nz, b_n_oz);

      ierr = MatSeqBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? libmesh_nullptr : &b_n_nz[0]));
      LIBMESH_CHKERR(ierr);

      ierr = MatMPIBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? libmesh_nullptr : &b_n_nz[0]),
                                         0,
                                         numeric_petsc_cast(b_n_oz.empty() ? libmesh_nullptr : &b_n_oz[0]));
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {

      ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
      LIBMESH_CHKERR(ierr);
      ierr = MatSeqAIJSetPreallocation (_mat,
                                        0,
                                        numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]));
      LIBMESH_CHKERR(ierr);
      ierr = MatMPIAIJSetPreallocation (_mat,
                                        0,
                                        numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]),
                                        0,
                                        numeric_petsc_cast(n_oz.empty() ? libmesh_nullptr : &n_oz[0]));
      LIBMESH_CHKERR(ierr);
    }

  // Is prefix information available somewhere? Perhaps pass in the system name?
  ierr = MatSetOptionsPrefix(_mat, "");
  LIBMESH_CHKERR(ierr);
  ierr = MatSetFromOptions(_mat);
  LIBMESH_CHKERR(ierr);


  this->zero();
}


template <typename T>
void PetscMatrix<T>::init ()
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

#ifdef LIBMESH_ENABLE_BLOCKED_STORAGE
  PetscInt blocksize  = static_cast<PetscInt>(this->_dof_map->block_size());
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
      ierr = MatSetBlockSize(_mat, blocksize);
      LIBMESH_CHKERR(ierr);

      // transform the per-entry n_nz and n_oz arrays into their block counterparts.
      std::vector<numeric_index_type> b_n_nz, b_n_oz;

      transform_preallocation_arrays (blocksize,
                                      n_nz, n_oz,
                                      b_n_nz, b_n_oz);

      ierr = MatSeqBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? libmesh_nullptr : &b_n_nz[0]));
      LIBMESH_CHKERR(ierr);

      ierr = MatMPIBAIJSetPreallocation (_mat,
                                         blocksize,
                                         0,
                                         numeric_petsc_cast(b_n_nz.empty() ? libmesh_nullptr : &b_n_nz[0]),
                                         0,
                                         numeric_petsc_cast(b_n_oz.empty() ? libmesh_nullptr : &b_n_oz[0]));
      LIBMESH_CHKERR(ierr);
    }
  else
#endif
    {
      // no block storage case
      ierr = MatSetType(_mat, MATAIJ); // Automatically chooses seqaij or mpiaij
      LIBMESH_CHKERR(ierr);

      ierr = MatSeqAIJSetPreallocation (_mat,
                                        0,
                                        numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]));
      LIBMESH_CHKERR(ierr);
      ierr = MatMPIAIJSetPreallocation (_mat,
                                        0,
                                        numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]),
                                        0,
                                        numeric_petsc_cast(n_oz.empty() ? libmesh_nullptr : &n_oz[0]));
      LIBMESH_CHKERR(ierr);
    }

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
  libmesh_assert(this->_dof_map);
  libmesh_assert(this->initialized());

  const std::vector<numeric_index_type> & n_nz = this->_dof_map->get_n_nz();
  const std::vector<numeric_index_type> & n_oz = this->_dof_map->get_n_oz();

  PetscErrorCode ierr = 0;

  {
#if !PETSC_VERSION_LESS_THAN(3,5,0)
    ierr = (*_mat->ops->destroy)(_mat);
    LIBMESH_CHKERR(ierr);
    _mat->ops->destroy = libmesh_nullptr;
    _mat->preallocated = PETSC_FALSE;
    _mat->assembled    = PETSC_FALSE;
    _mat->was_assembled = PETSC_FALSE;
    ++_mat->nonzerostate;
    ierr = PetscObjectStateIncrease((PetscObject)_mat);
    LIBMESH_CHKERR(ierr);
#else
    libmesh_error_msg("PetscMatrix::update_preallocation_and_zero() requires PETSc 3.5.0 or greater to work correctly.");
#endif

    ierr = MatSetType(_mat,MATAIJ);
    LIBMESH_CHKERR(ierr);

    ierr = MatSeqAIJSetPreallocation (_mat,
                                      0,
                                      numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]));
    LIBMESH_CHKERR(ierr);
    ierr = MatMPIAIJSetPreallocation (_mat,
                                      0,
                                      numeric_petsc_cast(n_nz.empty() ? libmesh_nullptr : &n_nz[0]),
                                      0,
                                      numeric_petsc_cast(n_oz.empty() ? libmesh_nullptr : &n_oz[0]));
    LIBMESH_CHKERR(ierr);
  }

  this->zero();
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

#if PETSC_RELEASE_LESS_THAN(3,1,1)
  if(!rows.empty())
    ierr = MatZeroRows(_mat, rows.size(),
                       numeric_petsc_cast(&rows[0]), diag_value);
  else
    ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value);
#else
  // As of petsc-dev at the time of 3.1.0, MatZeroRows now takes two additional
  // optional arguments.  The optional arguments (x,b) can be used to specify the
  // solutions for the zeroed rows (x) and right hand side (b) to update.
  // Could be useful for setting boundary conditions...
  if(!rows.empty())
    ierr = MatZeroRows(_mat, cast_int<PetscInt>(rows.size()),
                       numeric_petsc_cast(&rows[0]), diag_value,
                       PETSC_NULL, PETSC_NULL);
  else
    ierr = MatZeroRows(_mat, 0, PETSC_NULL, diag_value, PETSC_NULL,
                       PETSC_NULL);
#endif

  LIBMESH_CHKERR(ierr);
}

template <typename T>
void PetscMatrix<T>::clear ()
{
  PetscErrorCode ierr=0;

  if ((this->initialized()) && (this->_destroy_mat_on_exit))
    {
      semiparallel_only();

      ierr = LibMeshMatDestroy (&_mat);
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

  // libmesh_assert (this->closed());
  this->close();

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
  ierr = LibMeshPetscViewerDestroy (&petsc_viewer);
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
  this->close();

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

  // These casts are required for PETSc <= 2.1.5
  ierr = MatSetValues(_mat,
                      n_rows, numeric_petsc_cast(&rows[0]),
                      n_cols, numeric_petsc_cast(&cols[0]),
                      const_cast<PetscScalar *>(&dm.get_values()[0]),
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
                             n_brows, numeric_petsc_cast(&brows[0]),
                             n_bcols, numeric_petsc_cast(&bcols[0]),
                             const_cast<PetscScalar *>(&dm.get_values()[0]),
                             ADD_VALUES);
  LIBMESH_CHKERR(ierr);
}





template <typename T>
void PetscMatrix<T>::_get_submatrix(SparseMatrix<T> & submatrix,
                                    const std::vector<numeric_index_type> & rows,
                                    const std::vector<numeric_index_type> & cols,
                                    const bool reuse_submatrix) const
{
  // Can only extract submatrices from closed matrices
  this->close();

  // Make sure the SparseMatrix passed in is really a PetscMatrix
  PetscMatrix<T> * petsc_submatrix = cast_ptr<PetscMatrix<T> *>(&submatrix);

  // If we're not reusing submatrix and submatrix is already initialized
  // then we need to clear it, otherwise we get a memory leak.
  if( !reuse_submatrix && submatrix.initialized() )
    submatrix.clear();

  // Construct row and column index sets.
  PetscErrorCode ierr=0;
  IS isrow, iscol;

  ierr = ISCreateLibMesh(this->comm().get(),
                         rows.size(),
                         numeric_petsc_cast(&rows[0]),
                         PETSC_USE_POINTER,
                         &isrow); LIBMESH_CHKERR(ierr);

  ierr = ISCreateLibMesh(this->comm().get(),
                         cols.size(),
                         numeric_petsc_cast(&cols[0]),
                         PETSC_USE_POINTER,
                         &iscol); LIBMESH_CHKERR(ierr);

  // Extract submatrix
  ierr = MatGetSubMatrix(_mat,
                         isrow,
                         iscol,
#if PETSC_RELEASE_LESS_THAN(3,0,1)
                         PETSC_DECIDE,
#endif
                         (reuse_submatrix ? MAT_REUSE_MATRIX : MAT_INITIAL_MATRIX),
                         &(petsc_submatrix->_mat));  LIBMESH_CHKERR(ierr);

  // Specify that the new submatrix is initialized and close it.
  petsc_submatrix->_is_initialized = true;
  petsc_submatrix->close();

  // Clean up PETSc data structures
  ierr = LibMeshISDestroy(&isrow); LIBMESH_CHKERR(ierr);
  ierr = LibMeshISDestroy(&iscol); LIBMESH_CHKERR(ierr);
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
  if(&petsc_dest != this)
    dest.clear();

  PetscErrorCode ierr;
#if PETSC_VERSION_LESS_THAN(3,0,0)
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,PETSC_NULL);
  else
    ierr = MatTranspose(_mat,&petsc_dest._mat);
  LIBMESH_CHKERR(ierr);
#else
  // FIXME - we can probably use MAT_REUSE_MATRIX in more situations
  if (&petsc_dest == this)
    ierr = MatTranspose(_mat,MAT_REUSE_MATRIX,&petsc_dest._mat);
  else
    ierr = MatTranspose(_mat,MAT_INITIAL_MATRIX,&petsc_dest._mat);
  LIBMESH_CHKERR(ierr);
#endif

  // Specify that the transposed matrix is initialized and close it.
  petsc_dest._is_initialized = true;
  petsc_dest.close();
}





template <typename T>
void PetscMatrix<T>::close () const
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
void PetscMatrix<T>::add (const T a_in, SparseMatrix<T> & X_in)
{
  libmesh_assert (this->initialized());

  // sanity check. but this cannot avoid
  // crash due to incompatible sparsity structure...
  libmesh_assert_equal_to (this->m(), X_in.m());
  libmesh_assert_equal_to (this->n(), X_in.n());

  PetscScalar a = static_cast<PetscScalar>      (a_in);
  PetscMatrix<T> * X = cast_ptr<PetscMatrix<T> *> (&X_in);

  libmesh_assert (X);

  PetscErrorCode ierr=0;

  // the matrix from which we copy the values has to be assembled/closed
  libmesh_assert(X->closed());

  semiparallel_only();

  ierr = MatAXPY(_mat, a, X->_mat, DIFFERENT_NONZERO_PATTERN);
  LIBMESH_CHKERR(ierr);
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
    std::equal_range (&petsc_cols[0], &petsc_cols[0] + ncols, j_val);

  // Found an entry for j_val
  if (p.first != p.second)
    {
      // The entry in the contiguous row corresponding
      // to the j_val column of interest
      const std::size_t j =
        std::distance (const_cast<PetscInt *>(&petsc_cols[0]),
                       const_cast<PetscInt *>(p.first));

      libmesh_assert_less (static_cast<PetscInt>(j), ncols);
      libmesh_assert_equal_to (petsc_cols[j], j_val);

      value = static_cast<T> (petsc_row[j]);
    }

  ierr  = MatRestoreRow(_mat, i_val,
                        &ncols, &petsc_cols, &petsc_row);
  LIBMESH_CHKERR(ierr);

  return value;
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
