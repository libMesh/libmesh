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


#include "libmesh/libmesh_common.h"
#include "libmesh/parallel.h"

// C++ Includes
#include <cstdio> // for std::sprintf
#include <set>
#include <numeric> // for std::partial_sum

// Local Include
#include "libmesh/libmesh_version.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
//#include "libmesh/mesh_tools.h"
#include "libmesh/elem.h"
#include "libmesh/xdr_cxx.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"



// Anonymous namespace for implementation details.
namespace {

using libMesh::DofObject;
using libMesh::Number;
using libMesh::cast_int;

// Comments:
// ---------
// - The max_io_blksize governs how many nodes or elements will be
// treated as a single block when performing parallel IO on large
// systems.
// - This parameter only loosely affects the size of the actual IO
// buffer as this depends on the number of components a given
// variable has for the nodes/elements in the block.
// - When reading/writing each processor uses an ID map which is
// 3*io_blksize*sizeof(dof_id_type) bytes long, so with unsigned int
// and // io_blksize=256000 we would expect that buffer alone to be
// ~3Mb.
// - In general, an increase in max_io_blksize should increase the
// efficiency of large parallel read/writes by reducing the number
// of MPI messages at the expense of memory.
// - If the library exhausts memory during IO you might reduce this
// parameter.

const std::size_t max_io_blksize = 256000;


/**
 * Comparison object to use with DofObject pointers.  This sorts by id(),
 * so when we iterate over a set of DofObjects we visit the objects in
 * order of increasing ID.
 */
struct CompareDofObjectsByID
{
  bool operator()(const DofObject * a,
                  const DofObject * b) const
  {
    libmesh_assert (a);
    libmesh_assert (b);

    return a->id() < b->id();
  }
};

/**
 *
 */
template <typename InValType>
class ThreadedIO
{
private:
  libMesh::Xdr & _io;
  std::vector<InValType> & _data;

public:
  ThreadedIO (libMesh::Xdr & io, std::vector<InValType> & data) :
    _io(io),
    _data(data)
  {}

  void operator()()
  {
    if (_data.empty()) return;
    _io.data_stream (&_data[0], cast_int<unsigned int>(_data.size()));
  }
};
}


namespace libMesh
{


// ------------------------------------------------------------
// System class implementation
void System::read_header (Xdr & io,
                          const std::string & version,
                          const bool read_header_in,
                          const bool read_additional_data,
                          const bool read_legacy_format)
{
  // This method implements the input of a
  // System object, embedded in the output of
  // an EquationSystems<T_sys>.  This warrants some
  // documentation.  The output file essentially
  // consists of 5 sections:
  //
  // for this system
  //
  //   5.) The number of variables in the system (unsigned int)
  //
  //   for each variable in the system
  //
  //     6.) The name of the variable (string)
  //
  //     6.1.) Variable subdmains
  //
  //     7.) Combined in an FEType:
  //         - The approximation order(s) of the variable
  //           (Order Enum, cast to int/s)
  //         - The finite element family/ies of the variable
  //           (FEFamily Enum, cast to int/s)
  //
  //   end variable loop
  //
  //   8.) The number of additional vectors (unsigned int),
  //
  //     for each additional vector in the system object
  //
  //     9.) the name of the additional vector  (string)
  //
  // end system
  libmesh_assert (io.reading());

  // Possibly clear data structures and start from scratch.
  if (read_header_in)
    this->clear ();

  // Figure out if we need to read infinite element information.
  // This will be true if the version string contains " with infinite elements"
  const bool read_ifem_info =
    (version.rfind(" with infinite elements") < version.size()) ||
    libMesh::on_command_line ("--read_ifem_systems");


  {
    // 5.)
    // Read the number of variables in the system
    unsigned int nv=0;
    if (this->processor_id() == 0)
      io.data (nv);
    this->comm().broadcast(nv);

    _written_var_indices.clear();
    _written_var_indices.resize(nv, 0);

    for (unsigned int var=0; var<nv; var++)
      {
        // 6.)
        // Read the name of the var-th variable
        std::string var_name;
        if (this->processor_id() == 0)
          io.data (var_name);
        this->comm().broadcast(var_name);

        // 6.1.)
        std::set<subdomain_id_type> domains;
        if (io.version() >= LIBMESH_VERSION_ID(0,7,2))
          {
            std::vector<subdomain_id_type> domain_array;
            if (this->processor_id() == 0)
              io.data (domain_array);
            for (std::vector<subdomain_id_type>::iterator it = domain_array.begin(); it != domain_array.end(); ++it)
              domains.insert(*it);
          }
        this->comm().broadcast(domains);

        // 7.)
        // Read the approximation order(s) of the var-th variable
        int order=0;
        if (this->processor_id() == 0)
          io.data (order);
        this->comm().broadcast(order);


        // do the same for infinite element radial_order
        int rad_order=0;
        if (read_ifem_info)
          {
            if (this->processor_id() == 0)
              io.data(rad_order);
            this->comm().broadcast(rad_order);
          }

        // Read the finite element type of the var-th variable
        int fam=0;
        if (this->processor_id() == 0)
          io.data (fam);
        this->comm().broadcast(fam);
        FEType type;
        type.order  = static_cast<Order>(order);
        type.family = static_cast<FEFamily>(fam);

        // Check for incompatibilities.  The shape function indexing was
        // changed for the monomial and xyz finite element families to
        // simplify extension to arbitrary p.  The consequence is that
        // old restart files will not be read correctly.  This is expected
        // to be an unlikely occurance, but catch it anyway.
        if (read_legacy_format)
          if ((type.family == MONOMIAL || type.family == XYZ) &&
              ((type.order.get_order() > 2 && this->get_mesh().mesh_dimension() == 2) ||
               (type.order.get_order() > 1 && this->get_mesh().mesh_dimension() == 3)))
            {
              libmesh_here();
              libMesh::out << "*****************************************************************\n"
                           << "* WARNING: reading a potentially incompatible restart file!!!   *\n"
                           << "*  contact libmesh-users@lists.sourceforge.net for more details *\n"
                           << "*****************************************************************"
                           << std::endl;
            }

        // Read additional information for infinite elements
        int radial_fam=0;
        int i_map=0;
        if (read_ifem_info)
          {
            if (this->processor_id() == 0)
              io.data (radial_fam);
            this->comm().broadcast(radial_fam);
            if (this->processor_id() == 0)
              io.data (i_map);
            this->comm().broadcast(i_map);
          }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

        type.radial_order  = static_cast<Order>(rad_order);
        type.radial_family = static_cast<FEFamily>(radial_fam);
        type.inf_map       = static_cast<InfMapType>(i_map);

#endif

        if (read_header_in)
          {
            if (domains.empty())
              _written_var_indices[var] = this->add_variable (var_name, type);
            else
              _written_var_indices[var] = this->add_variable (var_name, type, &domains);
          }
        else
          _written_var_indices[var] = this->variable_number(var_name);
      }
  }

  // 8.)
  // Read the number of additional vectors.
  unsigned int nvecs=0;
  if (this->processor_id() == 0)
    io.data (nvecs);
  this->comm().broadcast(nvecs);

  // If nvecs > 0, this means that write_additional_data
  // was true when this file was written.  We will need to
  // make use of this fact later.
  this->_additional_data_written = nvecs;

  for (unsigned int vec=0; vec<nvecs; vec++)
    {
      // 9.)
      // Read the name of the vec-th additional vector
      std::string vec_name;
      if (this->processor_id() == 0)
        io.data (vec_name);
      this->comm().broadcast(vec_name);

      if (read_additional_data)
        {
          // Systems now can handle adding post-initialization vectors
          //  libmesh_assert(this->_can_add_vectors);
          // Some systems may have added their own vectors already
          //  libmesh_assert_equal_to (this->_vectors.count(vec_name), 0);

          this->add_vector(vec_name);
        }
    }
}



void System::read_legacy_data (Xdr & io,
                               const bool read_additional_data)
{
  libmesh_deprecated();

  // This method implements the output of the vectors
  // contained in this System object, embedded in the
  // output of an EquationSystems<T_sys>.
  //
  //   10.) The global solution vector, re-ordered to be node-major
  //       (More on this later.)
  //
  //      for each additional vector in the object
  //
  //      11.) The global additional vector, re-ordered to be
  //           node-major (More on this later.)
  libmesh_assert (io.reading());

  // read and reordering buffers
  std::vector<Number> global_vector;
  std::vector<Number> reordered_vector;

  // 10.)
  // Read and set the solution vector
  {
    if (this->processor_id() == 0)
      io.data (global_vector);
    this->comm().broadcast(global_vector);

    // Remember that the stored vector is node-major.
    // We need to put it into whatever application-specific
    // ordering we may have using the dof_map.
    reordered_vector.resize(global_vector.size());

    //libMesh::out << "global_vector.size()=" << global_vector.size() << std::endl;
    //libMesh::out << "this->n_dofs()=" << this->n_dofs() << std::endl;

    libmesh_assert_equal_to (global_vector.size(), this->n_dofs());

    dof_id_type cnt=0;

    const unsigned int sys = this->number();
    const unsigned int nv  = cast_int<unsigned int>
      (this->_written_var_indices.size());
    libmesh_assert_less_equal (nv, this->n_vars());

    for (unsigned int data_var=0; data_var<nv; data_var++)
      {
        const unsigned int var = _written_var_indices[data_var];

        // First reorder the nodal DOF values
        {
          MeshBase::node_iterator
            it  = this->get_mesh().nodes_begin(),
            end = this->get_mesh().nodes_end();

          for (; it != end; ++it)
            for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
              {
                libmesh_assert_not_equal_to ((*it)->dof_number(sys, var, index),
                                             DofObject::invalid_id);

                libmesh_assert_less (cnt, global_vector.size());

                reordered_vector[(*it)->dof_number(sys, var, index)] =
                  global_vector[cnt++];
              }
        }

        // Then reorder the element DOF values
        {
          MeshBase::element_iterator
            it  = this->get_mesh().active_elements_begin(),
            end = this->get_mesh().active_elements_end();

          for (; it != end; ++it)
            for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
              {
                libmesh_assert_not_equal_to ((*it)->dof_number(sys, var, index),
                                             DofObject::invalid_id);

                libmesh_assert_less (cnt, global_vector.size());

                reordered_vector[(*it)->dof_number(sys, var, index)] =
                  global_vector[cnt++];
              }
        }
      }

    *(this->solution) = reordered_vector;
  }

  // For each additional vector, simply go through the list.
  // ONLY attempt to do this IF additional data was actually
  // written to the file for this system (controlled by the
  // _additional_data_written flag).
  if (this->_additional_data_written)
    {
      const std::size_t nvecs = this->_vectors.size();

      // If the number of additional vectors written is non-zero, and
      // the number of additional vectors we have is non-zero, and
      // they don't match, then something is wrong and we can't be
      // sure we're reading data into the correct places.
      if (read_additional_data && nvecs &&
          nvecs != this->_additional_data_written)
        libmesh_error_msg
          ("Additional vectors in file do not match system");

      std::map<std::string, NumericVector<Number> * >::iterator
        pos = this->_vectors.begin();

      for (std::size_t i = 0; i != this->_additional_data_written; ++i)
        {
          // 11.)
          // Read the values of the vec-th additional vector.
          // Prior do _not_ clear, but fill with zero, since the
          // additional vectors _have_ to have the same size
          // as the solution vector
          std::fill (global_vector.begin(), global_vector.end(), libMesh::zero);

          if (this->processor_id() == 0)
            io.data (global_vector);

          // If read_additional_data==true and we have additional vectors,
          // then we will keep this vector data; otherwise we are going to
          // throw it away.
          if (read_additional_data && nvecs)
            {
              this->comm().broadcast(global_vector);

              // Remember that the stored vector is node-major.
              // We need to put it into whatever application-specific
              // ordering we may have using the dof_map.
              std::fill (reordered_vector.begin(),
                         reordered_vector.end(),
                         libMesh::zero);

              reordered_vector.resize(global_vector.size());

              libmesh_assert_equal_to (global_vector.size(), this->n_dofs());

              dof_id_type cnt=0;

              const unsigned int sys = this->number();
              const unsigned int nv  = cast_int<unsigned int>
                (this->_written_var_indices.size());
              libmesh_assert_less_equal (nv, this->n_vars());

              for (unsigned int data_var=0; data_var<nv; data_var++)
                {
                  const unsigned int var = _written_var_indices[data_var];
                  // First reorder the nodal DOF values
                  {
                    MeshBase::node_iterator
                      it  = this->get_mesh().nodes_begin(),
                      end = this->get_mesh().nodes_end();

                    for (; it!=end; ++it)
                      for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
                        {
                          libmesh_assert_not_equal_to ((*it)->dof_number(sys, var, index),
                                                       DofObject::invalid_id);

                          libmesh_assert_less (cnt, global_vector.size());

                          reordered_vector[(*it)->dof_number(sys, var, index)] =
                            global_vector[cnt++];
                        }
                  }

                  // Then reorder the element DOF values
                  {
                    MeshBase::element_iterator
                      it  = this->get_mesh().active_elements_begin(),
                      end = this->get_mesh().active_elements_end();

                    for (; it!=end; ++it)
                      for (unsigned int index=0; index<(*it)->n_comp(sys,var); index++)
                        {
                          libmesh_assert_not_equal_to ((*it)->dof_number(sys, var, index),
                                                       DofObject::invalid_id);

                          libmesh_assert_less (cnt, global_vector.size());

                          reordered_vector[(*it)->dof_number(sys, var, index)] =
                            global_vector[cnt++];
                        }
                  }
                }

              // use the overloaded operator=(std::vector) to assign the values
              *(pos->second) = reordered_vector;
            }

          // If we've got vectors then we need to be iterating through
          // those too
          if (pos != this->_vectors.end())
            ++pos;
        }
    } // end if (_additional_data_written)
}



template <typename InValType>
void System::read_parallel_data (Xdr & io,
                                 const bool read_additional_data)
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the
   * output of an EquationSystems<T_sys>.
   *
   *   9.) The global solution vector, re-ordered to be node-major
   *       (More on this later.)
   *
   *      for each additional vector in the object
   *
   *      10.) The global additional vector, re-ordered to be
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class
   * (to be renamed later?) which provides a uniform interface to
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */
  // PerfLog pl("IO Performance",false);
  // pl.push("read_parallel_data");
  dof_id_type total_read_size = 0;

  libmesh_assert (io.reading());
  libmesh_assert (io.is_open());

  // build the ordered nodes and element maps.
  // when writing/reading parallel files we need to iterate
  // over our nodes/elements in order of increasing global id().
  // however, this is not guaranteed to be ordering we obtain
  // by using the node_iterators/element_iterators directly.
  // so build a set, sorted by id(), that provides the ordering.
  // further, for memory economy build the set but then transfer
  // its contents to vectors, which will be sorted.
  std::vector<const DofObject *> ordered_nodes, ordered_elements;
  {
    std::set<const DofObject *, CompareDofObjectsByID>
      ordered_nodes_set (this->get_mesh().local_nodes_begin(),
                         this->get_mesh().local_nodes_end());

    ordered_nodes.insert(ordered_nodes.end(),
                         ordered_nodes_set.begin(),
                         ordered_nodes_set.end());
  }
  {
    std::set<const DofObject *, CompareDofObjectsByID>
      ordered_elements_set (this->get_mesh().local_elements_begin(),
                            this->get_mesh().local_elements_end());

    ordered_elements.insert(ordered_elements.end(),
                            ordered_elements_set.begin(),
                            ordered_elements_set.end());
  }

  //  std::vector<Number> io_buffer;
  std::vector<InValType> io_buffer;

  // 9.)
  //
  // Actually read the solution components
  // for the ith system to disk
  io.data(io_buffer);

  total_read_size += cast_int<dof_id_type>(io_buffer.size());

  const unsigned int sys_num = this->number();
  const unsigned int nv      = cast_int<unsigned int>
    (this->_written_var_indices.size());
  libmesh_assert_less_equal (nv, this->n_vars());

  dof_id_type cnt=0;

  // Loop over each non-SCALAR variable and each node, and read out the value.
  for (unsigned int data_var=0; data_var<nv; data_var++)
    {
      const unsigned int var = _written_var_indices[data_var];
      if(this->variable(var).type().family != SCALAR)
        {
          // First read the node DOF values
          for (std::vector<const DofObject *>::const_iterator
                 it = ordered_nodes.begin(); it != ordered_nodes.end(); ++it)
            for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
              {
                libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                             DofObject::invalid_id);
                libmesh_assert_less (cnt, io_buffer.size());
                this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
              }

          // Then read the element DOF values
          for (std::vector<const DofObject *>::const_iterator
                 it = ordered_elements.begin(); it != ordered_elements.end(); ++it)
            for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
              {
                libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                             DofObject::invalid_id);
                libmesh_assert_less (cnt, io_buffer.size());
                this->solution->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
              }
        }
    }

  // Finally, read the SCALAR variables on the last processor
  for (unsigned int data_var=0; data_var<nv; data_var++)
    {
      const unsigned int var = _written_var_indices[data_var];
      if(this->variable(var).type().family == SCALAR)
        {
          if (this->processor_id() == (this->n_processors()-1))
            {
              const DofMap & dof_map = this->get_dof_map();
              std::vector<dof_id_type> SCALAR_dofs;
              dof_map.SCALAR_dof_indices(SCALAR_dofs, var);

              for (std::size_t i=0; i<SCALAR_dofs.size(); i++)
                this->solution->set(SCALAR_dofs[i], io_buffer[cnt++]);
            }
        }
    }

  // And we're done setting solution entries
  this->solution->close();

  // For each additional vector, simply go through the list.
  // ONLY attempt to do this IF additional data was actually
  // written to the file for this system (controlled by the
  // _additional_data_written flag).
  if (this->_additional_data_written)
    {
      const std::size_t nvecs = this->_vectors.size();

      // If the number of additional vectors written is non-zero, and
      // the number of additional vectors we have is non-zero, and
      // they don't match, then something is wrong and we can't be
      // sure we're reading data into the correct places.
      if (read_additional_data && nvecs &&
          nvecs != this->_additional_data_written)
        libmesh_error_msg
          ("Additional vectors in file do not match system");

      std::map<std::string, NumericVector<Number> * >::const_iterator
        pos = _vectors.begin();

      for (std::size_t i = 0; i != this->_additional_data_written; ++i)
        {
          cnt=0;
          io_buffer.clear();

          // 10.)
          //
          // Actually read the additional vector components
          // for the ith system from disk
          io.data(io_buffer);

          total_read_size += cast_int<dof_id_type>(io_buffer.size());

          // If read_additional_data==true and we have additional vectors,
          // then we will keep this vector data; otherwise we are going to
          // throw it away.
          if (read_additional_data && nvecs)
            {
              // Loop over each non-SCALAR variable and each node, and read out the value.
              for (unsigned int data_var=0; data_var<nv; data_var++)
                {
                  const unsigned int var = _written_var_indices[data_var];
                  if(this->variable(var).type().family != SCALAR)
                    {
                      // First read the node DOF values
                      for (std::vector<const DofObject *>::const_iterator
                             it = ordered_nodes.begin(); it != ordered_nodes.end(); ++it)
                        for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
                          {
                            libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                                         DofObject::invalid_id);
                            libmesh_assert_less (cnt, io_buffer.size());
                            pos->second->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
                          }

                      // Then read the element DOF values
                      for (std::vector<const DofObject *>::const_iterator
                             it = ordered_elements.begin(); it != ordered_elements.end(); ++it)
                        for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
                          {
                            libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                                         DofObject::invalid_id);
                            libmesh_assert_less (cnt, io_buffer.size());
                            pos->second->set((*it)->dof_number(sys_num, var, comp), io_buffer[cnt++]);
                          }
                    }
                }

              // Finally, read the SCALAR variables on the last processor
              for (unsigned int data_var=0; data_var<nv; data_var++)
                {
                  const unsigned int var = _written_var_indices[data_var];
                  if(this->variable(var).type().family == SCALAR)
                    {
                      if (this->processor_id() == (this->n_processors()-1))
                        {
                          const DofMap & dof_map = this->get_dof_map();
                          std::vector<dof_id_type> SCALAR_dofs;
                          dof_map.SCALAR_dof_indices(SCALAR_dofs, var);

                          for (std::size_t i=0; i<SCALAR_dofs.size(); i++)
                            pos->second->set(SCALAR_dofs[i], io_buffer[cnt++]);
                        }
                    }
                }

              // And we're done setting entries for this variable
              pos->second->close();
            }

          // If we've got vectors then we need to be iterating through
          // those too
          if (pos != this->_vectors.end())
            ++pos;
        }
    }

  // const Real
  //   dt   = pl.get_elapsed_time(),
  //   rate = total_read_size*sizeof(Number)/dt;

  // libMesh::err << "Read " << total_read_size << " \"Number\" values\n"
  //     << " Elapsed time = " << dt << '\n'
  //     << " Rate = " << rate/1.e6 << "(MB/sec)\n\n";

  // pl.pop("read_parallel_data");
}


template <typename InValType>
void System::read_serialized_data (Xdr & io,
                                   const bool read_additional_data)
{
  // This method implements the input of the vectors
  // contained in this System object, embedded in the
  // output of an EquationSystems<T_sys>.
  //
  //   10.) The global solution vector, re-ordered to be node-major
  //       (More on this later.)
  //
  //      for each additional vector in the object
  //
  //      11.) The global additional vector, re-ordered to be
  //          node-major (More on this later.)
  parallel_object_only();
  std::string comment;

  // PerfLog pl("IO Performance",false);
  // pl.push("read_serialized_data");
  // std::size_t total_read_size = 0;

  // 10.)
  // Read the global solution vector
  {
    // total_read_size +=
    this->read_serialized_vector<InValType>(io, this->solution.get());

    // get the comment
    if (this->processor_id() == 0)
      io.comment (comment);
  }

  // 11.)
  // Only read additional vectors if data is available, and only use
  // that data to fill our vectors if the user requested it.
  if (this->_additional_data_written)
    {
      const std::size_t nvecs = this->_vectors.size();

      // If the number of additional vectors written is non-zero, and
      // the number of additional vectors we have is non-zero, and
      // they don't match, then we can't read additional vectors
      // and be sure we're reading data into the correct places.
      if (read_additional_data && nvecs &&
          nvecs != this->_additional_data_written)
        libmesh_error_msg
          ("Additional vectors in file do not match system");

      std::map<std::string, NumericVector<Number> * >::const_iterator
        pos = _vectors.begin();

      for (std::size_t i = 0; i != this->_additional_data_written; ++i)
        {
          // Read data, but only put it into a vector if we've been
          // asked to and if we have a corresponding vector to read.

          // total_read_size +=
          this->read_serialized_vector<InValType>(io,
            (read_additional_data && nvecs) ? pos->second : libmesh_nullptr);

          // get the comment
          if (this->processor_id() == 0)
            io.comment (comment);


          // If we've got vectors then we need to be iterating through
          // those too
          if (pos != this->_vectors.end())
            ++pos;
        }
    }

  // const Real
  //   dt   = pl.get_elapsed_time(),
  //   rate = total_read_size*sizeof(Number)/dt;

  // libMesh::out << "Read " << total_read_size << " \"Number\" values\n"
  //     << " Elapsed time = " << dt << '\n'
  //     << " Rate = " << rate/1.e6 << "(MB/sec)\n\n";

  // pl.pop("read_serialized_data");
}



template <typename iterator_type, typename InValType>
std::size_t System::read_serialized_blocked_dof_objects (const dof_id_type n_objs,
                                                         const iterator_type begin,
                                                         const iterator_type end,
                                                         const InValType ,
                                                         Xdr & io,
                                                         const std::vector<NumericVector<Number> *> & vecs,
                                                         const unsigned int var_to_read) const
{
  //-------------------------------------------------------
  // General order: (IO format 0.7.4 & greater)
  //
  // for (objects ...)
  //   for (vecs ....)
  //     for (vars ....)
  //       for (comps ...)
  //
  // where objects are nodes or elements, sorted to be
  // partition independent,
  // vecs are one or more *identically distributed* solution
  // coefficient vectors, vars are one or more variables
  // to write, and comps are all the components for said
  // vars on the object.

  typedef std::vector<NumericVector<Number> *>::const_iterator vec_iterator_type;

  // variables to read.  Unless specified otherwise, defaults to _written_var_indices.
  std::vector<unsigned int> vars_to_read (_written_var_indices);

  if (var_to_read != libMesh::invalid_uint)
    {
      vars_to_read.clear();
      vars_to_read.push_back(var_to_read);
    }

  const unsigned int
    sys_num    = this->number(),
    num_vecs   = cast_int<unsigned int>(vecs.size());
  const dof_id_type
    io_blksize = cast_int<dof_id_type>(std::min(max_io_blksize, static_cast<std::size_t>(n_objs))),
    num_blks   = cast_int<unsigned int>(std::ceil(static_cast<double>(n_objs)/
                                                  static_cast<double>(io_blksize)));

  libmesh_assert_less_equal (_written_var_indices.size(), this->n_vars());

  std::size_t n_read_values=0;

  std::vector<std::vector<dof_id_type> > xfer_ids(num_blks);  // The global IDs and # of components for the local objects in all blocks
  std::vector<std::vector<Number> >      recv_vals(num_blks); // The raw values for the local objects in all blocks
  std::vector<Parallel::Request>
    id_requests(num_blks), val_requests(num_blks);

  // ------------------------------------------------------
  // First pass - count the number of objects in each block
  // traverse all the objects and figure out which block they
  // will utlimately live in.
  std::vector<std::size_t>
    xfer_ids_size  (num_blks,0),
    recv_vals_size (num_blks,0);


  for (iterator_type it=begin; it!=end; ++it)
    {
      const dof_id_type
        id    = (*it)->id(),
        block = id/io_blksize;

      libmesh_assert_less (block, num_blks);

      xfer_ids_size[block] += 2; // for each object, we send its id, as well as the total number of components for all variables

      dof_id_type n_comp_tot=0;
      for (std::vector<unsigned int>::const_iterator var_it=vars_to_read.begin();
           var_it!=vars_to_read.end(); ++var_it)
        n_comp_tot += (*it)->n_comp(sys_num, *var_it); // for each variable, we will receive the nonzero components

      recv_vals_size[block] += n_comp_tot*num_vecs;
    }

  // knowing the recv_vals_size[block] for each processor allows
  // us to sum them and find the global size for each block.
  std::vector<std::size_t> tot_vals_size(recv_vals_size);
  this->comm().sum (tot_vals_size);


  //------------------------------------------
  // Collect the ids & number of values needed
  // for all local objects, binning them into
  // 'blocks' that will be sent to processor 0
  for (dof_id_type blk=0; blk<num_blks; blk++)
    {
      // Each processor should build up its transfer buffers for its
      // local objects in [first_object,last_object).
      const dof_id_type
        first_object = blk*io_blksize,
        last_object  = std::min(cast_int<dof_id_type>((blk+1)*io_blksize), n_objs);

      // convenience
      std::vector<dof_id_type> & ids (xfer_ids[blk]);
      std::vector<Number> & vals (recv_vals[blk]);

      // we now know the number of values we will store for each block,
      // so we can do efficient preallocation
      ids.clear(); /**/ ids.reserve (xfer_ids_size[blk]);
      vals.resize(recv_vals_size[blk]);

      if (recv_vals_size[blk] != 0) // only if there are nonzero values to receive
        for (iterator_type it=begin; it!=end; ++it)
          if (((*it)->id() >= first_object) && // object in [first_object,last_object)
              ((*it)->id() <   last_object))
            {
              ids.push_back((*it)->id());

              unsigned int n_comp_tot=0;

              for (std::vector<unsigned int>::const_iterator var_it=vars_to_read.begin();
                   var_it!=vars_to_read.end(); ++var_it)
                n_comp_tot += (*it)->n_comp(sys_num,*var_it);

              ids.push_back (n_comp_tot*num_vecs);
            }

#ifdef LIBMESH_HAVE_MPI
      Parallel::MessageTag id_tag  = this->comm().get_unique_tag(100*num_blks + blk);
      Parallel::MessageTag val_tag = this->comm().get_unique_tag(200*num_blks + blk);

      // nonblocking send the data for this block
      this->comm().send (0, ids,  id_requests[blk],  id_tag);

      // Go ahead and post the receive too
      this->comm().receive (0, vals, val_requests[blk], val_tag);
#endif
    }

  //---------------------------------------------------
  // Here processor 0 will read and distribute the data.
  // We have to do this block-wise to ensure that we
  // do not exhaust memory on processor 0.

  // give these variables scope outside the block to avoid reallocation
  std::vector<std::vector<dof_id_type> > recv_ids       (this->n_processors());
  std::vector<std::vector<Number> >      send_vals      (this->n_processors());
  std::vector<Parallel::Request>         reply_requests (this->n_processors());
  std::vector<unsigned int>              obj_val_offsets;          // map to traverse entry-wise rather than processor-wise
  std::vector<Number>                    input_vals;               // The input buffer for the current block
  std::vector<InValType>                 input_vals_tmp;               // The input buffer for the current block

  for (dof_id_type blk=0; blk<num_blks; blk++)
    {
      // Each processor should build up its transfer buffers for its
      // local objects in [first_object,last_object).
      const dof_id_type
        first_object  = blk*io_blksize,
        last_object   = std::min(cast_int<dof_id_type>((blk+1)*io_blksize), n_objs),
        n_objects_blk = last_object - first_object;

      // Processor 0 has a special job.  It needs to gather the requested indices
      // in [first_object,last_object) from all processors, read the data from
      // disk, and reply
      if (this->processor_id() == 0)
        {
          // we know the input buffer size for this block and can begin reading it now
          input_vals.resize(tot_vals_size[blk]);
          input_vals_tmp.resize(tot_vals_size[blk]);

          // a ThreadedIO object to perform asychronous file IO
          ThreadedIO<InValType> threaded_io(io, input_vals_tmp);
          Threads::Thread async_io(threaded_io);

          Parallel::MessageTag id_tag  = this->comm().get_unique_tag(100*num_blks + blk);
          Parallel::MessageTag val_tag = this->comm().get_unique_tag(200*num_blks + blk);

          // offset array. this will define where each object's values
          // map into the actual input_vals buffer.  this must get
          // 0-initialized because 0-component objects are not actually sent
          obj_val_offsets.resize (n_objects_blk); /**/ std::fill (obj_val_offsets.begin(), obj_val_offsets.end(), 0);
          recv_vals_size.resize(this->n_processors()); // reuse this to count how many values are going to each processor

#ifndef NDEBUG
          std::size_t n_vals_blk = 0;
#endif

          // loop over all processors and process their index request
          for (unsigned int comm_step=0; comm_step<this->n_processors(); comm_step++)
            {
#ifdef LIBMESH_HAVE_MPI
              // blocking receive indices for this block, imposing no particular order on processor
              Parallel::Status id_status (this->comm().probe (Parallel::any_source, id_tag));
              std::vector<dof_id_type> & ids (recv_ids[id_status.source()]);
              std::size_t & n_vals_proc (recv_vals_size[id_status.source()]);
              this->comm().receive (id_status.source(), ids, id_tag);
#else
              // straight copy without MPI
              std::vector<dof_id_type> & ids (recv_ids[0]);
              std::size_t & n_vals_proc (recv_vals_size[0]);
              ids = xfer_ids[blk];
#endif

              n_vals_proc = 0;

              // note its possible we didn't receive values for objects in
              // this block if they have no components allocated.
              for (std::size_t idx=0; idx<ids.size(); idx+=2)
                {
                  const dof_id_type
                    local_idx          = ids[idx+0]-first_object,
                    n_vals_tot_allvecs = ids[idx+1];

                  libmesh_assert_less (local_idx, n_objects_blk);

                  obj_val_offsets[local_idx] = n_vals_tot_allvecs;
                  n_vals_proc += n_vals_tot_allvecs;
                }

#ifndef NDEBUG
              n_vals_blk += n_vals_proc;
#endif
            }

          // We need the offests into the input_vals vector for each object.
          // fortunately, this is simply the partial sum of the total number
          // of components for each object
          std::partial_sum(obj_val_offsets.begin(), obj_val_offsets.end(),
                           obj_val_offsets.begin());

          libmesh_assert_equal_to (n_vals_blk, obj_val_offsets.back());
          libmesh_assert_equal_to (n_vals_blk, tot_vals_size[blk]);

          // Wait for read completion
          async_io.join();
          // now copy the values back to the main vector for transfer
          for (std::size_t i_val=0; i_val<input_vals.size(); i_val++)
            input_vals[i_val] = input_vals_tmp[i_val];

          n_read_values += input_vals.size();

          // pack data replies for each processor
          for (processor_id_type proc=0; proc<this->n_processors(); proc++)
            {
              const std::vector<dof_id_type> & ids (recv_ids[proc]);
              std::vector<Number> & vals (send_vals[proc]);
              const std::size_t & n_vals_proc (recv_vals_size[proc]);

              vals.clear(); /**/ vals.reserve(n_vals_proc);

              for (std::size_t idx=0; idx<ids.size(); idx+=2)
                {
                  const dof_id_type
                    local_idx          = ids[idx+0]-first_object,
                    n_vals_tot_allvecs = ids[idx+1];

                  std::vector<Number>::const_iterator in_vals(input_vals.begin());
                  if (local_idx != 0)
                    std::advance (in_vals, obj_val_offsets[local_idx-1]);

                  for (unsigned int val=0; val<n_vals_tot_allvecs; val++, ++in_vals)
                    {
                      libmesh_assert (in_vals != input_vals.end());
                      //libMesh::out << "*in_vals=" << *in_vals << '\n';
                      vals.push_back(*in_vals);
                    }
                }

#ifdef LIBMESH_HAVE_MPI
              // send the relevant values to this processor
              this->comm().send (proc, vals, reply_requests[proc], val_tag);
#else
              recv_vals[blk] = vals;
#endif
            }
        } // end processor 0 read/reply

      // all processors complete the (already posted) read for this block
      {
        Parallel::wait (val_requests[blk]);

        const std::vector<Number> & vals (recv_vals[blk]);
        std::vector<Number>::const_iterator val_it(vals.begin());

        if (!recv_vals[blk].empty()) // nonzero values to receive
          for (iterator_type it=begin; it!=end; ++it)
            if (((*it)->id() >= first_object) && // object in [first_object,last_object)
                ((*it)->id() <   last_object))
              // unpack & set the values
              for (vec_iterator_type vec_it=vecs.begin(); vec_it!=vecs.end(); ++vec_it)
                {
                  NumericVector<Number> * vec(*vec_it);

                  for (std::vector<unsigned int>::const_iterator var_it=vars_to_read.begin();
                       var_it!=vars_to_read.end(); ++var_it)
                    {
                      const unsigned int n_comp = (*it)->n_comp(sys_num,*var_it);

                      for (unsigned int comp=0; comp<n_comp; comp++, ++val_it)
                        {
                          const dof_id_type dof_index = (*it)->dof_number (sys_num, *var_it, comp);
                          libmesh_assert (val_it != vals.end());
                          if (vec)
                            {
                              libmesh_assert_greater_equal (dof_index, vec->first_local_index());
                              libmesh_assert_less (dof_index, vec->last_local_index());
                              //libMesh::out << "dof_index, *val_it = \t" << dof_index << ", " << *val_it << '\n';
                              vec->set (dof_index, *val_it);
                            }
                        }
                    }
                }
      }

      // processor 0 needs to make sure all replies have been handed off
      if (this->processor_id () == 0)
        Parallel::wait(reply_requests);
    }

  return n_read_values;
}



unsigned int System::read_SCALAR_dofs (const unsigned int var,
                                       Xdr & io,
                                       NumericVector<Number> * vec) const
{
  unsigned int n_assigned_vals = 0; // the number of values assigned, this will be returned.

  // Processor 0 will read the block from the buffer stream and send it to the last processor
  const unsigned int n_SCALAR_dofs = this->variable(var).type().order.get_order();
  std::vector<Number> input_buffer(n_SCALAR_dofs);
  if (this->processor_id() == 0)
    {
      io.data_stream(&input_buffer[0], n_SCALAR_dofs);
    }

#ifdef LIBMESH_HAVE_MPI
  if ( this->n_processors() > 1 )
    {
      const Parallel::MessageTag val_tag = this->comm().get_unique_tag(321);

      // Post the receive on the last processor
      if (this->processor_id() == (this->n_processors()-1))
        this->comm().receive(0, input_buffer, val_tag);

      // Send the data to processor 0
      if (this->processor_id() == 0)
        this->comm().send(this->n_processors()-1, input_buffer, val_tag);
    }
#endif

  // Finally, set the SCALAR values
  if (this->processor_id() == (this->n_processors()-1))
    {
      const DofMap & dof_map = this->get_dof_map();
      std::vector<dof_id_type> SCALAR_dofs;
      dof_map.SCALAR_dof_indices(SCALAR_dofs, var);

      for (std::size_t i=0; i<SCALAR_dofs.size(); i++)
        {
          if (vec)
            vec->set (SCALAR_dofs[i], input_buffer[i]);
          ++n_assigned_vals;
        }
    }

  return n_assigned_vals;
}


template <typename InValType>
numeric_index_type System::read_serialized_vector (Xdr & io,
                                                   NumericVector<Number> * vec)
{
  parallel_object_only();

#ifndef NDEBUG
  // In parallel we better be reading a parallel vector -- if not
  // we will not set all of its components below!!
  if (this->n_processors() > 1 && vec)
    {
      libmesh_assert (vec->type() == PARALLEL ||
                      vec->type() == GHOSTED);
    }
#endif

  libmesh_assert (io.reading());

  // vector length
  unsigned int vector_length=0; // FIXME?  size_t would break binary compatibility...
#ifndef NDEBUG
  std::size_t n_assigned_vals=0;
#endif

  // Get the buffer size
  if (this->processor_id() == 0)
    io.data(vector_length, "# vector length");
  this->comm().broadcast(vector_length);

  const unsigned int nv = cast_int<unsigned int>
    (this->_written_var_indices.size());
  const dof_id_type
    n_nodes = this->get_mesh().n_nodes(),
    n_elem  = this->get_mesh().n_elem();

  libmesh_assert_less_equal (nv, this->n_vars());

  // for newer versions, read variables node/elem major
  if (io.version() >= LIBMESH_VERSION_ID(0,7,4))
    {
      //---------------------------------
      // Collect the values for all nodes
#ifndef NDEBUG
      n_assigned_vals +=
#endif
        this->read_serialized_blocked_dof_objects (n_nodes,
                                                   this->get_mesh().local_nodes_begin(),
                                                   this->get_mesh().local_nodes_end(),
                                                   InValType(),
                                                   io,
                                                   std::vector<NumericVector<Number> *> (1,vec));


      //------------------------------------
      // Collect the values for all elements
#ifndef NDEBUG
      n_assigned_vals +=
#endif
        this->read_serialized_blocked_dof_objects (n_elem,
                                                   this->get_mesh().local_elements_begin(),
                                                   this->get_mesh().local_elements_end(),
                                                   InValType(),
                                                   io,
                                                   std::vector<NumericVector<Number> *> (1,vec));
    }

  // for older versions, read variables var-major
  else
    {
      // Loop over each variable in the system, and then each node/element in the mesh.
      for (unsigned int data_var=0; data_var<nv; data_var++)
        {
          const unsigned int var = _written_var_indices[data_var];
          if(this->variable(var).type().family != SCALAR)
            {
              //---------------------------------
              // Collect the values for all nodes
#ifndef NDEBUG
              n_assigned_vals +=
#endif
                this->read_serialized_blocked_dof_objects (n_nodes,
                                                           this->get_mesh().local_nodes_begin(),
                                                           this->get_mesh().local_nodes_end(),
                                                           InValType(),
                                                           io,
                                                           std::vector<NumericVector<Number> *> (1,vec),
                                                           var);


              //------------------------------------
              // Collect the values for all elements
#ifndef NDEBUG
              n_assigned_vals +=
#endif
                this->read_serialized_blocked_dof_objects (n_elem,
                                                           this->get_mesh().local_elements_begin(),
                                                           this->get_mesh().local_elements_end(),
                                                           InValType(),
                                                           io,
                                                           std::vector<NumericVector<Number> *> (1,vec),
                                                           var);
            } // end variable loop
        }
    }

  //-------------------------------------------
  // Finally loop over all the SCALAR variables
  for (unsigned int data_var=0; data_var<nv; data_var++)
    {
      const unsigned int var = _written_var_indices[data_var];
      if(this->variable(var).type().family == SCALAR)
        {
#ifndef NDEBUG
          n_assigned_vals +=
#endif
            this->read_SCALAR_dofs (var, io, vec);
        }
    }

  if (vec)
    vec->close();

#ifndef NDEBUG
  this->comm().sum (n_assigned_vals);
  libmesh_assert_equal_to (n_assigned_vals, vector_length);
#endif

  return vector_length;
}



void System::write_header (Xdr & io,
                           const std::string & /* version is currently unused */,
                           const bool write_additional_data) const
{
  /**
   * This method implements the output of a
   * System object, embedded in the output of
   * an EquationSystems<T_sys>.  This warrants some
   * documentation.  The output of this part
   * consists of 5 sections:
   *
   * for this system
   *
   *   5.) The number of variables in the system (unsigned int)
   *
   *   for each variable in the system
   *
   *     6.) The name of the variable (string)
   *
   *     6.1.) subdomain where the variable lives
   *
   *     7.) Combined in an FEType:
   *         - The approximation order(s) of the variable
   *           (Order Enum, cast to int/s)
   *         - The finite element family/ies of the variable
   *           (FEFamily Enum, cast to int/s)
   *
   *   end variable loop
   *
   *   8.) The number of additional vectors (unsigned int),
   *
   *     for each additional vector in the system object
   *
   *     9.) the name of the additional vector  (string)
   *
   * end system
   */
  libmesh_assert (io.writing());


  // Only write the header information
  // if we are processor 0.
  if (this->get_mesh().processor_id() != 0)
    return;

  std::string comment;
  char buf[80];

  // 5.)
  // Write the number of variables in the system

  {
    // set up the comment
    comment = "# No. of Variables in System \"";
    comment += this->name();
    comment += "\"";

    unsigned int nv = this->n_vars();
    io.data (nv, comment.c_str());
  }


  for (unsigned int var=0; var<this->n_vars(); var++)
    {
      // 6.)
      // Write the name of the var-th variable
      {
        // set up the comment
        comment  = "#   Name, Variable No. ";
        std::sprintf(buf, "%u", var);
        comment += buf;
        comment += ", System \"";
        comment += this->name();
        comment += "\"";

        std::string var_name = this->variable_name(var);
        io.data (var_name, comment.c_str());
      }

      // 6.1.) Variable subdomains
      {
        // set up the comment
        comment  = "#     Subdomains, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        const std::set<subdomain_id_type> & domains = this->variable(var).active_subdomains();
        std::vector<subdomain_id_type> domain_array;
        domain_array.assign(domains.begin(), domains.end());
        io.data (domain_array, comment.c_str());
      }

      // 7.)
      // Write the approximation order of the var-th variable
      // in this system
      {
        // set up the comment
        comment = "#     Approximation Order, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        int order = static_cast<int>(this->variable_type(var).order);
        io.data (order, comment.c_str());
      }


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

      // do the same for radial_order
      {
        comment = "#     Radial Approximation Order, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        int rad_order = static_cast<int>(this->variable_type(var).radial_order);
        io.data (rad_order, comment.c_str());
      }

#endif

      // Write the Finite Element type of the var-th variable
      // in this System
      {
        // set up the comment
        comment = "#     FE Family, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        const FEType & type = this->variable_type(var);
        int fam = static_cast<int>(type.family);
        io.data (fam, comment.c_str());

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

        comment = "#     Radial FE Family, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        int radial_fam = static_cast<int>(type.radial_family);
        io.data (radial_fam, comment.c_str());

        comment = "#     Infinite Mapping Type, Variable \"";
        std::sprintf(buf, "%s", this->variable_name(var).c_str());
        comment += buf;
        comment += "\", System \"";
        comment += this->name();
        comment += "\"";

        int i_map = static_cast<int>(type.inf_map);
        io.data (i_map, comment.c_str());
#endif
      }
    } // end of the variable loop

  // 8.)
  // Write the number of additional vectors in the System.
  // If write_additional_data==false, then write zero for
  // the number of additional vectors.
  {
    {
      // set up the comment
      comment = "# No. of Additional Vectors, System \"";
      comment += this->name();
      comment += "\"";

      unsigned int nvecs = write_additional_data ? this->n_vectors () : 0;
      io.data (nvecs, comment.c_str());
    }

    if (write_additional_data)
      {
        std::map<std::string, NumericVector<Number> * >::const_iterator
          vec_pos = this->_vectors.begin();
        unsigned int cnt=0;

        for (; vec_pos != this->_vectors.end(); ++vec_pos)
          {
            // 9.)
            // write the name of the cnt-th additional vector
            comment =  "# Name of ";
            std::sprintf(buf, "%d", cnt++);
            comment += buf;
            comment += "th vector";
            std::string vec_name = vec_pos->first;

            io.data (vec_name, comment.c_str());
          }
      }
  }
}



void System::write_parallel_data (Xdr & io,
                                  const bool write_additional_data) const
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the
   * output of an EquationSystems<T_sys>.
   *
   *   9.) The global solution vector, re-ordered to be node-major
   *       (More on this later.)
   *
   *      for each additional vector in the object
   *
   *      10.) The global additional vector, re-ordered to be
   *           node-major (More on this later.)
   *
   * Note that the actual IO is handled through the Xdr class
   * (to be renamed later?) which provides a uniform interface to
   * both the XDR (eXternal Data Representation) interface and standard
   * ASCII output.  Thus this one section of code will read XDR or ASCII
   * files with no changes.
   */
  // PerfLog pl("IO Performance",false);
  // pl.push("write_parallel_data");
  // std::size_t total_written_size = 0;

  std::string comment;

  libmesh_assert (io.writing());

  std::vector<Number> io_buffer; io_buffer.reserve(this->solution->local_size());

  // build the ordered nodes and element maps.
  // when writing/reading parallel files we need to iterate
  // over our nodes/elements in order of increasing global id().
  // however, this is not guaranteed to be ordering we obtain
  // by using the node_iterators/element_iterators directly.
  // so build a set, sorted by id(), that provides the ordering.
  // further, for memory economy build the set but then transfer
  // its contents to vectors, which will be sorted.
  std::vector<const DofObject *> ordered_nodes, ordered_elements;
  {
    std::set<const DofObject *, CompareDofObjectsByID>
      ordered_nodes_set (this->get_mesh().local_nodes_begin(),
                         this->get_mesh().local_nodes_end());

    ordered_nodes.insert(ordered_nodes.end(),
                         ordered_nodes_set.begin(),
                         ordered_nodes_set.end());
  }
  {
    std::set<const DofObject *, CompareDofObjectsByID>
      ordered_elements_set (this->get_mesh().local_elements_begin(),
                            this->get_mesh().local_elements_end());

    ordered_elements.insert(ordered_elements.end(),
                            ordered_elements_set.begin(),
                            ordered_elements_set.end());
  }

  const unsigned int sys_num = this->number();
  const unsigned int nv      = this->n_vars();

  // Loop over each non-SCALAR variable and each node, and write out the value.
  for (unsigned int var=0; var<nv; var++)
    if (this->variable(var).type().family != SCALAR)
      {
        // First write the node DOF values
        for (std::vector<const DofObject *>::const_iterator
               it = ordered_nodes.begin(); it != ordered_nodes.end(); ++it)
          for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
            {
              //libMesh::out << "(*it)->id()=" << (*it)->id() << std::endl;
              libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                           DofObject::invalid_id);

              io_buffer.push_back((*this->solution)((*it)->dof_number(sys_num, var, comp)));
            }

        // Then write the element DOF values
        for (std::vector<const DofObject *>::const_iterator
               it = ordered_elements.begin(); it != ordered_elements.end(); ++it)
          for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
            {
              libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                           DofObject::invalid_id);

              io_buffer.push_back((*this->solution)((*it)->dof_number(sys_num, var, comp)));
            }
      }

  // Finally, write the SCALAR data on the last processor
  for (unsigned int var=0; var<this->n_vars(); var++)
    if(this->variable(var).type().family == SCALAR)
      {
        if (this->processor_id() == (this->n_processors()-1))
          {
            const DofMap & dof_map = this->get_dof_map();
            std::vector<dof_id_type> SCALAR_dofs;
            dof_map.SCALAR_dof_indices(SCALAR_dofs, var);

            for (std::size_t i=0; i<SCALAR_dofs.size(); i++)
              io_buffer.push_back((*this->solution)(SCALAR_dofs[i]));
          }
      }

  // 9.)
  //
  // Actually write the reordered solution vector
  // for the ith system to disk

  // set up the comment
  {
    comment = "# System \"";
    comment += this->name();
    comment += "\" Solution Vector";
  }

  io.data (io_buffer, comment.c_str());

  // total_written_size += io_buffer.size();

  // Only write additional vectors if wanted
  if (write_additional_data)
    {
      std::map<std::string, NumericVector<Number> *>::const_iterator
        pos = _vectors.begin();

      for(; pos != this->_vectors.end(); ++pos)
        {
          io_buffer.clear(); io_buffer.reserve( pos->second->local_size());

          // Loop over each non-SCALAR variable and each node, and write out the value.
          for (unsigned int var=0; var<nv; var++)
            if(this->variable(var).type().family != SCALAR)
              {
                // First write the node DOF values
                for (std::vector<const DofObject *>::const_iterator
                       it = ordered_nodes.begin(); it != ordered_nodes.end(); ++it)
                  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
                    {
                      libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                                   DofObject::invalid_id);

                      io_buffer.push_back((*pos->second)((*it)->dof_number(sys_num, var, comp)));
                    }

                // Then write the element DOF values
                for (std::vector<const DofObject *>::const_iterator
                       it = ordered_elements.begin(); it != ordered_elements.end(); ++it)
                  for (unsigned int comp=0; comp<(*it)->n_comp(sys_num, var); comp++)
                    {
                      libmesh_assert_not_equal_to ((*it)->dof_number(sys_num, var, comp),
                                                   DofObject::invalid_id);

                      io_buffer.push_back((*pos->second)((*it)->dof_number(sys_num, var, comp)));
                    }
              }

          // Finally, write the SCALAR data on the last processor
          for (unsigned int var=0; var<this->n_vars(); var++)
            if(this->variable(var).type().family == SCALAR)
              {
                if (this->processor_id() == (this->n_processors()-1))
                  {
                    const DofMap & dof_map = this->get_dof_map();
                    std::vector<dof_id_type> SCALAR_dofs;
                    dof_map.SCALAR_dof_indices(SCALAR_dofs, var);

                    for (std::size_t i=0; i<SCALAR_dofs.size(); i++)
                      io_buffer.push_back((*pos->second)(SCALAR_dofs[i]));
                  }
              }

          // 10.)
          //
          // Actually write the reordered additional vector
          // for this system to disk

          // set up the comment
          {
            comment = "# System \"";
            comment += this->name();
            comment += "\" Additional Vector \"";
            comment += pos->first;
            comment += "\"";
          }

          io.data (io_buffer, comment.c_str());

          // total_written_size += io_buffer.size();
        }
    }

  // const Real
  //   dt   = pl.get_elapsed_time(),
  //   rate = total_written_size*sizeof(Number)/dt;

  // libMesh::err << "Write " << total_written_size << " \"Number\" values\n"
  //     << " Elapsed time = " << dt << '\n'
  //     << " Rate = " << rate/1.e6 << "(MB/sec)\n\n";

  // pl.pop("write_parallel_data");
}



void System::write_serialized_data (Xdr & io,
                                    const bool write_additional_data) const
{
  /**
   * This method implements the output of the vectors
   * contained in this System object, embedded in the
   * output of an EquationSystems<T_sys>.
   *
   *   9.) The global solution vector, re-ordered to be node-major
   *       (More on this later.)
   *
   *      for each additional vector in the object
   *
   *      10.) The global additional vector, re-ordered to be
   *          node-major (More on this later.)
   */
  parallel_object_only();
  std::string comment;

  // PerfLog pl("IO Performance",false);
  // pl.push("write_serialized_data");
  // std::size_t total_written_size = 0;

  // total_written_size +=
  this->write_serialized_vector(io, *this->solution);

  // set up the comment
  if (this->processor_id() == 0)
    {
      comment = "# System \"";
      comment += this->name();
      comment += "\" Solution Vector";

      io.comment (comment);
    }

  // Only write additional vectors if wanted
  if (write_additional_data)
    {
      std::map<std::string, NumericVector<Number> *>::const_iterator
        pos = _vectors.begin();

      for(; pos != this->_vectors.end(); ++pos)
        {
          // total_written_size +=
          this->write_serialized_vector(io, *pos->second);

          // set up the comment
          if (this->processor_id() == 0)
            {
              comment = "# System \"";
              comment += this->name();
              comment += "\" Additional Vector \"";
              comment += pos->first;
              comment += "\"";
              io.comment (comment);
            }
        }
    }

  // const Real
  //   dt   = pl.get_elapsed_time(),
  //   rate = total_written_size*sizeof(Number)/dt;

  // libMesh::out << "Write " << total_written_size << " \"Number\" values\n"
  //     << " Elapsed time = " << dt << '\n'
  //     << " Rate = " << rate/1.e6 << "(MB/sec)\n\n";

  // pl.pop("write_serialized_data");




  // // test the new method
  // {
  //   std::vector<std::string> names;
  //   std::vector<NumericVector<Number> *> vectors_to_write;

  //   names.push_back("Solution Vector");
  //   vectors_to_write.push_back(this->solution.get());

  //   // Only write additional vectors if wanted
  //   if (write_additional_data)
  //     {
  // std::map<std::string, NumericVector<Number> * >::const_iterator
  //   pos = _vectors.begin();

  // for(; pos != this->_vectors.end(); ++pos)
  //   {
  //     names.push_back("Additional Vector " + pos->first);
  //     vectors_to_write.push_back(pos->second);
  //   }
  //     }

  //   total_written_size =
  //     this->write_serialized_vectors (io, names, vectors_to_write);

  //   const Real
  //     dt2   = pl.get_elapsed_time(),
  //     rate2 = total_written_size*sizeof(Number)/(dt2-dt);

  //   libMesh::out << "Write (new) " << total_written_size << " \"Number\" values\n"
  //       << " Elapsed time = " << (dt2-dt) << '\n'
  //       << " Rate = " << rate2/1.e6 << "(MB/sec)\n\n";

  // }
}



template <typename iterator_type>
std::size_t System::write_serialized_blocked_dof_objects (const std::vector<const NumericVector<Number> *> & vecs,
                                                          const dof_id_type n_objs,
                                                          const iterator_type begin,
                                                          const iterator_type end,
                                                          Xdr & io,
                                                          const unsigned int var_to_write) const
{
  //-------------------------------------------------------
  // General order: (IO format 0.7.4 & greater)
  //
  // for (objects ...)
  //   for (vecs ....)
  //     for (vars ....)
  //       for (comps ...)
  //
  // where objects are nodes or elements, sorted to be
  // partition independent,
  // vecs are one or more *identically distributed* solution
  // coefficient vectors, vars are one or more variables
  // to write, and comps are all the components for said
  // vars on the object.

  typedef std::vector<const NumericVector<Number> *>::const_iterator vec_iterator_type;

  // We will write all variables unless requested otherwise.
  std::vector<unsigned int> vars_to_write(1, var_to_write);

  if (var_to_write == libMesh::invalid_uint)
    {
      vars_to_write.clear(); /**/ vars_to_write.reserve(this->n_vars());
      for (unsigned int var=0; var<this->n_vars(); var++)
        vars_to_write.push_back(var);
    }

  const dof_id_type io_blksize = cast_int<dof_id_type>
    (std::min(max_io_blksize, static_cast<std::size_t>(n_objs)));

  const unsigned int
    sys_num  = this->number(),
    num_vecs = cast_int<unsigned int>(vecs.size()),
    num_blks = cast_int<unsigned int>(std::ceil(static_cast<double>(n_objs)/
                                                static_cast<double>(io_blksize)));

  // libMesh::out << "io_blksize = "    << io_blksize
  //     << ", num_objects = " << n_objs
  //     << ", num_blks = "    << num_blks
  //     << std::endl;

  dof_id_type written_length=0;                                   // The numer of values written.  This will be returned
  std::vector<std::vector<dof_id_type> > xfer_ids(num_blks);      // The global IDs and # of components for the local objects in all blocks
  std::vector<std::vector<Number> >      send_vals(num_blks);     // The raw values for the local objects in all blocks
  std::vector<Parallel::Request>
    id_requests(num_blks), val_requests(num_blks);                 // send request handle for each block

  // ------------------------------------------------------
  // First pass - count the number of objects in each block
  // traverse all the objects and figure out which block they
  // will utlimately live in.
  std::vector<unsigned int>
    xfer_ids_size  (num_blks,0),
    send_vals_size (num_blks,0);

  for (iterator_type it=begin; it!=end; ++it)
    {
      const dof_id_type
        id    = (*it)->id(),
        block = id/io_blksize;

      libmesh_assert_less (block, num_blks);

      xfer_ids_size[block] += 2; // for each object, we store its id, as well as the total number of components for all variables

      unsigned int n_comp_tot=0;

      for (std::vector<unsigned int>::const_iterator var_it=vars_to_write.begin();
           var_it!=vars_to_write.end(); ++var_it)
        n_comp_tot += (*it)->n_comp(sys_num, *var_it); // for each variable, we will store the nonzero components

      send_vals_size[block] += n_comp_tot*num_vecs;
    }

  //-----------------------------------------
  // Collect the values for all local objects,
  // binning them into 'blocks' that will be
  // sent to processor 0
  for (unsigned int blk=0; blk<num_blks; blk++)
    {
      // libMesh::out << "Writing object block " << blk << std::endl;

      // Each processor should build up its transfer buffers for its
      // local objects in [first_object,last_object).
      const dof_id_type
        first_object = blk*io_blksize,
        last_object  = std::min(cast_int<dof_id_type>((blk+1)*io_blksize), n_objs);

      // convenience
      std::vector<dof_id_type> & ids  (xfer_ids[blk]);
      std::vector<Number>      & vals (send_vals[blk]);

      // we now know the number of values we will store for each block,
      // so we can do efficient preallocation
      ids.clear();  /**/ ids.reserve  (xfer_ids_size[blk]);
      vals.clear(); /**/ vals.reserve (send_vals_size[blk]);

      if (send_vals_size[blk] != 0) // only send if we have nonzero components to write
        for (iterator_type it=begin; it!=end; ++it)
          if (((*it)->id() >= first_object) && // object in [first_object,last_object)
              ((*it)->id() <   last_object))
            {
              ids.push_back((*it)->id());

              // count the total number of nonzeros transferred for this object
              {
                unsigned int n_comp_tot=0;

                for (std::vector<unsigned int>::const_iterator var_it=vars_to_write.begin();
                     var_it!=vars_to_write.end(); ++var_it)
                  n_comp_tot += (*it)->n_comp(sys_num,*var_it);

                ids.push_back (n_comp_tot*num_vecs); // even if 0 - processor 0 has no way of knowing otherwise...
              }

              // pack the values to send
              for (vec_iterator_type vec_it=vecs.begin(); vec_it!=vecs.end(); ++vec_it)
                {
                  const NumericVector<Number> & vec(**vec_it);

                  for (std::vector<unsigned int>::const_iterator var_it=vars_to_write.begin();
                       var_it!=vars_to_write.end(); ++var_it)
                    {
                      const unsigned int n_comp = (*it)->n_comp(sys_num,*var_it);

                      for (unsigned int comp=0; comp<n_comp; comp++)
                        {
                          libmesh_assert_greater_equal ((*it)->dof_number(sys_num, *var_it, comp), vec.first_local_index());
                          libmesh_assert_less ((*it)->dof_number(sys_num, *var_it, comp), vec.last_local_index());
                          vals.push_back(vec((*it)->dof_number(sys_num, *var_it, comp)));
                        }
                    }
                }
            }

#ifdef LIBMESH_HAVE_MPI
      Parallel::MessageTag id_tag  = this->comm().get_unique_tag(100*num_blks + blk);
      Parallel::MessageTag val_tag = this->comm().get_unique_tag(200*num_blks + blk);

      // nonblocking send the data for this block
      this->comm().send (0, ids,  id_requests[blk],  id_tag);
      this->comm().send (0, vals, val_requests[blk], val_tag);
#endif
    }


  if (this->processor_id() == 0)
    {
      std::vector<std::vector<dof_id_type> > recv_ids  (this->n_processors());
      std::vector<std::vector<Number> >      recv_vals (this->n_processors());
      std::vector<unsigned int> obj_val_offsets;          // map to traverse entry-wise rather than processor-wise
      std::vector<Number>       output_vals;              // The output buffer for the current block

      // a ThreadedIO object to perform asychronous file IO
      ThreadedIO<Number> threaded_io(io, output_vals);
      UniquePtr<Threads::Thread> async_io;

      for (unsigned int blk=0; blk<num_blks; blk++)
        {
          // Each processor should build up its transfer buffers for its
          // local objects in [first_object,last_object).
          const dof_id_type
            first_object  = cast_int<dof_id_type>(blk*io_blksize),
            last_object   = std::min(cast_int<dof_id_type>((blk+1)*io_blksize), n_objs),
            n_objects_blk = last_object - first_object;

          // offset array. this will define where each object's values
          // map into the actual output_vals buffer.  this must get
          // 0-initialized because 0-component objects are not actually sent
          obj_val_offsets.resize (n_objects_blk); /**/ std::fill (obj_val_offsets.begin(), obj_val_offsets.end(), 0);

          std::size_t n_val_recvd_blk=0;

          // tags to select data received
          Parallel::MessageTag id_tag  (this->comm().get_unique_tag(100*num_blks + blk));
          Parallel::MessageTag val_tag (this->comm().get_unique_tag(200*num_blks + blk));

          // receive this block of data from all processors.
          for (unsigned int comm_step=0; comm_step<this->n_processors(); comm_step++)
            {
#ifdef LIBMESH_HAVE_MPI
              // blocking receive indices for this block, imposing no particular order on processor
              Parallel::Status id_status (this->comm().probe (Parallel::any_source, id_tag));
              std::vector<dof_id_type> & ids (recv_ids[id_status.source()]);
              this->comm().receive (id_status.source(), ids,  id_tag);
#else
              std::vector<dof_id_type> & ids (recv_ids[0]);
              ids = xfer_ids[blk];
#endif

              // note its possible we didn't receive values for objects in
              // this block if they have no components allocated.
              for (std::size_t idx=0; idx<ids.size(); idx+=2)
                {
                  const dof_id_type
                    local_idx          = ids[idx+0]-first_object,
                    n_vals_tot_allvecs = ids[idx+1];

                  libmesh_assert_less (local_idx, n_objects_blk);
                  libmesh_assert_less (local_idx, obj_val_offsets.size());

                  obj_val_offsets[local_idx] = n_vals_tot_allvecs;
                }

#ifdef LIBMESH_HAVE_MPI
              // blocking receive values for this block, imposing no particular order on processor
              Parallel::Status val_status  (this->comm().probe (Parallel::any_source, val_tag));
              std::vector<Number> & vals    (recv_vals[val_status.source()]);
              this->comm().receive (val_status.source(), vals, val_tag);
#else
              // straight copy without MPI
              std::vector<Number> & vals (recv_vals[0]);
              vals = send_vals[blk];
#endif

              n_val_recvd_blk += vals.size();
            }

          // We need the offests into the output_vals vector for each object.
          // fortunately, this is simply the partial sum of the total number
          // of components for each object
          std::partial_sum(obj_val_offsets.begin(), obj_val_offsets.end(),
                           obj_val_offsets.begin());

          // wait on any previous asynchronous IO - this *must* complete before
          // we start messing with the output_vals buffer!
          if (async_io.get()) async_io->join();

          // this is the actual output buffer that will be written to disk.
          // at ths point we finally know wha size it will be.
          output_vals.resize(n_val_recvd_blk);

          // pack data from all processors into output values
          for (unsigned int proc=0; proc<this->n_processors(); proc++)
            {
              const std::vector<dof_id_type> & ids (recv_ids [proc]);
              const std::vector<Number>      & vals(recv_vals[proc]);
              std::vector<Number>::const_iterator proc_vals(vals.begin());

              for (std::size_t idx=0; idx<ids.size(); idx+=2)
                {
                  const dof_id_type
                    local_idx          = ids[idx+0]-first_object,
                    n_vals_tot_allvecs = ids[idx+1];

                  // put this object's data into the proper location
                  // in  the output buffer
                  std::vector<Number>::iterator out_vals(output_vals.begin());
                  if (local_idx != 0)
                    std::advance (out_vals, obj_val_offsets[local_idx-1]);

                  for (unsigned int val=0; val<n_vals_tot_allvecs; val++, ++out_vals, ++proc_vals)
                    {
                      libmesh_assert (out_vals  != output_vals.end());
                      libmesh_assert (proc_vals != vals.end());
                      *out_vals = *proc_vals;
                    }
                }
            }

          // output_vals buffer is now filled for this block.
          // write it to disk
          async_io.reset(new Threads::Thread(threaded_io));
          written_length += output_vals.size();
        }

      // wait on any previous asynchronous IO - this *must* complete before
      // our stuff goes out of scope
      async_io->join();
    }

  Parallel::wait(id_requests);
  Parallel::wait(val_requests);

  // we need some synchronization here.  Because this method
  // can be called for a range of nodes, then a range of elements,
  // we need some mechanism to prevent processors from racing past
  // to the next range and overtaking ongoing communication.  one
  // approach would be to figure out unique tags for each range,
  // but for now we just impose a barrier here.  And might as
  // well have it do some useful work.
  this->comm().broadcast(written_length);

  return written_length;
}



unsigned int System::write_SCALAR_dofs (const NumericVector<Number> & vec,
                                        const unsigned int var,
                                        Xdr & io) const
{
  unsigned int written_length=0;
  std::vector<Number> vals; // The raw values for the local objects in the current block
  // Collect the SCALARs for the current variable
  if (this->processor_id() == (this->n_processors()-1))
    {
      const DofMap & dof_map = this->get_dof_map();
      std::vector<dof_id_type> SCALAR_dofs;
      dof_map.SCALAR_dof_indices(SCALAR_dofs, var);
      const unsigned int n_scalar_dofs = cast_int<unsigned int>
        (SCALAR_dofs.size());

      for(unsigned int i=0; i<n_scalar_dofs; i++)
        {
          vals.push_back( vec(SCALAR_dofs[i]) );
        }
    }

#ifdef LIBMESH_HAVE_MPI
  if ( this->n_processors() > 1 )
    {
      const Parallel::MessageTag val_tag(1);

      // Post the receive on processor 0
      if ( this->processor_id() == 0 )
        {
          this->comm().receive(this->n_processors()-1, vals, val_tag);
        }

      // Send the data to processor 0
      if (this->processor_id() == (this->n_processors()-1))
        {
          this->comm().send(0, vals, val_tag);
        }
    }
#endif

  // -------------------------------------------------------
  // Write the output on processor 0.
  if (this->processor_id() == 0)
    {
      const unsigned int vals_size =
        cast_int<unsigned int>(vals.size());
      io.data_stream (&vals[0], vals_size);
      written_length += vals_size;
    }

  return written_length;
}



dof_id_type System::write_serialized_vector (Xdr & io,
                                             const NumericVector<Number> & vec) const
{
  parallel_object_only();

  libmesh_assert (io.writing());

  dof_id_type vec_length = vec.size();
  if (this->processor_id() == 0) io.data (vec_length, "# vector length");

  dof_id_type written_length = 0;

  //---------------------------------
  // Collect the values for all nodes
  written_length += cast_int<dof_id_type>
    (this->write_serialized_blocked_dof_objects (std::vector<const NumericVector<Number> *>(1,&vec),
                                                 this->get_mesh().n_nodes(),
                                                 this->get_mesh().local_nodes_begin(),
                                                 this->get_mesh().local_nodes_end(),
                                                 io));

  //------------------------------------
  // Collect the values for all elements
  written_length += cast_int<dof_id_type>
    (this->write_serialized_blocked_dof_objects (std::vector<const NumericVector<Number> *>(1,&vec),
                                                 this->get_mesh().n_elem(),
                                                 this->get_mesh().local_elements_begin(),
                                                 this->get_mesh().local_elements_end(),
                                                 io));

  //-------------------------------------------
  // Finally loop over all the SCALAR variables
  for (unsigned int var=0; var<this->n_vars(); var++)
    if(this->variable(var).type().family == SCALAR)
      {
        written_length +=
          this->write_SCALAR_dofs (vec, var, io);
      }

  if (this->processor_id() == 0)
    libmesh_assert_equal_to (written_length, vec_length);

  return written_length;
}


template <typename InValType>
std::size_t System::read_serialized_vectors (Xdr & io,
                                             const std::vector<NumericVector<Number> *> & vectors) const
{
  parallel_object_only();

  // Error checking
  // #ifndef NDEBUG
  //   // In parallel we better be reading a parallel vector -- if not
  //   // we will not set all of its components below!!
  //   if (this->n_processors() > 1)
  //     {
  //       libmesh_assert (vec.type() == PARALLEL ||
  //       vec.type() == GHOSTED);
  //     }
  // #endif

  libmesh_assert (io.reading());

  if (this->processor_id() == 0)
    {
      // sizes
      unsigned int num_vecs=0;
      dof_id_type vector_length=0;

      // Get the number of vectors
      io.data(num_vecs);
      // Get the buffer size
      io.data(vector_length);

      libmesh_assert_equal_to (num_vecs, vectors.size());

      if (num_vecs != 0)
        {
          libmesh_assert_not_equal_to (vectors[0], 0);
          libmesh_assert_equal_to     (vectors[0]->size(), vector_length);
        }
    }

  // no need to actually communicate these.
  // this->comm().broadcast(num_vecs);
  // this->comm().broadcast(vector_length);

  // Cache these - they are not free!
  const dof_id_type
    n_nodes = this->get_mesh().n_nodes(),
    n_elem  = this->get_mesh().n_elem();

  std::size_t read_length = 0;

  //---------------------------------
  // Collect the values for all nodes
  read_length +=
    this->read_serialized_blocked_dof_objects (n_nodes,
                                               this->get_mesh().local_nodes_begin(),
                                               this->get_mesh().local_nodes_end(),
                                               InValType(),
                                               io,
                                               vectors);

  //------------------------------------
  // Collect the values for all elements
  read_length +=
    this->read_serialized_blocked_dof_objects (n_elem,
                                               this->get_mesh().local_elements_begin(),
                                               this->get_mesh().local_elements_end(),
                                               InValType(),
                                               io,
                                               vectors);

  //-------------------------------------------
  // Finally loop over all the SCALAR variables
  for (std::size_t vec=0; vec<vectors.size(); vec++)
    for (unsigned int var=0; var<this->n_vars(); var++)
      if(this->variable(var).type().family == SCALAR)
        {
          libmesh_assert_not_equal_to (vectors[vec], 0);

          read_length +=
            this->read_SCALAR_dofs (var, io, vectors[vec]);
        }

  //---------------------------------------
  // last step - must close all the vectors
  for (std::size_t vec=0; vec<vectors.size(); vec++)
    {
      libmesh_assert_not_equal_to (vectors[vec], 0);
      vectors[vec]->close();
    }

  return read_length;
}



std::size_t System::write_serialized_vectors (Xdr & io,
                                              const std::vector<const NumericVector<Number> *> & vectors) const
{
  parallel_object_only();

  libmesh_assert (io.writing());

  // Cache these - they are not free!
  const dof_id_type
    n_nodes       = this->get_mesh().n_nodes(),
    n_elem        = this->get_mesh().n_elem();

  std::size_t written_length = 0;

  if (this->processor_id() == 0)
    {
      unsigned int
        n_vec    = cast_int<unsigned int>(vectors.size());
      dof_id_type
        vec_size = vectors.empty() ? 0 : vectors[0]->size();
      // Set the number of vectors
      io.data(n_vec, "# number of vectors");
      // Set the buffer size
      io.data(vec_size, "# vector length");
    }

  //---------------------------------
  // Collect the values for all nodes
  written_length +=
    this->write_serialized_blocked_dof_objects (vectors,
                                                n_nodes,
                                                this->get_mesh().local_nodes_begin(),
                                                this->get_mesh().local_nodes_end(),
                                                io);

  //------------------------------------
  // Collect the values for all elements
  written_length +=
    this->write_serialized_blocked_dof_objects (vectors,
                                                n_elem,
                                                this->get_mesh().local_elements_begin(),
                                                this->get_mesh().local_elements_end(),
                                                io);

  //-------------------------------------------
  // Finally loop over all the SCALAR variables
  for (std::size_t vec=0; vec<vectors.size(); vec++)
    for (unsigned int var=0; var<this->n_vars(); var++)
      if(this->variable(var).type().family == SCALAR)
        {
          libmesh_assert_not_equal_to (vectors[vec], 0);

          written_length +=
            this->write_SCALAR_dofs (*vectors[vec], var, io);
        }

  return written_length;
}




template void System::read_parallel_data<Number> (Xdr & io, const bool read_additional_data);
template void System::read_serialized_data<Number> (Xdr & io, const bool read_additional_data);
template numeric_index_type System::read_serialized_vector<Number> (Xdr & io, NumericVector<Number> * vec);
template std::size_t System::read_serialized_vectors<Number> (Xdr & io, const std::vector<NumericVector<Number> *> & vectors) const;
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template void System::read_parallel_data<Real> (Xdr & io, const bool read_additional_data);
template void System::read_serialized_data<Real> (Xdr & io, const bool read_additional_data);
template numeric_index_type System::read_serialized_vector<Real> (Xdr & io, NumericVector<Number> * vec);
template std::size_t System::read_serialized_vectors<Real> (Xdr & io, const std::vector<NumericVector<Number> *> & vectors) const;
#endif

} // namespace libMesh
