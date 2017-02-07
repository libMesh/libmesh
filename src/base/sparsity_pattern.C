// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public  License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// Local Includes
#include "libmesh/coupling_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/ghosting_functor.h"
#include "libmesh/sparsity_pattern.h"


namespace libMesh
{
namespace SparsityPattern
{

//-------------------------------------------------------
// we need to implement these constructors here so that
// a full DofMap definition is available.
Build::Build (const MeshBase & mesh_in,
              const DofMap & dof_map_in,
              const CouplingMatrix * dof_coupling_in,
              std::set<GhostingFunctor *> coupling_functors_in,
              const bool implicit_neighbor_dofs_in,
              const bool need_full_sparsity_pattern_in) :
  ParallelObject(dof_map_in),
  mesh(mesh_in),
  dof_map(dof_map_in),
  dof_coupling(dof_coupling_in),
  coupling_functors(coupling_functors_in),
  implicit_neighbor_dofs(implicit_neighbor_dofs_in),
  need_full_sparsity_pattern(need_full_sparsity_pattern_in),
  sparsity_pattern(),
  nonlocal_pattern(),
  n_nz(),
  n_oz()
{}



Build::Build (Build & other, Threads::split) :
  ParallelObject(other),
  mesh(other.mesh),
  dof_map(other.dof_map),
  dof_coupling(other.dof_coupling),
  coupling_functors(other.coupling_functors),
  implicit_neighbor_dofs(other.implicit_neighbor_dofs),
  need_full_sparsity_pattern(other.need_full_sparsity_pattern),
  sparsity_pattern(),
  nonlocal_pattern(),
  n_nz(),
  n_oz()
{}



#if defined(__GNUC__) && (__GNUC__ < 4) && !defined(__INTEL_COMPILER)

void _dummy_function(void) {}

#endif



void Build::handle_vi_vj
  (const Elem * partner,
   const std::vector<dof_id_type> & element_dofs_i,
   unsigned int vj)
{
  const unsigned int n_dofs_on_element_i =
    cast_int<unsigned int>(element_dofs_i.size());

  const processor_id_type proc_id           = mesh.processor_id();
  const dof_id_type first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type end_dof_on_proc   = dof_map.end_dof(proc_id);

  std::vector<dof_id_type> element_dofs_j,
                           dofs_to_add;

  // Find element dofs for variable vj
  dof_map.dof_indices (partner, element_dofs_j, vj);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  dof_map.find_connected_dofs (element_dofs_j);
#endif

  // We can be more efficient if we sort the element DOFs
  // into increasing order
  std::sort (element_dofs_j.begin(), element_dofs_j.end());

  const unsigned int n_dofs_on_element_j =
    cast_int<unsigned int>(element_dofs_j.size());

  // there might be 0 dofs for the other variable on the same element
  // (when subdomain variables do not overlap) and that's when we do
  // not do anything
  if (n_dofs_on_element_j > 0)
    {
      for (unsigned int i=0; i<n_dofs_on_element_i; i++)
        {
          const dof_id_type ig = element_dofs_i[i];

          SparsityPattern::Row * row;

          // We save non-local row components for now so we can
          // communicate them to other processors later.

          if ((ig >= first_dof_on_proc) &&
              (ig <  end_dof_on_proc))
            {
              // This is what I mean
              // libmesh_assert_greater_equal ((ig - first_dof_on_proc), 0);
              // but do the test like this because ig and
              // first_dof_on_proc are unsigned ints
              libmesh_assert_greater_equal (ig, first_dof_on_proc);
              libmesh_assert_less (ig, (sparsity_pattern.size() +
                                        first_dof_on_proc));

              row = &sparsity_pattern[ig - first_dof_on_proc];
            }
          else
            {
              row = &nonlocal_pattern[ig];
            }

          // If the row is empty we will add *all*
          // the element j DOFs, so just do that.
          if (row->empty())
            {
              row->insert(row->end(),
                          element_dofs_j.begin(),
                          element_dofs_j.end());
            }
          else
            {
              // Build a list of the DOF indices not found in the
              // sparsity pattern
              dofs_to_add.clear();

              // Cache iterators.  Low will move forward, subsequent
              // searches will be on smaller ranges
              SparsityPattern::Row::iterator
                low  = std::lower_bound
                (row->begin(), row->end(), element_dofs_j.front()),
                high = std::upper_bound
                (low,          row->end(), element_dofs_j.back());

              for (unsigned int j=0; j<n_dofs_on_element_j; j++)
                {
                  const dof_id_type jg = element_dofs_j[j];

                  // See if jg is in the sorted range
                  std::pair<SparsityPattern::Row::iterator,
                            SparsityPattern::Row::iterator>
                    pos = std::equal_range (low, high, jg);

                  // Must add jg if it wasn't found
                  if (pos.first == pos.second)
                    dofs_to_add.push_back(jg);

                  // pos.first is now a valid lower bound for any
                  // remaining element j DOFs. (That's why we sorted them.)
                  // Use it for the next search
                  low = pos.first;
                }

              // Add to the sparsity pattern
              if (!dofs_to_add.empty())
                {
                  const std::size_t old_size = row->size();

                  row->insert (row->end(),
                               dofs_to_add.begin(),
                               dofs_to_add.end());

                  SparsityPattern::sort_row
                    (row->begin(), row->begin()+old_size,
                     row->end());
                }
            }
        } // End dofs-of-var-i loop
    } // End if-dofs-of-var-j
}



void Build::operator()(const ConstElemRange & range)
{
  // Compute the sparsity structure of the global matrix.  This can be
  // fed into a PetscMatrix to allocate exacly the number of nonzeros
  // necessary to store the matrix.  This algorithm should be linear
  // in the (# of elements)*(# nodes per element)
  const processor_id_type proc_id           = mesh.processor_id();
  const dof_id_type n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const dof_id_type first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type end_dof_on_proc   = dof_map.end_dof(proc_id);

  sparsity_pattern.resize(n_dofs_on_proc);

  // Handle dof coupling specified by library and user coupling functors
    {
      const unsigned int n_var = dof_map.n_variables();

      std::vector<dof_id_type> element_dofs_i;

      std::vector<const Elem *> coupled_neighbors;
      for (ConstElemRange::const_iterator elem_it = range.begin() ; elem_it != range.end(); ++elem_it)
        {
          const Elem * const elem = *elem_it;

          // Make some fakey element iterators defining a range
          // pointing to only this element.
          Elem * const * elempp = const_cast<Elem * const *>(&elem);
          Elem * const * elemend = elempp+1;
          const MeshBase::const_element_iterator fake_elem_it =
            MeshBase::const_element_iterator(elempp, elemend, Predicates::NotNull<Elem * const *>());
          const MeshBase::const_element_iterator fake_elem_end =
            MeshBase::const_element_iterator(elemend, elemend, Predicates::NotNull<Elem * const *>());

          GhostingFunctor::map_type elements_to_couple;

          // Man, I wish we had guaranteed unique_ptr availability...
          std::set<CouplingMatrix*> temporary_coupling_matrices;

          dof_map.merge_ghost_functor_outputs
            (elements_to_couple,
             temporary_coupling_matrices,
             dof_map.coupling_functors_begin(),
             dof_map.coupling_functors_end(),
             fake_elem_it, fake_elem_end, DofObject::invalid_processor_id);

          for (unsigned int vi=0; vi<n_var; vi++)
            {
              // Find element dofs for variable vi
              dof_map.dof_indices (elem, element_dofs_i, vi);
#ifdef LIBMESH_ENABLE_CONSTRAINTS
              dof_map.find_connected_dofs (element_dofs_i);
#endif

              // We can be more efficient if we sort the element DOFs
              // into increasing order
              std::sort(element_dofs_i.begin(), element_dofs_i.end());

              GhostingFunctor::map_type::iterator        etg_it = elements_to_couple.begin();
              const GhostingFunctor::map_type::iterator etg_end = elements_to_couple.end();
              for (; etg_it != etg_end; ++etg_it)
                {
                  const Elem * const partner = etg_it->first;
                  const CouplingMatrix *ghost_coupling = etg_it->second;

                  // Loop over coupling matrix row variables if we have a
                  // coupling matrix, or all variables if not.
                  if (ghost_coupling)
                    {
                      libmesh_assert_equal_to (ghost_coupling->size(), n_var);
                      ConstCouplingRow ccr(vi, *ghost_coupling);

                      for (ConstCouplingRow::const_iterator  it = ccr.begin(),
                                                            end = ccr.end();
                           it != end; ++it)
                        this->handle_vi_vj(partner, element_dofs_i, *it);
                    }
                  else
                    {
                      for (unsigned int vj = 0; vj != n_var; ++vj)
                        this->handle_vi_vj(partner, element_dofs_i, vj);
                    }
                } // End ghosted element loop
            } // End vi loop

          for (std::set<CouplingMatrix*>::iterator
                 it  = temporary_coupling_matrices.begin(),
                 end = temporary_coupling_matrices.begin();
               it != end; ++it)
            delete *it;

        } // End range element loop
    } // End ghosting functor section

  // Now a new chunk of sparsity structure is built for all of the
  // DOFs connected to our rows of the matrix.

  // If we're building a full sparsity pattern, then we've got
  // complete rows to work with, so we can just count them from
  // scratch.
  if (need_full_sparsity_pattern)
    {
      n_nz.clear();
      n_oz.clear();
    }

  n_nz.resize (n_dofs_on_proc, 0);
  n_oz.resize (n_dofs_on_proc, 0);

  for (dof_id_type i=0; i<n_dofs_on_proc; i++)
    {
      // Get the row of the sparsity pattern
      SparsityPattern::Row & row = sparsity_pattern[i];

      for (std::size_t j=0; j<row.size(); j++)
        if ((row[j] < first_dof_on_proc) || (row[j] >= end_dof_on_proc))
          n_oz[i]++;
        else
          n_nz[i]++;

      // If we're not building a full sparsity pattern, then we want
      // to avoid overcounting these entries as much as possible.
      if (!need_full_sparsity_pattern)
        row.clear();
    }
}



void Build::join (const SparsityPattern::Build & other)
{
  const processor_id_type proc_id           = mesh.processor_id();
  const dof_id_type       n_global_dofs     = dof_map.n_dofs();
  const dof_id_type       n_dofs_on_proc    = dof_map.n_dofs_on_processor(proc_id);
  const dof_id_type       first_dof_on_proc = dof_map.first_dof(proc_id);
  const dof_id_type       end_dof_on_proc   = dof_map.end_dof(proc_id);

  libmesh_assert_equal_to (sparsity_pattern.size(), other.sparsity_pattern.size());
  libmesh_assert_equal_to (n_nz.size(), sparsity_pattern.size());
  libmesh_assert_equal_to (n_oz.size(), sparsity_pattern.size());

  for (dof_id_type r=0; r<n_dofs_on_proc; r++)
    {
      // increment the number of on and off-processor nonzeros in this row
      // (note this will be an upper bound unless we need the full sparsity pattern)
      if (need_full_sparsity_pattern)
        {
          SparsityPattern::Row       & my_row    = sparsity_pattern[r];
          const SparsityPattern::Row & their_row = other.sparsity_pattern[r];

          // simple copy if I have no dofs
          if (my_row.empty())
            my_row = their_row;

          // otherwise add their DOFs to mine, resort, and re-unique the row
          else if (!their_row.empty()) // do nothing for the trivial case where
            {                          // their row is empty
              my_row.insert (my_row.end(),
                             their_row.begin(),
                             their_row.end());

              // We cannot use SparsityPattern::sort_row() here because it expects
              // the [begin,middle) [middle,end) to be non-overlapping.  This is not
              // necessarily the case here, so use std::sort()
              std::sort (my_row.begin(), my_row.end());

              my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
            }

          // fix the number of on and off-processor nonzeros in this row
          n_nz[r] = n_oz[r] = 0;

          for (std::size_t j=0; j<my_row.size(); j++)
            if ((my_row[j] < first_dof_on_proc) || (my_row[j] >= end_dof_on_proc))
              n_oz[r]++;
            else
              n_nz[r]++;
        }
      else
        {
          n_nz[r] += other.n_nz[r];
          n_nz[r] = std::min(n_nz[r], n_dofs_on_proc);
          n_oz[r] += other.n_oz[r];
          n_oz[r] =std::min(n_oz[r], static_cast<dof_id_type>(n_global_dofs-n_nz[r]));
        }
    }

  // Move nonlocal row information to ourselves; the other thread
  // won't need it in the map after that.
  NonlocalGraph::const_iterator it = other.nonlocal_pattern.begin();
  for (; it != other.nonlocal_pattern.end(); ++it)
    {
#ifndef NDEBUG
      const dof_id_type dof_id = it->first;

      processor_id_type dbg_proc_id = 0;
      while (dof_id >= dof_map.end_dof(dbg_proc_id))
        dbg_proc_id++;
      libmesh_assert (dbg_proc_id != this->processor_id());
#endif

      const SparsityPattern::Row & their_row = it->second;

      // We should have no empty values in a map
      libmesh_assert (!their_row.empty());

      NonlocalGraph::iterator my_it = nonlocal_pattern.find(it->first);
      if (my_it == nonlocal_pattern.end())
        {
          //          nonlocal_pattern[it->first].swap(their_row);
          nonlocal_pattern[it->first] = their_row;
        }
      else
        {
          SparsityPattern::Row & my_row = my_it->second;

          my_row.insert (my_row.end(),
                         their_row.begin(),
                         their_row.end());

          // We cannot use SparsityPattern::sort_row() here because it expects
          // the [begin,middle) [middle,end) to be non-overlapping.  This is not
          // necessarily the case here, so use std::sort()
          std::sort (my_row.begin(), my_row.end());

          my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
        }
    }
}



void Build::parallel_sync ()
{
  parallel_object_only();
  libmesh_assert(this->comm().verify(need_full_sparsity_pattern));

  const dof_id_type n_global_dofs   = dof_map.n_dofs();
  const dof_id_type n_dofs_on_proc  = dof_map.n_dofs_on_processor(this->processor_id());
  const dof_id_type local_first_dof = dof_map.first_dof();
  const dof_id_type local_end_dof   = dof_map.end_dof();

  // Trade sparsity rows with other processors
  for (processor_id_type p=1; p != this->n_processors(); ++p)
    {
      // Push to processor procup while receiving from procdown
      processor_id_type procup =
        cast_int<processor_id_type>((this->processor_id() + p) %
                                    this->n_processors());
      processor_id_type procdown =
        cast_int<processor_id_type>((this->n_processors() + this->processor_id() - p) %
                                    this->n_processors());

      // Pack the sparsity pattern rows to push to procup
      std::vector<dof_id_type> pushed_row_ids,
        pushed_row_ids_to_me;
      std::vector<std::vector<dof_id_type> > pushed_rows,
        pushed_rows_to_me;

      // Move nonlocal row information to a structure to send it from;
      // we don't need it in the map after that.
      NonlocalGraph::iterator it = nonlocal_pattern.begin();
      while (it != nonlocal_pattern.end())
        {
          const dof_id_type dof_id = it->first;
          processor_id_type proc_id = 0;
          while (dof_id >= dof_map.end_dof(proc_id))
            proc_id++;

          libmesh_assert (proc_id != this->processor_id());

          if (proc_id == procup)
            {
              pushed_row_ids.push_back(dof_id);

              // We can't just do the swap trick here, thanks to the
              // differing vector allocators?
              pushed_rows.push_back(std::vector<dof_id_type>());
              pushed_rows.back().assign
                (it->second.begin(), it->second.end());

              nonlocal_pattern.erase(it++);
            }
          else
            ++it;
        }

      this->comm().send_receive(procup, pushed_row_ids,
                                procdown, pushed_row_ids_to_me);
      this->comm().send_receive(procup, pushed_rows,
                                procdown, pushed_rows_to_me);
      pushed_row_ids.clear();
      pushed_rows.clear();

      const std::size_t n_rows = pushed_row_ids_to_me.size();
      for (std::size_t i=0; i != n_rows; ++i)
        {
          const dof_id_type r = pushed_row_ids_to_me[i];
          const dof_id_type my_r = r - local_first_dof;

          std::vector<dof_id_type> & their_row = pushed_rows_to_me[i];

          if (need_full_sparsity_pattern)
            {
              SparsityPattern::Row & my_row =
                sparsity_pattern[my_r];

              // They wouldn't have sent an empty row
              libmesh_assert(!their_row.empty());

              // We can end up with an empty row on a dof that touches our
              // inactive elements but not our active ones
              if (my_row.empty())
                {
                  my_row.assign (their_row.begin(),
                                 their_row.end());
                }
              else
                {
                  my_row.insert (my_row.end(),
                                 their_row.begin(),
                                 their_row.end());

                  // We cannot use SparsityPattern::sort_row() here because it expects
                  // the [begin,middle) [middle,end) to be non-overlapping.  This is not
                  // necessarily the case here, so use std::sort()
                  std::sort (my_row.begin(), my_row.end());

                  my_row.erase(std::unique (my_row.begin(), my_row.end()), my_row.end());
                }

              // fix the number of on and off-processor nonzeros in this row
              n_nz[my_r] = n_oz[my_r] = 0;

              for (std::size_t j=0; j<my_row.size(); j++)
                if ((my_row[j] < local_first_dof) || (my_row[j] >= local_end_dof))
                  n_oz[my_r]++;
                else
                  n_nz[my_r]++;
            }
          else
            {
              for (std::size_t j=0; j<their_row.size(); j++)
                if ((their_row[j] < local_first_dof) || (their_row[j] >= local_end_dof))
                  n_oz[my_r]++;
                else
                  n_nz[my_r]++;

              n_nz[my_r] = std::min(n_nz[my_r], n_dofs_on_proc);
              n_oz[my_r] = std::min(n_oz[my_r],
                                    static_cast<dof_id_type>(n_global_dofs-n_nz[my_r]));
            }
        }
    }

  // We should have sent everything at this point.
  libmesh_assert (nonlocal_pattern.empty());
}


} // namespace SparsityPattern
} // namespace libMesh
