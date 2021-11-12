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

#include "libmesh/memory_history_data.h"

namespace libMesh
{
    void MemoryHistoryData::store_initial_solution()
    {
     // The initial data should only be stored once.
     libmesh_assert(previously_stored == false);

     time_stamp = 0;

     deltat_at = std::numeric_limits<double>::signaling_NaN();

     store_vectors();

     previously_stored = true;
    }

    void MemoryHistoryData::store_primal_solution(stored_data_iterator stored_datum)
    {
     stored_data_iterator stored_datum_last = stored_datum;
     stored_datum_last--;

     time_stamp = (stored_datum_last->second)->get_time_stamp() + 1;

     // For the current time instant, we dont know yet what timestep the solver might decide, so a placeholder NaN for now.
     deltat_at = std::numeric_limits<double>::signaling_NaN();

     (stored_datum_last->second)->set_deltat_at(_system.time_solver->TimeSolver::last_completed_timestep_size());

     store_vectors();

     previously_stored = true;
    }

    void MemoryHistoryData::store_adjoint_solution()
    {
     libmesh_error_msg("For MemorySolutionHistory, primal and adjoints are stored in the same container.");
    }

    void MemoryHistoryData::rewrite_stored_solution()
    {
      // We are rewriting.
      libmesh_assert(previously_stored == true);

      store_vectors();
    }

    void MemoryHistoryData::retrieve_primal_solution()
    {
     retrieve_vectors();
    }

    void MemoryHistoryData::retrieve_adjoint_solution()
    {
     retrieve_vectors();
    }

    void MemoryHistoryData::store_vectors()
    {
     // Now save all the preserved vectors in stored_datum
     // Loop over all the system vectors
     for (System::vectors_iterator vec = _system.vectors_begin(),
          vec_end = _system.vectors_end(); vec != vec_end; ++vec)
     {
      // The name of this vector
      const std::string & vec_name = vec->first;

      // Store the vector if it is to be preserved
      if (_system.vector_preservation(vec_name))
        {
         stored_vecs[vec_name] = vec->second->clone();
        }
     }

     // Of course, we will usually save the actual solution
     std::string _solution("_solution");
     if (_system.project_solution_on_reinit())
     {
      stored_vecs[_solution] = _system.solution->clone();
     }
    }

    void MemoryHistoryData::retrieve_vectors()
    {
      // We are reading, hopefully something has been written before
      libmesh_assert(previously_stored == true);

      map_type::iterator vec = stored_vecs.begin();
      map_type::iterator vec_end = stored_vecs.end();

      // Loop over all the saved vectors
      for (; vec != vec_end; ++vec)
      {
       // The name of this vector
       const std::string & vec_name = vec->first;

       // Get the vec_name entry in the saved vectors map and set the
       // current system vec[vec_name] entry to it
       if (vec_name != "_solution")
         _system.get_vector(vec_name) = *(vec->second);
      }

      std::string _solution("_solution");
      *(_system.solution) = *(stored_vecs[_solution]);

    }
}