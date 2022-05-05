// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/file_history_data.h"

#include "libmesh/enum_xdr_mode.h"
#include "libmesh/equation_systems.h"

namespace libMesh
{
  void FileHistoryData::store_initial_solution()
  {
    // The initial data should only be stored once.
    libmesh_assert(previously_stored == false);

    time_stamp = 0;

    primal_filename = "primal.out.xda.";

    primal_filename += std::to_string(time_stamp);

    _system.get_equation_systems().write (primal_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);

    // We wont know the deltat taken at this timestep until the solve is completed, which is done after the store operation.
    deltat_at = std::numeric_limits<double>::signaling_NaN();

    previously_stored = true;
  }

  void FileHistoryData::store_primal_solution(stored_data_iterator stored_datum)
  {
    // An iterator to the datum being stored has been passed. This is taken to imply
    // that we are in a sequential time stepping sequence.

    stored_data_iterator stored_datum_last = stored_datum;
    stored_datum_last--;

    time_stamp = (stored_datum_last->second)->get_time_stamp() + 1;

    primal_filename = "primal.out.xda.";

    primal_filename += std::to_string(time_stamp);

    _system.get_equation_systems().write (primal_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);

    // We dont know the deltat that will be taken at the current timestep.
    deltat_at = std::numeric_limits<double>::signaling_NaN();

    // But we know what deltat got us to the current time.
    (stored_datum_last->second)->set_deltat_at(_system.time_solver->TimeSolver::last_completed_timestep_size());

    previously_stored = true;
  }

  void FileHistoryData::store_adjoint_solution()
    {
     adjoint_filename = "adjoint.out.xda.";

     adjoint_filename += std::to_string(time_stamp);

     _system.get_equation_systems().write(adjoint_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);
    }

  void FileHistoryData::rewrite_stored_solution()
    {
     // We are rewriting.
     libmesh_assert(previously_stored == true);

     _system.get_equation_systems().write(primal_filename, WRITE, EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA);
    }

   void FileHistoryData::retrieve_primal_solution()
    {
     // Read in the primal solution stored at the current recovery time from the disk
     _system.get_equation_systems().read(primal_filename, READ, EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA);
    }

   void FileHistoryData::retrieve_adjoint_solution()
    {
     // Read in the adjoint solution stored at the current recovery time from the disk
     _system.get_equation_systems().read(adjoint_filename, READ, EquationSystems::READ_DATA | EquationSystems::READ_ADDITIONAL_DATA);
    }
}
