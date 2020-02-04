// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local includes
#include "libmesh/mesh_iterators_impl.h"

namespace libMesh
{

// Instantiate various element iterator accessor functions.
INSTANTIATE_ELEM_ACCESSORS(Real, elements,                        NotNull,              EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, active_elements,                 Active,               EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, not_active_elements,             NotActive,            EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, ancestor_elements,               Ancestor,             EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, not_ancestor_elements,           NotAncestor,          EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, subactive_elements,              SubActive,            EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, not_subactive_elements,          NotSubActive,         EMPTY,                          EMPTY)
INSTANTIATE_ELEM_ACCESSORS(Real, local_elements,                  Local,                EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, semilocal_elements,              ActiveSemiLocal,      EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, active_semilocal_elements,       ActiveSemiLocal,      EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, facelocal_elements,              FaceLocal,            EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, not_local_elements,              NotLocal,             EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, active_local_elements,           ActiveLocal,          EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, active_not_local_elements,       ActiveNotLocal,       EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, level_elements,                  Level,                unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(Real, not_level_elements,              NotLevel,             unsigned int level,             level)
INSTANTIATE_ELEM_ACCESSORS(Real, pid_elements,                    PID,                  processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(Real, type_elements,                   Type,                 ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(Real, active_type_elements,            ActiveType,           ElemType type,                  type)
INSTANTIATE_ELEM_ACCESSORS(Real, active_pid_elements,             ActivePID,            processor_id_type proc_id,      proc_id)
INSTANTIATE_ELEM_ACCESSORS(Real, active_subdomain_elements,       ActiveSubdomain,      subdomain_id_type subdomain_id, subdomain_id)
INSTANTIATE_ELEM_ACCESSORS(Real, active_subdomain_set_elements,   ActiveSubdomainSet,   std::set<subdomain_id_type> ss, ss)
INSTANTIATE_ELEM_ACCESSORS(Real, ghost_elements,                  Ghost,                EMPTY,                          this->processor_id())
INSTANTIATE_ELEM_ACCESSORS(Real, evaluable_elements,              Evaluable,            const DofMap & dof_map LIBMESH_COMMA unsigned int var_num, dof_map, var_num)
INSTANTIATE_ELEM_ACCESSORS(Real, unpartitioned_elements,          PID,                  EMPTY,                          DofObject::invalid_processor_id)
INSTANTIATE_ELEM_ACCESSORS(Real, active_unpartitioned_elements,   ActivePID,            EMPTY,                          DofObject::invalid_processor_id)

#ifdef LIBMESH_ENABLE_AMR
INSTANTIATE_ELEM_ACCESSORS(Real, flagged_elements,                Flagged,              unsigned char rflag,            rflag)
INSTANTIATE_ELEM_ACCESSORS(Real, flagged_pid_elements,            FlaggedPID,           unsigned char rflag LIBMESH_COMMA processor_id_type pid,    rflag, pid)
#endif

INSTANTIATE_ELEM_ACCESSORS(Real, local_level_elements,            LocalLevel,           unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(Real, local_not_level_elements,        LocalNotLevel,        unsigned int level,             this->processor_id(), level)
INSTANTIATE_ELEM_ACCESSORS(Real, active_local_subdomain_elements, ActiveLocalSubdomain, subdomain_id_type subdomain_id, this->processor_id(), subdomain_id)

// Instantiate various node iterator accessor functions.
INSTANTIATE_NODE_ACCESSORS(Real, nodes,        NotNull, EMPTY,                               EMPTY)
INSTANTIATE_NODE_ACCESSORS(Real, active_nodes, Active,  EMPTY,                               EMPTY)
INSTANTIATE_NODE_ACCESSORS(Real, local_nodes,  Local,   EMPTY,                               this->processor_id())
INSTANTIATE_NODE_ACCESSORS(Real, pid_nodes,    PID,     processor_id_type proc_id,           proc_id)
INSTANTIATE_NODE_ACCESSORS(Real, bnd_nodes,    BND,     EMPTY,                               this->get_boundary_info())
INSTANTIATE_NODE_ACCESSORS(Real, bid_nodes,    BID,     boundary_id_type bndry_id, bndry_id, this->get_boundary_info())
INSTANTIATE_NODE_ACCESSORS(Real, evaluable_nodes, Evaluable, const DofMap & dof_map LIBMESH_COMMA unsigned int var_num, dof_map, var_num)

} // namespace libMesh
