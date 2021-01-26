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



#ifndef LIBMESH_INTER_MESH_PROJECTION_H
#define LIBMESH_INTER_MESH_PROJECTION_H

// Local includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh_function.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/fem_function_base.h"

// C++ includes
#include <cstddef>
#include <vector>

namespace libMesh
{

    // Forward declarations

    /**
     * This class implements inter mesh projection, i.e. projection of
     * vectors defined on a given mesh (from_mesh associated with from_system)
     * to another mesh (to_mesh of to_system).
     */

    class InterMeshProjection
    {
        public:

        // Constructor, specifies the _from_system whose vectors will be
        // projected onto the _to_mesh
        InterMeshProjection(System & _from_system, System & _to_mesh);

        // Projects from_system vectors onto the to_mesh
        void project_system_vectors();

        static Number fptr(const Point & p, const Parameters &, const std::string & libmesh_dbg_var(sys_name), const std::string & unknown_name);

        static Gradient gptr(const Point & p, const Parameters &, const std::string & libmesh_dbg_var(sys_name), const std::string & unknown_name);

        private:

        // Local copy of the _from_system
        System & from_system;

        // Local copy of the _to_system
        System & to_system;

    };

    // This class provides the functor we will supply to System::project_vector
    // inside InterMeshProjection::project_system_vectors.
    // Object is constructed by passing in a pointer to the mesh function whose
    // gradient we want to shim via operator().
    class GradientMeshFunction : public FunctionBase<Gradient>
    {
        public:
        // Constructor
        GradientMeshFunction(MeshFunction * _mesh_function);

        // Destructor
        virtual ~GradientMeshFunction () { }

        virtual void init () { }

        virtual std::unique_ptr<FunctionBase<Gradient>> clone () const
        {
          //return libmesh_make_unique<GradientMeshFunction>(dynamic_cast<MeshFunction *>(mesh_function.get()));
          return libmesh_make_unique<GradientMeshFunction>(mesh_function.get());
        }

        virtual Gradient operator() (const Point & , const Real)
        { libmesh_not_implemented(); }

        virtual void operator() (const Point & p, const Real, DenseVector<Gradient> & output);

        private:

        // Local copy of the passed in mesh function.
        std::unique_ptr<MeshFunction> mesh_function;

    };
}

#endif // LIBMESH_INTER_MESH_PROJECTION_H
