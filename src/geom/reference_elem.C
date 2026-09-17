// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/reference_elem.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/threads.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fuzzy_equals.h"

// C++ includes
#include <map>
#include <sstream>
#include <memory> // std::unique_ptr
#include <cmath>  // std::sqrt, std::cbrt
#include <utility> // std::pair
#include <vector>


//-----------------------------------------------
// anonymous namespace for implementation details
namespace
{
using namespace libMesh;

namespace ElemDataStrings
{
// GCC 5.2.0 warns about overlength strings in the auto-generated
// reference_elem.data file.
#pragma GCC diagnostic ignored "-Woverlength-strings"
#include "reference_elem.data"
#pragma GCC diagnostic warning "-Woverlength-strings"
}

typedef Threads::spin_mutex InitMutex;

// Mutex for thread safety.
InitMutex init_mtx;

// map from ElemType to reference element file system object name
typedef std::map<ElemType, const char *> FileMapType;
FileMapType ref_elem_file;
Elem * ref_elem_map[INVALID_ELEM];



class SingletonCache : public libMesh::Singleton
{
public:
  virtual ~SingletonCache() = default;

  std::vector<std::unique_ptr<Node>> node_list;
  std::vector<std::unique_ptr<Elem>> elem_list;
};

// From [0], regarding the lifetime of the singleton_cache variable:
//
//     "All variables at namespace level, including the anonymous
//     namespace and function local static variable have static
//     storage duration unless they are declared thread_local."
//
// Variables with static storage duration are destroyed at the end of
// the program execution. From [1],
//
//     "If it is a pointer to the data which is static ... then like all
//     other dynamically allocated data, it will only be destructed when
//     you delete it.  There are two frequent solutions:
//     * use a smart pointer, which has a destructor which deletes it, or
//     * don't delete it; in most cases, there's really no reason to call the
//       destructor, and if you happen to use the instance in the destructors
//       of other static objects, you'll run into an order of destruction
//       problem."
//
// I tried making the singleton_cache a std::unique_ptr, but this resulted
// in a segfault during program shutdown which appeared to come from the
// single_cache unique_ptr's destructor. I didn't investigate whether the
// issue was caused by a double deletion or what, but it appears that the
// first suggestion above may not be valid in general.
//
// [0]: https://stackoverflow.com/questions/24342393/how-anonymous-namespaces-avoids-making-global-static-variable
// [1]: https://stackoverflow.com/questions/6850009/c-deleting-static-data
SingletonCache * singleton_cache = nullptr;



void read_ref_elem (const ElemType type_in,
                    std::istream & in)
{
  libmesh_assert (singleton_cache != nullptr);

  std::string dummy;
  unsigned int n_elem, n_nodes, elem_type_read, nn;
  double x, y, z;

  in >> dummy;
  in >> n_elem;  /**/ std::getline (in, dummy); libmesh_assert_equal_to (n_elem, 1);
  in >> n_nodes; /**/ std::getline (in, dummy);
  in >> dummy;   /**/ std::getline (in, dummy);
  in >> dummy;   /**/ std::getline (in, dummy);
  in >> dummy;   /**/ std::getline (in, dummy);
  in >> dummy;   /**/ std::getline (in, dummy);
  in >> n_elem;  /**/ std::getline (in, dummy); libmesh_assert_equal_to (n_elem, 1);

  in >> elem_type_read;

  libmesh_assert_less (elem_type_read, INVALID_ELEM);
  libmesh_assert_equal_to (elem_type_read, static_cast<unsigned int>(type_in));
  libmesh_assert_equal_to (n_nodes, Elem::type_to_n_nodes_map[elem_type_read]);

  // Construct elem of appropriate type, store in the elem_list
  auto & uelem = singleton_cache->elem_list.emplace_back(Elem::build(type_in));

  // We are expecting an identity map, so assert it!
  for (unsigned int n=0; n<n_nodes; n++)
    {
      in >> nn;
      libmesh_assert_equal_to (n,nn);
    }

  for (unsigned int n=0; n<n_nodes; n++)
    {
      in >> x >> y >> z;

      auto & new_node =
        singleton_cache->node_list.emplace_back(Node::build(x,y,z,n));

      uelem->set_node(n, new_node.get());
    }

  // it is entirely possible we ran out of file or encountered
  // another error.  If so, throw an error.
  libmesh_error_msg_if(!in, "ERROR while creating element singleton!");

  // Also store a pointer to the newly created Elem in the ref_elem_map array.
  ref_elem_map[type_in] = uelem.get();
}



void init_ref_elem_table()
{
  // outside mutex - if this pointer is set, we can trust it.
  if (singleton_cache != nullptr)
    return;

  // playing with fire here - lock before touching shared
  // data structures
  InitMutex::scoped_lock lock(init_mtx);

  // inside mutex - pointer may have changed while waiting
  // for the lock to acquire, check it again.
  if (singleton_cache != nullptr)
    return;

  // OK, if we get here we have the lock and we are not
  // initialized.  populate singleton. Note that we do not
  // use a smart pointer to manage the singleton_cache variable
  // since variables with static storage duration are destroyed
  // automatically at the end of program execution.
  singleton_cache = new SingletonCache;

  // initialize the reference file table
  {
    ref_elem_file.clear();

    // 0D elements
    ref_elem_file[NODEELEM] = ElemDataStrings::one_nodeelem;

    // 1D elements
    ref_elem_file[EDGE2]    = ElemDataStrings::one_edge;
    ref_elem_file[EDGE3]    = ElemDataStrings::one_edge3;
    ref_elem_file[EDGE4]    = ElemDataStrings::one_edge4;

    // 2D elements
    ref_elem_file[TRI3]     = ElemDataStrings::one_tri;
    ref_elem_file[TRI6]     = ElemDataStrings::one_tri6;
    ref_elem_file[TRI7]     = ElemDataStrings::one_tri7;

    ref_elem_file[QUAD4]    = ElemDataStrings::one_quad;
    ref_elem_file[QUAD8]    = ElemDataStrings::one_quad8;
    ref_elem_file[QUAD9]    = ElemDataStrings::one_quad9;

    // 3D elements
    ref_elem_file[HEX8]     = ElemDataStrings::one_hex;
    ref_elem_file[HEX20]    = ElemDataStrings::one_hex20;
    ref_elem_file[HEX27]    = ElemDataStrings::one_hex27;

    ref_elem_file[TET4]     = ElemDataStrings::one_tet;
    ref_elem_file[TET10]    = ElemDataStrings::one_tet10;
    ref_elem_file[TET14]    = ElemDataStrings::one_tet14;

    ref_elem_file[PRISM6]   = ElemDataStrings::one_prism;
    ref_elem_file[PRISM15]  = ElemDataStrings::one_prism15;
    ref_elem_file[PRISM18]  = ElemDataStrings::one_prism18;
    ref_elem_file[PRISM20]  = ElemDataStrings::one_prism20;
    ref_elem_file[PRISM21]  = ElemDataStrings::one_prism21;

    ref_elem_file[PYRAMID5] = ElemDataStrings::one_pyramid;
    ref_elem_file[PYRAMID13] = ElemDataStrings::one_pyramid13;
    ref_elem_file[PYRAMID14] = ElemDataStrings::one_pyramid14;
    ref_elem_file[PYRAMID18] = ElemDataStrings::one_pyramid18;

    // No entry here for C0POLYGON or C0POLYHEDRON - the reference
    // element there depends on the number of nodes and can't be
    // precomputed.
  }

  // Read'em
  for (const auto & [elem_type, filename] : ref_elem_file)
    {
      std::istringstream stream(filename);
      read_ref_elem(elem_type, stream);
    }
}


// no reason to do this at startup -
// data structures will get initialized *if*
// ReferenceElem::get() is ever called.
// // Class to setup singleton data
// class ReferenceElemSetup : public Singleton::Setup
// {
//   void setup ()
//   {
//     init_ref_elem_table();
//   }
// } reference_elem_setup;

} // anonymous namespace



//----------------------------------------------------------------------------
// external API Implementation
namespace libMesh
{
namespace ReferenceElem
{
const Elem & get (const ElemType type_in)
{
  ElemType base_type = type_in;

  // For shell elements, use non shell type as the base type
  if (type_in == TRISHELL3)
    base_type = TRI3;

  if (type_in == QUADSHELL4)
    base_type = QUAD4;

  if (type_in == QUADSHELL8)
    base_type = QUAD8;

  if (type_in == QUADSHELL9)
    base_type = QUAD9;

  init_ref_elem_table();

  // Throw an error if the user asked for an ElemType that we don't
  // have a reference element for.
  libmesh_error_msg_if(ref_elem_map[base_type] == nullptr || type_in == INVALID_ELEM,
                       "No reference elem data available for ElemType " << type_in
                       << " = " << Utility::enum_to_string(type_in) << ".");

  return *ref_elem_map[base_type];
}

std::pair<std::unique_ptr<Elem>, std::vector<std::unique_ptr<Node>>>
ideal_target (const ElemType type)
{
  // Build target element
  auto target_elem = Elem::build(type);

  // Volume of reference element
  const auto ref_vol = target_elem->reference_elem()->volume();

  // Update the nodes of the target element, depending on type
  const Real sqrt_2 = std::sqrt(Real(2));
  const Real sqrt_3 = std::sqrt(Real(3));
  std::vector<std::unique_ptr<Node>> owned_nodes;

  const auto type_str = Utility::enum_to_string(type);

  // Elems deriving from Tri
  if (type_str.compare(0, 3, "TRI") == 0)
    {

      // The target element will be an equilateral triangle with area equal to
      // the area of the reference element.

      // Equilateral triangle side length preserving area of the reference element
      const auto side_length = std::sqrt(4. / sqrt_3 * ref_vol);

      // Define the nodal locations of the vertices
      const auto & s = side_length;
      //                                         x        y                  node_id
      owned_nodes.emplace_back(Node::build(Point(0.,      0.),               0));
      owned_nodes.emplace_back(Node::build(Point(s,       0.),               1));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * sqrt_3 * s), 2));

      switch (type)
        {
            case TRI3: {
              // Nothing to do here, vertices already added above
              break;
            }

            case TRI6: {
              // Define the midpoint nodes of the equilateral triangle
              //                                         x         y                   node_id
              owned_nodes.emplace_back(Node::build(Point(0.50 * s, 0.00),              3));
              owned_nodes.emplace_back(Node::build(Point(0.75 * s, 0.25 * sqrt_3 * s), 4));
              owned_nodes.emplace_back(Node::build(Point(0.25 * s, 0.25 * sqrt_3 * s), 5));

              break;
            }

          default:
            libmesh_error_msg("Unsupported triangular element: " << type_str);
            break;
        }
    } // if Tri

  // Elems deriving from Prism
  else if (type_str.compare(0, 5, "PRISM") == 0)
    {

      // The target element will be a prism with an equilateral triangular
      // base with volume equal to the volume of the reference element.

      // For an equilateral triangular base with side length s, the
      // base area is s^2 * sqrt(3) / 4.
      // The prism height that will result in equal face areas is
      // s * sqrt(3) / 4. We choose s such that the target element has
      // the same volume as the reference element:
      // v = (s^2 * sqrt(3) / 4) * (s * sqrt(3) / 4) = 3 * s^3 / 4
      // --> s = (16 * v / 3)^(1/3)
      // I have no particular motivation for imposing equal face areas,
      // so this can be updated if a more `optimal` target prism is
      // identified.

      // Side length that preserves the volume of the reference element
      const auto side_length = std::cbrt(16. * ref_vol / 3.);
      // Prism height with the property that all faces have equal area
      const auto target_height = 0.25 * side_length * sqrt_3;

      const auto & s = side_length;
      const auto & h = target_height;
      //                                         x        y                 z    node_id
      owned_nodes.emplace_back(Node::build(Point(0.,      0.,               0.), 0));
      owned_nodes.emplace_back(Node::build(Point(s,       0.,               0.), 1));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * sqrt_3 * s, 0.), 2));
      owned_nodes.emplace_back(Node::build(Point(0.,      0.,               h),  3));
      owned_nodes.emplace_back(Node::build(Point(s,       0.,               h),  4));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * sqrt_3 * s, h),  5));

      if (type == PRISM15 || type == PRISM18 || type == PRISM20 || type == PRISM21)
        {
          // Define the edge midpoint nodes of the prism
          const auto & on = owned_nodes;
          owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1]) / 2.), 6));
          owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[2]) / 2.), 7));
          owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[0]) / 2.), 8));
          owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[3]) / 2.), 9));
          owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[4]) / 2.), 10));
          owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[5]) / 2.), 11));
          owned_nodes.emplace_back(Node::build(Point((*on[3] + *on[4]) / 2.), 12));
          owned_nodes.emplace_back(Node::build(Point((*on[4] + *on[5]) / 2.), 13));
          owned_nodes.emplace_back(Node::build(Point((*on[5] + *on[3]) / 2.), 14));

          if (type == PRISM18 || type == PRISM20 || type == PRISM21)
            {
              // Define the rectangular face midpoint nodes of the prism
              owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1] + *on[3] + *on[4]) / 4.), 15));
              owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[2] + *on[4] + *on[5]) / 4.), 16));
              owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[2] + *on[3] + *on[5]) / 4.), 17));

              if (type == PRISM20 || type == PRISM21)
                {
                  // Define the triangular face midpoint nodes of the prism
                  owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1] + *on[2]) / 3.), 18));
                  owned_nodes.emplace_back(Node::build(Point((*on[3] + *on[4] + *on[5]) / 3.), 19));

                  if (type == PRISM21)
                    // Define the interior point of the prism
                    owned_nodes.emplace_back(Node::build(Point((*on[9] + *on[10] + *on[11]) / 3.), 20));

                }
            }
        }

      else if (type != PRISM6)
        libmesh_error_msg("Unsupported prism element: " << type_str);

    } // if Prism

  // Elems deriving from Pyramid
  else if (type_str.compare(0, 7, "PYRAMID") == 0)
    {

      // The target element is a pyramid with an square base and
      // equilateral triangular sides with volume equal to the volume of the
      // reference element.

      // A pyramid with square base sidelength s and equilateral triangular
      // sides has height h = s / sqrt(2).
      // The volume is v = s^2 h / 3 = s^3 / ( 3 sqrt(2)).
      // Solving for s: s = (3 sqrt(2) v)^(1/3), where v is the volume of the
      // non-optimal reference element.

      // Side length that preserves the volume of the reference element
      const auto side_length = std::cbrt(3. * sqrt_2 * ref_vol);
      // Pyramid height with the property that all faces are equilateral triangles
      const auto target_height = side_length / sqrt_2;

      const auto & s = side_length;
      const auto & h = target_height;

      //                                         x        y        z    node_id
      owned_nodes.emplace_back(Node::build(Point(0.,      0.,      0.), 0));
      owned_nodes.emplace_back(Node::build(Point(s,       0.,      0.), 1));
      owned_nodes.emplace_back(Node::build(Point(s,       s,       0.), 2));
      owned_nodes.emplace_back(Node::build(Point(0.,      s,       0.), 3));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * s, h),  4));

      if (type == PYRAMID13 || type == PYRAMID14 || type == PYRAMID18)
        {
          const auto & on = owned_nodes;
          // Define the edge midpoint nodes of the pyramid

          // Base node to base node midpoint nodes
          owned_nodes.emplace_back(Node::build((*on[0] + *on[1]) / 2., 5));
          owned_nodes.emplace_back(Node::build((*on[1] + *on[2]) / 2., 6));
          owned_nodes.emplace_back(Node::build((*on[2] + *on[3]) / 2., 7));
          owned_nodes.emplace_back(Node::build((*on[3] + *on[0]) / 2., 8));

          // Base node to apex node midpoint nodes
          owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[4]) / 2.), 9));
          owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[4]) / 2.), 10));
          owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[4]) / 2.), 11));
          owned_nodes.emplace_back(Node::build(Point((*on[3] + *on[4]) / 2.), 12));

          if (type == PYRAMID14 || type == PYRAMID18)
            {
              // Define the square face midpoint node of the pyramid
              owned_nodes.emplace_back(
                  Node::build(Point((*on[0] + *on[1] + *on[2] + *on[3]) / 4.), 13));

              if (type == PYRAMID18)
                {
                  // Define the triangular face nodes
                  owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1] + *on[4]) / 3.), 14));
                  owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[2] + *on[4]) / 3.), 15));
                  owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[3] + *on[4]) / 3.), 16));
                  owned_nodes.emplace_back(Node::build(Point((*on[3] + *on[0] + *on[4]) / 3.), 17));
                }
            }
        }

      else if (type != PYRAMID5)
        libmesh_error_msg("Unsupported pyramid element: " << type_str);

    } // if Pyramid

  // Elems deriving from Tet
  else if (type_str.compare(0, 3, "TET") == 0)
    {

      // The ideal target element is a a regular tet with equilateral
      // triangles for all faces, with volume equal to the volume of the
      // reference element.

      // The volume of a tet is given by v = b * h / 3, where b is the area of
      // the base face and h is the height of the apex node. The area of an
      // equilateral triangle with side length s is b = sqrt(3) s^2 / 4.
      // For all faces to have side length s, the height of the apex node is
      // h = sqrt(2/3) * s. Then the volume is v = sqrt(2) * s^3 / 12.
      // Solving for s, the side length that will preserve the volume of the
      // reference element is s = (6 * sqrt(2) * v)^(1/3), where v is the volume
      // of the non-optimal reference element (i.e., a right tet).

      // Side length that preserves the volume of the reference element
      const auto side_length = std::cbrt(6. * sqrt_2 * ref_vol);
      // tet height with the property that all faces are equilateral triangles
      const auto target_height = sqrt_2 / sqrt_3 * side_length;

      const auto & s = side_length;
      const auto & h = target_height;

      // For regular tet
      //                                         x        y                z     node_id
      owned_nodes.emplace_back(Node::build(Point(0.,      0.,               0.), 0));
      owned_nodes.emplace_back(Node::build(Point(s,       0.,               0.), 1));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, 0.5 * sqrt_3 * s, 0.), 2));
      owned_nodes.emplace_back(Node::build(Point(0.5 * s, sqrt_3 / 6. * s,  h),  3));

      if (type == TET10 || type == TET14)
        {
          const auto & on = owned_nodes;
          // Define the edge midpoint nodes of the tet

          // Base node to base node midpoint nodes
          owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1]) / 2.), 4));
          owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[2]) / 2.), 5));
          owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[0]) / 2.), 6));
          // Base node to apex node midpoint nodes
          owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[3]) / 2.), 7));
          owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[3]) / 2.), 8));
          owned_nodes.emplace_back(Node::build(Point((*on[2] + *on[3]) / 2.), 9));

          if (type == TET14)
            {
              // Define the face midpoint nodes of the tet
              owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1] + *on[2]) / 3.), 10));
              owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[1] + *on[3]) / 3.), 11));
              owned_nodes.emplace_back(Node::build(Point((*on[1] + *on[2] + *on[3]) / 3.), 12));
              owned_nodes.emplace_back(Node::build(Point((*on[0] + *on[2] + *on[3]) / 3.), 13));
            }
        }

      else if (type != TET4)
        libmesh_error_msg("Unsupported tet element: " << type_str);

    } // if Tet

  // Set the target_elem equal to the reference elem
  else
    for (const auto & node : target_elem->reference_elem()->node_ref_range())
      owned_nodes.emplace_back(Node::build(node, node.id()));

  // Set nodes of target element
  for (const auto & node_ptr : owned_nodes)
    target_elem->set_node(node_ptr->id(), node_ptr.get());

  libmesh_assert(relative_fuzzy_equals(target_elem->volume(), ref_vol, TOLERANCE));

  return std::make_pair(std::move(target_elem), std::move(owned_nodes));
}
} // namespace ReferenceElem
} // namespace libMesh
