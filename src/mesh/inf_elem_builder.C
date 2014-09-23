// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes
#include "libmesh/inf_elem_builder.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/cell_inf_prism12.h"
#include "libmesh/cell_inf_hex8.h"
#include "libmesh/cell_inf_hex16.h"
#include "libmesh/cell_inf_hex18.h"
#include "libmesh/mesh_base.h"

#ifdef DEBUG
#include "libmesh/parallel_mesh.h"
#endif

namespace libMesh
{

const Point InfElemBuilder::build_inf_elem(bool be_verbose)
{
  // determine origin automatically,
  // works only if the mesh has no symmetry planes.
  const MeshTools::BoundingBox b_box = MeshTools::bounding_box(_mesh);
  Point origin = (b_box.first + b_box.second) / 2;

  if (be_verbose && _mesh.processor_id() == 0)
    {
#ifdef DEBUG
      libMesh::out << " Determined origin for Infinite Elements:"
                   << std::endl
                   << "  ";
      origin.write_unformatted(libMesh::out);
      libMesh::out << std::endl;
#endif
    }

  // Call the protected implementation function with the
  // automatically determined origin.
  this->build_inf_elem(origin, false, false, false, be_verbose);

  // when finished with building the Ifems,
  // it remains to prepare the mesh for use:
  // find neighbors (again), partition (if needed)...
  this->_mesh.prepare_for_use (/*skip_renumber =*/ false);

  return origin;
}












const Point InfElemBuilder::build_inf_elem (const InfElemOriginValue& origin_x,
                                            const InfElemOriginValue& origin_y,
                                            const InfElemOriginValue& origin_z,
                                            const bool x_sym,
                                            const bool y_sym,
                                            const bool z_sym,
                                            const bool be_verbose,
                                            std::vector<const Node*>* inner_boundary_nodes)
{
  START_LOG("build_inf_elem()", "InfElemBuilder");

  // first determine the origin of the
  // infinite elements.  For this, the
  // origin defaults to the given values,
  // and may be overridden when the user
  // provided values
  Point origin(origin_x.second, origin_y.second, origin_z.second);

  // when only _one_ of the origin coordinates is _not_
  // given, we have to determine it on our own
  if ( !origin_x.first || !origin_y.first || !origin_z.first)
    {
      // determine origin
      const MeshTools::BoundingBox b_box = MeshTools::bounding_box(_mesh);
      const Point auto_origin = (b_box.first+b_box.second)/2;

      // override default values, if necessary
      if (!origin_x.first)
        origin(0) = auto_origin(0);
      if (!origin_y.first)
        origin(1) = auto_origin(1);
      if (!origin_z.first)
        origin(2) = auto_origin(2);

      if (be_verbose)
        {
          libMesh::out << " Origin for Infinite Elements:" << std::endl;

          if (!origin_x.first)
            libMesh::out << "  determined x-coordinate" << std::endl;
          if (!origin_y.first)
            libMesh::out << "  determined y-coordinate" << std::endl;
          if (!origin_z.first)
            libMesh::out << "  determined z-coordinate" << std::endl;

          libMesh::out << "  coordinates: ";
          origin.write_unformatted(libMesh::out);
          libMesh::out << std::endl;
        }
    }

  else if (be_verbose)

    {
      libMesh::out << " Origin for Infinite Elements:" << std::endl;
      libMesh::out << "  coordinates: ";
      origin.write_unformatted(libMesh::out);
      libMesh::out << std::endl;
    }



  // Now that we have the origin, check if the user provided an \p
  // inner_boundary_nodes.  If so, we pass a std::set to the actual
  // implementation of the build_inf_elem(), so that we can convert
  // this to the Node* vector
  if (inner_boundary_nodes != NULL)
    {
      // note that the std::set that we will get
      // from build_inf_elem() uses the index of
      // the element in this->_elements vector,
      // and the second entry is the side index
      // for this element.  Therefore, we do _not_
      // need to renumber nodes and elements
      // prior to building the infinite elements.
      //
      // However, note that this method here uses
      // node id's... Do we need to renumber?


      // Form the list of faces of elements which finally
      // will tell us which nodes should receive boundary
      // conditions (to form the std::vector<const Node*>)
      std::set< std::pair<dof_id_type,
        unsigned int> > inner_faces;


      // build infinite elements
      this->build_inf_elem(origin,
                           x_sym, y_sym, z_sym,
                           be_verbose,
                           &inner_faces);

      if (be_verbose)
        {
          this->_mesh.print_info();
          libMesh::out << "Data pre-processing:" << std::endl
                       << " convert the <int,int> list to a Node* list..."
                       << std::endl;
        }

      // First use a std::vector<dof_id_type> that holds
      // the global node numbers.  Then sort this vector,
      // so that it can be made unique (no multiple occurence
      // of a node), and then finally insert the Node* in
      // the vector inner_boundary_nodes.
      //
      // Reserve memory for the vector<> with
      // 4 times the size of the number of elements in the
      // std::set. This is a good bet for Quad4 face elements.
      // For higher-order elements, this probably _has_ to lead
      // to additional allocations...
      // Practice has to show how this affects performance.
      std::vector<dof_id_type> inner_boundary_node_numbers;
      inner_boundary_node_numbers.reserve(4*inner_faces.size());

      // Now transform the set of pairs to a list of (possibly
      // duplicate) global node numbers.
      std::set< std::pair<dof_id_type,unsigned int> >::iterator face_it = inner_faces.begin();
      const std::set< std::pair<dof_id_type,unsigned int> >::iterator face_end = inner_faces.end();
      for(; face_it!=face_end; ++face_it)
        {
          std::pair<dof_id_type,unsigned int> p = *face_it;

          // build a full-ordered side element to get _all_ the base nodes
          UniquePtr<Elem> side( this->_mesh.elem(p.first)->build_side(p.second) );

          // insert all the node numbers in inner_boundary_node_numbers
          for (unsigned int n=0; n< side->n_nodes(); n++)
            inner_boundary_node_numbers.push_back(side->node(n));
        }


      // inner_boundary_node_numbers now still holds multiple entries of
      // node numbers.  So first sort, then unique the vector.
      // Note that \p std::unique only puts the new ones in
      // front, while to leftovers are @e not deleted.  Instead,
      // it returns a pointer to the end of the unique range.
      //TODO:[BSK] int_ibn_size_before is not the same type as unique_size!
#ifndef NDEBUG
      const std::size_t ibn_size_before = inner_boundary_node_numbers.size();
#endif
      std::sort (inner_boundary_node_numbers.begin(), inner_boundary_node_numbers.end());
      std::vector<dof_id_type>::iterator unique_end =
        std::unique (inner_boundary_node_numbers.begin(), inner_boundary_node_numbers.end());

      std::size_t unique_size = std::distance(inner_boundary_node_numbers.begin(), unique_end);
      libmesh_assert_less_equal (unique_size, ibn_size_before);

      // Finally, create const Node* in the inner_boundary_nodes
      // vector.  Reserve, not resize (otherwise, the push_back
      // would append the interesting nodes, while NULL-nodes
      // live in the resize'd area...
      inner_boundary_nodes->reserve (unique_size);
      inner_boundary_nodes->clear();


      std::vector<dof_id_type>::iterator pos_it = inner_boundary_node_numbers.begin();
      for (; pos_it != unique_end; ++pos_it)
        {
          const Node& node = this->_mesh.node(*pos_it);
          inner_boundary_nodes->push_back(&node);
        }

      if (be_verbose)
        libMesh::out << "  finished identifying " << unique_size
                     << " target nodes." << std::endl;
    }

  else

    {
      // There are no inner boundary nodes, so simply build the infinite elements
      this->build_inf_elem(origin, x_sym, y_sym, z_sym, be_verbose);
    }


  STOP_LOG("build_inf_elem()", "InfElemBuilder");

  // when finished with building the Ifems,
  // it remains to prepare the mesh for use:
  // find neighbors again, partition (if needed)...
  this->_mesh.prepare_for_use (/*skip_renumber =*/ false);

  return origin;
}









// The actual implementation of building elements.
void InfElemBuilder::build_inf_elem(const Point& origin,
                                    const bool x_sym,
                                    const bool y_sym,
                                    const bool z_sym,
                                    const bool be_verbose,
                                    std::set< std::pair<dof_id_type,
                                    unsigned int> >* inner_faces)
{
  if (be_verbose)
    {
#ifdef DEBUG
      libMesh::out << " Building Infinite Elements:" << std::endl;
      libMesh::out << "  updating element neighbor tables..." << std::endl;
#else
      libMesh::out << " Verbose mode disabled in non-debug mode." << std::endl;
#endif
    }


  // update element neighbors
  this->_mesh.find_neighbors();

  START_LOG("build_inf_elem()", "InfElemBuilder");

  // A set for storing element number, side number pairs.
  // pair.first == element number, pair.second == side number
  std::set< std::pair<dof_id_type,unsigned int> > faces;
  std::set< std::pair<dof_id_type,unsigned int> > ofaces;

  // A set for storing node numbers on the outer faces.
  std::set<dof_id_type> onodes;

  // The distance to the farthest point in the mesh from the origin
  Real max_r=0.;

  // The index of the farthest point in the mesh from the origin
  int max_r_node = -1;

#ifdef DEBUG
  if (be_verbose)
    {
      libMesh::out << "  collecting boundary sides";
      if (x_sym || y_sym || z_sym)
        libMesh::out << ", skipping sides in symmetry planes..." << std::endl;
      else
        libMesh::out << "..." << std::endl;
    }
#endif

  // Iterate through all elements and sides, collect indices of all active
  // boundary sides in the faces set. Skip sides which lie in symmetry planes.
  // Later, sides of the inner boundary will be sorted out.
  {
    MeshBase::element_iterator       it  = this->_mesh.active_elements_begin();
    const MeshBase::element_iterator end = this->_mesh.active_elements_end();

    for(; it != end; ++it)
      {
        Elem* elem = *it;

        for (unsigned int s=0; s<elem->n_neighbors(); s++)
          {
            // check if elem(e) is on the boundary
            if (elem->neighbor(s) == NULL)
              {
                // note that it is safe to use the Elem::side() method,
                // which gives a non-full-ordered element
                UniquePtr<Elem> side(elem->build_side(s));

                // bool flags for symmetry detection
                bool sym_side=false;
                bool on_x_sym=true;
                bool on_y_sym=true;
                bool on_z_sym=true;


                // Loop over the nodes to check whether they are on the symmetry planes,
                // and therefore sufficient to use a non-full-ordered side element
                for(unsigned int n=0; n<side->n_nodes(); n++)
                  {
                    const Point dist_from_origin = this->_mesh.point(side->node(n)) - origin;

                    if(x_sym)
                      if( std::abs(dist_from_origin(0)) > 1.e-3 )
                        on_x_sym=false;

                    if(y_sym)
                      if( std::abs(dist_from_origin(1)) > 1.e-3 )
                        on_y_sym=false;

                    if(z_sym)
                      if( std::abs(dist_from_origin(2)) > 1.e-3 )
                        on_z_sym=false;

                    //       if(x_sym)
                    // if( std::abs(dist_from_origin(0)) > 1.e-6 )
                    //   on_x_sym=false;

                    //       if(y_sym)
                    // if( std::abs(dist_from_origin(1)) > 1.e-6 )
                    //   on_y_sym=false;

                    //       if(z_sym)
                    // if( std::abs(dist_from_origin(2)) > 1.e-6 )
                    //   on_z_sym=false;

                    //find the node most distant from origin

                    Real r = dist_from_origin.size();
                    if (r > max_r)
                      {
                        max_r = r;
                        max_r_node=side->node(n);
                      }

                  }

                sym_side = (x_sym && on_x_sym) || (y_sym && on_y_sym) || (z_sym && on_z_sym);

                if (!sym_side)
                  faces.insert( std::make_pair(elem->id(), s) );

              } // neighbor(s) == NULL
          } // sides
      } // elems
  }






  //  If a boundary side has one node on the outer boundary,
  //  all points of this side are on the outer boundary.
  //  Start with the node most distant from origin, which has
  //  to be on the outer boundary, then recursively find all
  //  sides and nodes connected to it. Found sides are moved
  //  from faces to ofaces, nodes are collected in onodes.
  //  Here, the search is done iteratively, because, depending on
  //  the mesh, a very high level of recursion might be necessary.
  if (max_r_node > 0)
    onodes.insert(max_r_node);


  {
    std::set< std::pair<dof_id_type,unsigned int> >::iterator face_it = faces.begin();
    unsigned int facesfound=0;
    while (face_it != faces.end()) {

      std::pair<dof_id_type, unsigned int> p;
      p = *face_it;

      // This has to be a full-ordered side element,
      // since we need the correct n_nodes,
      UniquePtr<Elem> side(this->_mesh.elem(p.first)->build_side(p.second));

      bool found=false;
      for(unsigned int sn=0; sn<side->n_nodes(); sn++)
        if(onodes.count(side->node(sn)))
          {
            found=true;
            break;
          }


      // If a new oface is found, include its nodes in onodes
      if(found)
        {
          for(unsigned int sn=0; sn<side->n_nodes(); sn++)
            onodes.insert(side->node(sn));

          ofaces.insert(p);
          ++face_it; // iteration is done here
          faces.erase(p);

          facesfound++;
        }

      else
        ++face_it; // iteration is done here

      // If at least one new oface was found in this cycle,
      // do another search cycle.
      if(facesfound>0 && face_it == faces.end())
        {
          facesfound = 0;
          face_it    = faces.begin();
        }

    }
  }


#ifdef DEBUG
  if (be_verbose)
    libMesh::out << "  found "
                 << faces.size()
                 << " inner and "
                 << ofaces.size()
                 << " outer boundary faces"
                 << std::endl;
#endif

  // When the user provided a non-null pointer to
  // inner_faces, that implies he wants to have
  // this std::set.  For now, simply copy the data.
  if (inner_faces != NULL)
    *inner_faces = faces;

  // free memory, clear our local variable, no need
  // for it any more.
  faces.clear();


  // outer_nodes maps onodes to their duplicates
  std::map<dof_id_type, Node *> outer_nodes;

  // We may need to pick our own object ids in parallel
  dof_id_type old_max_node_id = _mesh.max_node_id();
  dof_id_type old_max_elem_id = _mesh.max_elem_id();

  // for each boundary node, add an outer_node with
  // double distance from origin.
  std::set<dof_id_type>::iterator on_it = onodes.begin();
  for( ; on_it != onodes.end(); ++on_it)
    {
      Point p = (Point(this->_mesh.point(*on_it)) * 2) - origin;
      if (_mesh.is_serial())
        {
          // Add with a default id in serial
          outer_nodes[*on_it]=this->_mesh.add_point(p);
        }
      else
        {
          // Pick a unique id in parallel
          Node &bnode = _mesh.node(*on_it);
          dof_id_type new_id = bnode.id() + old_max_node_id;
          outer_nodes[*on_it] =
            this->_mesh.add_point(p, new_id,
                                  bnode.processor_id());
        }
    }


#ifdef DEBUG
  // for verbose, remember n_elem
  dof_id_type n_conventional_elem = this->_mesh.n_elem();
#endif


  // build Elems based on boundary side type
  std::set< std::pair<dof_id_type,unsigned int> >::iterator face_it = ofaces.begin();
  for( ; face_it != ofaces.end(); ++face_it)
    {
      // Shortcut to the pair being iterated over
      std::pair<dof_id_type,unsigned int> p = *face_it;

      // build a full-ordered side element to get the base nodes
      UniquePtr<Elem> side(this->_mesh.elem(p.first)->build_side(p.second));

      // create cell depending on side type, assign nodes,
      // use braces to force scope.
      bool is_higher_order_elem = false;

      Elem* el;
      switch(side->type())
        {
          // 3D infinite elements
          // TRIs
        case TRI3:
          el=new InfPrism6;
          break;

        case TRI6:
          el=new InfPrism12;
          is_higher_order_elem = true;
          break;

          // QUADs
        case QUAD4:
          el=new InfHex8;
          break;

        case QUAD8:
          el=new InfHex16;
          is_higher_order_elem = true;
          break;

        case QUAD9:
          el=new InfHex18;

          // the method of assigning nodes (which follows below)
          // omits in the case of QUAD9 the bubble node; therefore
          // we assign these first by hand here.
          el->set_node(16) = side->get_node(8);
          el->set_node(17) = outer_nodes[side->node(8)];
          is_higher_order_elem=true;
          break;

          // 2D infinite elements
        case EDGE2:
          el=new InfQuad4;
          break;

        case EDGE3:
          el=new InfQuad6;
          el->set_node(4) = side->get_node(2);
          break;

          // 1D infinite elements not supported
        default:
          libMesh::out << "InfElemBuilder::build_inf_elem(Point, bool, bool, bool, bool): "
                       << "invalid face element "
                       << std::endl;
          continue;
        }

      // In parallel, assign unique ids to the new element
      if (!_mesh.is_serial())
        {
          Elem *belem = _mesh.elem(p.first);
          el->processor_id() = belem->processor_id();
          // We'd better not have elements with more than 6 sides
          el->set_id (belem->id() * 6 + p.second + old_max_elem_id);
        }

      // assign vertices to the new infinite element
      const unsigned int n_base_vertices = side->n_vertices();
      for(unsigned int i=0; i<n_base_vertices; i++)
        {
          el->set_node(i                ) = side->get_node(i);
          el->set_node(i+n_base_vertices) = outer_nodes[side->node(i)];
        }


      // when this is a higher order element,
      // assign also the nodes in between
      if (is_higher_order_elem)
        {
          // n_safe_base_nodes is the number of nodes in \p side
          // that may be safely assigned using below for loop.
          // Actually, n_safe_base_nodes is _identical_ with el->n_vertices(),
          // since for QUAD9, the 9th node was already assigned above
          const unsigned int n_safe_base_nodes   = el->n_vertices();

          for(unsigned int i=n_base_vertices; i<n_safe_base_nodes; i++)
            {
              el->set_node(i+n_base_vertices)   = side->get_node(i);
              el->set_node(i+n_safe_base_nodes) = outer_nodes[side->node(i)];
            }
        }


      // add infinite element to mesh
      this->_mesh.add_elem(el);
    } // for


#ifdef DEBUG
  _mesh.libmesh_assert_valid_parallel_ids();

  if (be_verbose)
    libMesh::out << "  added "
                 << this->_mesh.n_elem() - n_conventional_elem
                 << " infinite elements and "
                 << onodes.size()
                 << " nodes to the mesh"
                 << std::endl
                 << std::endl;
#endif

  STOP_LOG("build_inf_elem()", "InfElemBuilder");
}

} // namespace libMesh





#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
