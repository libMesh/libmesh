// $Id: mesh_exodus_support.h,v 1.9 2003-10-02 03:39:25 jwpeterson Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002-2003  Benjamin S. Kirk, John W. Peterson
  
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



//--------------------------------------------------------
// Class for interfacing with Exodus II API

#ifndef __mesh_exodus_support_h__
#define __mesh_exodus_support_h__

#include "libmesh_common.h"

#ifdef HAVE_EXODUS_API

// C++ includes
#include <vector>

// local includes
#include "enum_elem_type.h"

namespace exII {
extern "C" {
#include "exodusII.h" // MAX_LINE_LENGTH and MAX_STR_LENGTH
}
}


  /**
   * This is the \p ExodusII class.
   * This class hides the implementation
   * details of interfacing with
   * the Exodus binary format.
   *
   * @author Johw W. Peterson, 2002.
   */

  
class ExodusII
{
 public:

  /**
   * Constructor. Automatically
   * initializes all the private
   * members of the class.  Also
   * allows you to set the verbosity
   * level to v=1 (on) or v=0 (off).
   */
  ExodusII(int v=0) : verbose(v),
                      comp_ws(sizeof(double)),
                      io_ws(0),
                      ex_id(0),
                      ex_err(0),
                      num_dim(0),
                      num_nodes(0),
                      num_elem(0),
                      num_elem_blk(0),
                      num_node_sets(0),
                      num_side_sets(0),
                      num_elem_this_blk(0),
                      num_nodes_per_elem(0),
                      num_attr(0),
                      req_info(0),
                      ret_int(0),
                      num_elem_all_sidesets(0),
                      ex_version(0.0),
                      ret_float(0.0),
                      ret_char(0),
                      title(new char[MAX_LINE_LENGTH]),
                      elem_type(new char[MAX_STR_LENGTH])
    {}

  /**
   * Destructor.  The only memory
   * allocated is for \p title and
   * \p elem_type.  This memory
   * is freed in the destructor.
   */
  ~ExodusII();

  
  /**
   * @returns the \p ExodusII
   * mesh dimension.
   */
  int get_num_dim()                const { return num_dim; }

  /**
   * @returns the total number of
   * nodes in the \p ExodusII mesh.
   */
  int get_num_nodes()              const { return num_nodes; }

  
  /**
   * @returns the total number of
   * elements in the \p ExodusII mesh.
   */
  int get_num_elem()               const { return num_elem; }

  /**
   * @returns the total number
   * of element blocks in
   * the \p ExodusII mesh.
   */
  int get_num_elem_blk()           const { return num_elem_blk; }

  /**
   * For a given block,
   * returns the total number
   * of elements.
   */
  int get_num_elem_this_blk()      const { return num_elem_this_blk; }

  /**
   * @returns the number of
   * nodes per element in
   * a given block. e.g.
   * for HEX27 it returns 27.
   */
  int get_num_nodes_per_elem()     const { return num_nodes_per_elem; }

  /**
   * @returns the total number
   * of sidesets in the \p ExodusII
   * mesh.  Each sideset contains
   * only one type of element.
   */
  int get_num_side_sets()          const { return num_side_sets; }

  /**
   * @returns the number of
   * elements in all the sidesets.
   * Effectively returns the
   * total number of elements
   * on the \p ExodusII mesh boundary.
   */
  int get_num_elem_all_sidesets()  const { return num_elem_all_sidesets; }

  /**
   * @returns the \f$ i^{th} \f$
   * node number in the
   * element connectivity
   * list for a given element.
   */
  int get_connect(int i)           const { return connect[i]; }

  /**
   * For a single sideset,
   * returns the total number of
   * elements in the sideset.
   */
  int get_num_sides_per_set(int i) const { return num_sides_per_set[i]; }

  /**
   * @returns the \f$ i^{th} \f$ entry
   * in the element list.
   * The element list contains
   * the numbers of all elements
   * on the boundary.
   */
  int get_elem_list(int i)         const { return elem_list[i]; }

  /**
   * @return a constant reference to the \p elem_list.   
   */
  const std::vector<int>& get_elem_list() const { return elem_list; }
  
  /**
   * @returns the \f$ i^{th} \f$ entry in
   * the side list.  This is
   * effectively the "side"
   * (face in 3D or edge in
   * 2D) number which lies
   * on the boundary.
   */
  int get_side_list(int i)         const { return side_list[i]; }

  /**
   * @return a constant reference to the \p side_list.   
   */
  const std::vector<int>& get_side_list() const { return side_list; }
  
  /**
   * @returns the \f$ i^{th} \f$ entry in
   * the id list.  This is the id
   * for the ith face on the boundary.
   */
  int get_id_list(int i)         const { return id_list[i]; }

  /**
   * @return a constant reference to the \p id_list.   
   */
  const std::vector<int>& get_id_list() const { return id_list; }
  
  /**
   * @returns the current
   * element type.  Note:
   * the default behavior
   * is for this value
   * to be in all capital
   * letters, e.g. \p HEX27.
   */
  char* get_elem_type()            const { return elem_type; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's x-coordinate.
   */
  double get_x(int i) const { return x[i]; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's y-coordinate.
   */
  double get_y(int i) const { return y[i]; }

  /**
   * @returns the \f$ i^{th} \f$
   * node's z-coordinate.
   */
  double get_z(int i) const { return z[i]; }

  /**
   * Opens an \p ExodusII mesh
   * file named \p filename
   * for reading.
   */
  void open(const char* filename);

  /**
   * Reads an \p ExodusII mesh
   * file header.
   */
  void read_header();

  /**
   * Prints the \p ExodusII
   * mesh file header,
   * which includes the
   * mesh title, the number
   * of nodes, number of
   * elements, mesh dimension,
   * and number of sidesets.
   */
  void print_header();

  /**
   * Reads the nodal data
   * (x,y,z coordinates)
   * from the \p ExodusII mesh
   * file.
   */
  void read_nodes();

  /**
   * Prints the nodal information
   * to \p std::cout.
   */
  void print_nodes();

  /**
   * Reads information for
   * all of the blocks in
   * the \p ExodusII mesh file.
   */
  void read_block_info();

  /**
   * Reads all of the element
   * connectivity for
   * block \p{block} in the
   * \p ExodusII mesh file.
   */
  void read_elem_in_block(int block);

  /**
   * Reads information about
   * all of the sidesets in
   * the \p ExodusII mesh file.
   */
  void read_sideset_info();

  /**
   * Reads information about
   * sideset \p id and
   * inserts it into the global
   * sideset array at the
   * position \p offset.
   */
  void read_sideset(int id, int offset);

  /**
   * Prints information
   * about all the sidesets.
   */
  void print_sideset_info();

  /**
   * Closes the \p ExodusII
   * mesh file.
   */
  void close();


  
  /**
   * This is the \p ExodusII
   * Conversion class.
   *
   * It provides a 
   * data structure
   * which contains
   * \p ExodusII node/edge maps
   * and name conversions.
   */

  class Conversion
    {
    public:

      /**
       * Constructor.  Initializes the const private member
       * variables.
       */
      Conversion(const int* nm, const int* sm, const ElemType ct) 
	: node_map(nm),       // Node map for this element
	side_map(sm),
	canonical_type(ct)    // Element type name in this code
	{}

      /**
       * Returns the ith component of the node map for this
       * element.  The node map maps the exodusII node numbering
       * format to this library's format.
       */
      int get_node_map(int i)          const { return node_map[i]; }

      /**
       * Returns the ith component of the side map for this
       * element.  The side map maps the exodusII side numbering
       * format to this library's format.
       */
      int get_side_map(int i)          const { return side_map[i]; }

      /**
       * Returns the canonical element type for this
       * element.  The canonical element type is the standard
       * element type understood by this library.
       */
      ElemType get_canonical_type()    const { return canonical_type; }

    private:
      /**
       * Pointer to the node map for this element.
       */
      const int* node_map;
      
      /**
       * Pointer to the side map for this element.
       */
      const int* side_map;

      /**
       * The canonical (i.e. standard for this library)
       * element type.
       */
      const ElemType canonical_type;
    };


      
  /**
   * This is the \p ExodusII
   * ElementMap class.
   *
   * It contains constant
   * maps between the \p ExodusII
   * naming/numbering schemes
   * and the canonical schemes
   * used in this code.
   */
  
  class ElementMaps
    {
    public:

      /**
       * Constructor.
       */
      ElementMaps() {}

      /**
       * 2D node maps.  These define
       * mappings from ExodusII-formatted
       * element numberings.
       */

      /**
       * The Quad4 node map.
       * Use this map for bi-linear
       * quadrilateral elements in 2D.
       */
      static const int quad4_node_map[4]; 

      /**
       * The Quad8 node map.
       * Use this map for serendipity
       * quadrilateral elements in 2D.
       */
      static const int quad8_node_map[8]; 

      /**
       * The Quad9 node map.
       * Use this map for bi-quadratic
       * quadrilateral elements in 2D.
       */
      static const int quad9_node_map[9]; 

      /**
       * The Tri3 node map.
       * Use this map for linear
       * triangles in 2D.
       */
      static const int tri3_node_map[3];     

      /**
       * The Tri6 node map.
       * Use this map for quadratic
       * triangular elements in 2D.
       */
      static const int tri6_node_map[6];     

      /**
       * 2D edge maps
       */

      /**
       * Maps the Exodus edge numbering for triangles.
       * Useful for reading sideset information.
       */ 
      static const int tri_edge_map[3]; 

      /**
       * Maps the Exodus edge numbering for quadrilaterals.
       * Useful for reading sideset information.
       */
      static const int quad_edge_map[4]; 

      
      
      /**
       * 3D maps.  These define
       * mappings from ExodusII-formatted
       * element numberings.
       */

      /**
       * The Hex8 node map.
       * Use this map for bi-linear
       * hexahedral elements in 3D.
       */
      static const int hex8_node_map[8]; 

      /**
       * The Hex20 node map.
       * Use this map for serendipity
       * hexahedral elements in 3D.
       */
      static const int hex20_node_map[20]; 

      /**
       * The Hex27 node map.
       * Use this map for bi-quadratic
       * hexahedral elements in 3D.
       */
      static const int hex27_node_map[27]; 

      /**
       * The Tet4 node map.
       * Use this map for linear
       * tetrahedral elements in 3D.
       */
      static const int tet4_node_map[4]; 


      /**
       * The Tet10 node map.
       * Use this map for quadratic
       * tetrahedral elements in 3D.
       */
      static const int tet10_node_map[10]; 


      /**
       * 3D face maps
       */
      
      /**
       * Maps the Exodus face numbering for general hexahedrals.
       * Useful for reading sideset information.
       */
      static const int hex_face_map[6];
      
      /**
       * Maps the Exodus face numbering for 27-noded hexahedrals.
       * Useful for reading sideset information.
       */
      static const int hex27_face_map[6];
      
      /**
       * Maps the Exodus face numbering for general tetrahedrals.
       * Useful for reading sideset information.
       */
      static const int tet_face_map[4];
      

      /**
       * @returns a conversion object given an element type name.
       */
      const Conversion assign_conversion(const std::string type);
      
      /**
       * @returns a conversion object given an element type.
       */
      const Conversion assign_conversion(const ElemType type);
    };


  
 private:

  
  /**
   * All of the \p ExodusII
   * API functions return
   * an \p{int} error value.
   * This function checks
   * to see if the error has
   * been set, and if it has,
   * prints the error message
   * contained in \p msg.
   */
  void check_err(const int error, const std::string msg);

  /**
   * Prints the message defined
   * in \p msg to \p std::cout.
   * Can be turned off if
   * verbosity is set to 0.
   */
  void message(const std::string msg);

  /**
   * Prints the message defined
   * in \p msg to \p std::cout
   * and appends the number
   * \p i to the end of the
   * message.  Useful for
   * printing messages in loops.
   * Can be turned off if
   * verbosity is set to 0.
   */
  void message(const std::string msg, int i);

  int   verbose;                       // On/Off message flag, 0 = no messages
  int   comp_ws;                       // ?
  int   io_ws;                         // ?
  int   ex_id;                         // File identification flag
  int   ex_err;                        // General error flag
  int   num_dim;                       // Number of dimensions in the mesh
  int   num_nodes;                     // Total number of nodes in the mesh
  int   num_elem;                      // Total number of elements in the mesh
  int   num_elem_blk;                  // Total number of element blocks
  int   num_node_sets;                 // Total number of node sets
  int   num_side_sets;                 // Total number of element sets
  int   num_elem_this_blk;             // Number of elements in this block
  int   num_nodes_per_elem;            // Number of nodes in each element
  int   num_attr;                      // Number of attributes for a given block
  int   req_info;                      // Generic required info tag
  int   ret_int;                       // Generic int returned by ex_inquire
  int   num_elem_all_sidesets;         // Total number of elements in all side sets
  std::vector<int> block_ids;          // Vector of the block identification numbers
  std::vector<int> connect;            // Vector of nodes in an element
  std::vector<int> ss_ids;             // Vector of the sideset IDs
  std::vector<int> num_sides_per_set;  // Number of sides (edges/faces) in current set
  std::vector<int> num_df_per_set;     // Number of distribution factors per set
  std::vector<int> elem_list;          // List of element numbers in all sidesets
  std::vector<int> side_list;          // Side (face/edge) number actually on the boundary 
  std::vector<int> id_list;            // Side (face/edge) id number
  float ex_version;                    // Version of Exodus you are using
  float ret_float;                     // Generic float returned by ex_inquire
  std::vector<double> x;               // x locations of node points
  std::vector<double> y;               // y locations of node points
  std::vector<double> z;               // z locations of node points
  char    ret_char;                    // Generic char returned by ex_inquire
  char*   title;                       // Problem title
  char*   elem_type;                   // Type of element in a given block

};




#endif
#endif
