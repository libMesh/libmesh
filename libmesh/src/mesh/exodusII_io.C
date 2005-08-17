// $Id: exodusII_io.C,v 1.10 2005-08-17 21:20:22 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <fstream>

// Local includes
#include "exodusII_io.h"
#include "boundary_info.h"
#include "mesh_base.h"
#include "enum_elem_type.h"
#include "elem.h"

// Wrap all the helper classes in an #ifdef to avoid excessive compilation
// time in the case of no ExodusII support

#ifdef HAVE_EXODUS_API

namespace exII {
  extern "C" {
#include "exodusII.h" // defines MAX_LINE_LENGTH, MAX_STR_LENGTH used later
  }
}

//-----------------------------------------------------------------------------
// Exodus class - private helper class defined in an anonymous namespace
namespace
{
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
    ExodusII(const bool v=false) : verbose(v),
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

//     /**
//      * @returns the number of
//      * elements in all the sidesets.
//      * Effectively returns the
//      * total number of elements
//      * on the \p ExodusII mesh boundary.
//      */
//     int get_num_elem_all_sidesets()  const { return num_elem_all_sidesets; }

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

//     /**
//      * @returns the \f$ i^{th} \f$ entry
//      * in the element list.
//      * The element list contains
//      * the numbers of all elements
//      * on the boundary.
//      */
//     int get_elem_list(int i)         const { return elem_list[i]; }

    /**
     * @return a constant reference to the \p elem_list.   
     */
    const std::vector<int>& get_elem_list() const { return elem_list; }
  
//     /**
//      * @returns the \f$ i^{th} \f$ entry in
//      * the side list.  This is
//      * effectively the "side"
//      * (face in 3D or edge in
//      * 2D) number which lies
//      * on the boundary.
//      */
//     int get_side_list(int i)         const { return side_list[i]; }

    /**
     * @return a constant reference to the \p side_list.   
     */
    const std::vector<int>& get_side_list() const { return side_list; }
  
//     /**
//      * @returns the \f$ i^{th} \f$ entry in
//      * the id list.  This is the id
//      * for the ith face on the boundary.
//      */
//     int get_id_list(int i)         const { return id_list[i]; }

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
     * block \p block in the
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


  
    //-------------------------------------------------------------------------
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


      
    //-------------------------------------------------------------------------
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
       * The Prism6 node map.
       * Use this map for quadratic
       * tetrahedral elements in 3D.
       */
      static const int prism6_node_map[6]; 


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
       * Maps the Exodus face numbering for general prisms.
       * Useful for reading sideset information.
       */
      static const int prism_face_map[5];
      

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
     * an \p int error value.
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

    const bool verbose;                  // On/Off message flag
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


  // ------------------------------------------------------------
  // ExodusII::ElementMaps static data

  // 2D node map definitions
  const int ExodusII::ElementMaps::quad4_node_map[4] = {0, 1, 2, 3};
  const int ExodusII::ElementMaps::quad8_node_map[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  const int ExodusII::ElementMaps::quad9_node_map[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  const int ExodusII::ElementMaps::tri3_node_map[3]  = {0, 1, 2};
  const int ExodusII::ElementMaps::tri6_node_map[6]  = {0, 1, 2, 3, 4, 5};

  // 2D edge map definitions
  const int ExodusII::ElementMaps::tri_edge_map[3] = {0, 1, 2};
  const int ExodusII::ElementMaps::quad_edge_map[4] = {0, 1, 2, 3};

  // 3D node map definitions
  const int ExodusII::ElementMaps::hex8_node_map[8]   = {0, 1, 2, 3, 4, 5, 6, 7};
  const int ExodusII::ElementMaps::hex20_node_map[20] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
							  10, 11, 12, 13, 14, 15, 16, 17, 18, 19};  
  const int ExodusII::ElementMaps::hex27_node_map[27] = { 1,  5,  6,  2,  0,  4,  7,  3, 13, 17, 14,  9,  8, 16,
							  18, 10, 12, 19, 15, 11, 24, 25, 22, 26, 21, 23, 20};  
  const int ExodusII::ElementMaps::tet4_node_map[4]   = {0, 1, 2, 3};
  const int ExodusII::ElementMaps::tet10_node_map[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  const int ExodusII::ElementMaps::prism6_node_map[6] = {0, 1, 2, 3, 4, 5};
  
  // 3D face map definitions
  const int ExodusII::ElementMaps::tet_face_map[4] =   {1, 2, 3, 0};
  const int ExodusII::ElementMaps::hex_face_map[6] =   {1, 2, 3, 4, 0, 5};
  const int ExodusII::ElementMaps::hex27_face_map[6] = {1, 0, 3, 5, 4, 2};
  const int ExodusII::ElementMaps::prism_face_map[5] = {-1,-1,-1,-1,-1}; // Not Implemented!


  // ------------------------------------------------------------
  // ExodusII class members
  ExodusII::~ExodusII()
  {
    delete [] title;
    delete [] elem_type;
  }



  void ExodusII::check_err(const int err, const std::string msg)
  {
    if (err < 0)
      {
	std::cout << msg << std::endl;
	error();
      }
  }



  void ExodusII::message(const std::string msg)
  {
    if (verbose) std::cout << msg << std::endl;
  }



  void ExodusII::message(const std::string msg, int i)
  {
    if (verbose) std::cout << msg << i << "." << std::endl;
  }



  void ExodusII::open(const char* filename)
  {
    ex_id = exII::ex_open(filename,
			  EX_READ,
			  &comp_ws,
			  &io_ws,
			  &ex_version);
  
    check_err(ex_id, "Error opening ExodusII mesh file.");
    if (verbose) std::cout << "File opened successfully." << std::endl;
  }



  void ExodusII::read_header()
  {
    ex_err = exII::ex_get_init(ex_id,
			       title,
			       &num_dim,
			       &num_nodes,
			       &num_elem,
			       &num_elem_blk,
			       &num_node_sets,
			       &num_side_sets);

    check_err(ex_err, "Error retrieving header info.");
    message("Exodus header info retrieved successfully.");
  }




  void ExodusII::print_header()
  {
    if (verbose)
      std::cout << "Title: \t" << title << std::endl
		<< "Mesh Dimension: \t"   << num_dim << std::endl
		<< "Number of Nodes: \t" << num_nodes << std::endl
		<< "Number of elements: \t" << num_elem << std::endl
		<< "Number of elt blocks: \t" << num_elem_blk << std::endl
		<< "Number of node sets: \t" << num_node_sets << std::endl
		<< "Number of side sets: \t" << num_side_sets << std::endl;
  }



  void ExodusII::read_nodes()
  {
    x.resize(num_nodes);
    y.resize(num_nodes); 
    z.resize(num_nodes); 

    ex_err = exII::ex_get_coord(ex_id,
				static_cast<void*>(&x[0]),
				static_cast<void*>(&y[0]),
				static_cast<void*>(&z[0]));
  
    check_err(ex_err, "Error retrieving nodal data.");
    message("Nodal data retrieved successfully."); 
  }



  void ExodusII::print_nodes()
  {
    for (int i=0; i<num_nodes; i++)
      {
	std::cout << "(" << x[i] << ", " << y[i] << ", " << z[i] << ")" << std::endl;
      }
  }



  void ExodusII::read_block_info()
  {
    block_ids.resize(num_elem_blk);
    ex_err = exII::ex_get_elem_blk_ids(ex_id,
				       &block_ids[0]); // Get all element block IDs.
    // Usually, there is only one
    // block since there is only
    // one type of element.
    // However, there could be more.

    check_err(ex_err, "Error getting block IDs.");
    message("All block IDs retrieved successfully."); 
  }



  void ExodusII::read_elem_in_block(int block)
  {
    ex_err = exII::ex_get_elem_block(ex_id,
				     block_ids[block],
				     elem_type,
				     &num_elem_this_blk,
				     &num_nodes_per_elem,
				     &num_attr);
    check_err(ex_err, "Error getting block info.");
    message("Info retrieved successfully for block: ", block); 
  
  
  
    // Read in the connectivity of the elements of this block
    connect.resize(num_nodes_per_elem*num_elem_this_blk);
    ex_err = exII::ex_get_elem_conn(ex_id,
				    block_ids[block],
				    &connect[0]);
  
    check_err(ex_err, "Error reading block connectivity.");
    message("Connectivity retrieved successfully for block: ", block); 
  }



  void ExodusII::read_sideset_info()
  {
    ss_ids.resize(num_side_sets);
    if (num_side_sets > 0)
      {
	ex_err = exII::ex_get_side_set_ids(ex_id,
					   &ss_ids[0]);
	check_err(ex_err, "Error retrieving sideset information.");
	message("All sideset information retrieved successfully."); 

	// Resize appropriate data structures -- only do this once outside the loop
	num_sides_per_set.resize(num_side_sets);
	num_df_per_set.resize(num_side_sets);
  
	req_info = EX_INQ_SS_ELEM_LEN; // Inquire about the length of the
	// concatenated side sets element list
	ex_err = exII::ex_inquire(ex_id,
				  req_info,
				  &ret_int,
				  &ret_float,
				  &ret_char);
	check_err(ex_err, "Error inquiring about side set element list length.");

	//std::cout << "Value returned by ex_inquire was: " << ret_int << std::endl;
	num_elem_all_sidesets = ret_int;	
	elem_list.resize (num_elem_all_sidesets);
	side_list.resize (num_elem_all_sidesets);
	id_list.resize   (num_elem_all_sidesets);
      }
  }


  void ExodusII::read_sideset(int id, int offset)
  {
    ex_err = exII::ex_get_side_set_param(ex_id,
					 ss_ids[id],
					 &num_sides_per_set[id],
					 &num_df_per_set[id]);
    check_err(ex_err, "Error retrieving sideset parameters.");
    message("Parameters retrieved successfully for sideset: ", id);

    ex_err = exII::ex_get_side_set(ex_id,
				   ss_ids[id],
				   &elem_list[offset],
				   &side_list[offset]);
    check_err(ex_err, "Error retrieving sideset data.");
    message("Data retrieved successfully for sideset: ", id);

    for (int i=0; i<num_sides_per_set[id]; i++)
      id_list[i+offset] = ss_ids[id];
  }



  void ExodusII::print_sideset_info()
  {
    for (int i=0; i<num_elem_all_sidesets; i++)
      {
	std::cout << elem_list[i] << " " << side_list[i] << std::endl;
      }
  }



  void ExodusII::close()
  {
    ex_err = exII::ex_close(ex_id);
    check_err(ex_err, "Error closing Exodus file.");
    message("Exodus file closed successfully."); 
  }



  // ------------------------------------------------------------
  // ExodusII::Conversion class members
  const ExodusII::Conversion ExodusII::ElementMaps::assign_conversion(const std::string type)
  {
    if (type == "QUAD4") 
      return assign_conversion(QUAD4);

    else if (type == "QUAD8")
      return assign_conversion(QUAD8);

    else if (type == "QUAD9")
      return assign_conversion(QUAD9);

    else if (type == "TRI3")
      return assign_conversion(TRI3);

    else if (type == "TRI6")
      return assign_conversion(TRI6);

    else if ((type == "HEX8") || (type == "HEX")) 
      return assign_conversion(HEX8);

    else if (type == "HEX20")
      return assign_conversion(HEX20);

    else if (type == "HEX27")
      return assign_conversion(HEX27);

    else if (type == "TETRA4")
      return assign_conversion(TET4);

    else if (type == "TETRA10")
      return assign_conversion(TET10);

    else if (type == "WEDGE")
      return assign_conversion(PRISM6);

    else
      {
	std::cerr << "ERROR! Unrecognized element type: " << type << std::endl;
	error();
      }

    error();
  
    const Conversion conv(tri3_node_map, tri_edge_map, TRI3); // dummy
    return conv;  
  }



  const ExodusII::Conversion ExodusII::ElementMaps::assign_conversion(const ElemType type)
  {
    switch (type)
      {

      case QUAD4:
	{
	  const Conversion conv(quad4_node_map, quad_edge_map, QUAD4);
	  return conv;
	}

      case QUAD8:
	{
	  const Conversion conv(quad8_node_map, quad_edge_map, QUAD8);
	  return conv;
	}
      
      case QUAD9:
	{
	  const Conversion conv(quad9_node_map, quad_edge_map, QUAD9);
	  return conv;
	}
      
      case TRI3:
	{
	  const Conversion conv(tri3_node_map, tri_edge_map, TRI3);
	  return conv;
	}
      
      case TRI6:
	{
	  const Conversion conv(tri6_node_map, tri_edge_map, TRI6);
	  return conv;
	}
      
      case HEX8:
	{
	  const Conversion conv(hex8_node_map, hex_face_map, HEX8);
	  return conv;
	}
      
      case HEX20:
	{
	  const Conversion conv(hex20_node_map, hex_face_map, HEX20);
	  return conv;
	}
      
      case HEX27:
	{
	  const Conversion conv(hex27_node_map, hex27_face_map, HEX27);
	  return conv;
	}
      
      case TET4:
	{
	  const Conversion conv(tet4_node_map, tet_face_map, TET4);
	  return conv;
	}
      
      case TET10:
	{
	  const Conversion conv(tet10_node_map, tet_face_map, TET10);
	  return conv;
	}

      case PRISM6:
	{
	  const Conversion conv(prism6_node_map, prism_face_map, PRISM6);
	  return conv;
	}
	
      default:
	error();
      }
    
    error();
    
    const Conversion conv(tri3_node_map, tri_edge_map, TRI3); // dummy
    return conv;  
  }
  
    
} // end anonymous namespace

#endif // #ifdef HAVE_EXODUS_API




// ------------------------------------------------------------
// ExodusII_IO class members
void ExodusII_IO::read (const std::string& fname)
{
#ifndef HAVE_EXODUS_API

  std::cerr <<  "ERROR, ExodusII API is not defined.\n"
	    << "Input file " << fname << " cannot be read"
	    << std::endl;
  error();
    
#else
  
  // Get a reference to the mesh we are reading
  MeshBase& mesh = MeshInput<MeshBase>::mesh();
  
  // Clear any existing mesh data
  mesh.clear();
  
  assert(mesh.mesh_dimension() != 1); // No support for 1D ExodusII meshes
  
#ifdef DEBUG
    this->verbose() = true;
#endif
  
  ExodusII ex(this->verbose()); // Instantiate ExodusII interface
  ExodusII::ElementMaps em;     // Instantiate the ElementMaps interface
    
  ex.open(fname.c_str());       // Open the exodus file, if possible
  ex.read_header();             // Get header information from exodus file
  ex.print_header();            // Print header information

  assert(static_cast<unsigned int>(ex.get_num_dim()) == mesh.mesh_dimension()); // Be sure number of dimensions
                                                                                // is equal to the number of 
                                                                                // dimensions in the mesh supplied.

  ex.read_nodes();                        // Read nodes from the exodus file
  mesh.reserve_nodes(ex.get_num_nodes()); // Reserve space for the nodes.
  
  // Loop over the nodes, create Nodes.
  for (int i=0; i<ex.get_num_nodes(); i++)
    mesh.add_point (Point(ex.get_x(i),
			  ex.get_y(i),
			  ex.get_z(i)));

  assert (static_cast<unsigned int>(ex.get_num_nodes()) == mesh.n_nodes());

  ex.read_block_info();                 // Get information about all the blocks
  mesh.reserve_elem(ex.get_num_elem()); // Reserve space for the elements
   

  // Read in the element connectivity for each block.
  int nelem_last_block = 0;

  // Loop over all the blocks
  for (int i=0; i<ex.get_num_elem_blk(); i++)
    {
      // Read the information for block i
      ex.read_elem_in_block (i);
      
      // Set any relevant node/edge maps for this element
      const std::string type (ex.get_elem_type());
      const ExodusII::Conversion conv = em.assign_conversion(type); 
      
      // Loop over all the faces in this block
      int jmax = nelem_last_block+ex.get_num_elem_this_blk();
      for (int j=nelem_last_block; j<jmax; j++)
	{
	  Elem* elem = mesh.add_elem (Elem::build (conv.get_canonical_type()).release());
	  assert (elem != NULL);
	    
	  // Set all the nodes for this element
	  for (int k=0; k<ex.get_num_nodes_per_elem(); k++)
	    {
	      int gi = (j-nelem_last_block)*ex.get_num_nodes_per_elem() + conv.get_node_map(k); // global index 
	      int node_number   = ex.get_connect(gi);             // Global node number (1-based)
	      elem->set_node(k) = mesh.node_ptr((node_number-1)); // Set node number
	                                                          // Subtract 1 since
		                                                  // exodus is internally 1-based
	    }
	}
      
      // running sum of # of elements per block,
      // (should equal total number of elements in the end)
      nelem_last_block += ex.get_num_elem_this_blk();
    }

  assert (static_cast<unsigned int>(nelem_last_block) == mesh.n_elem());
  
  // Read in sideset information -- this is useful for applying boundary conditions
  {
    ex.read_sideset_info(); // Get basic information about ALL sidesets
    int offset=0;
    for (int i=0; i<ex.get_num_side_sets(); i++)
      {
	offset += (i > 0 ? ex.get_num_sides_per_set(i-1) : 0); // Compute new offset
	ex.read_sideset (i, offset);
      }
    
    const std::vector<int>& elem_list = ex.get_elem_list();
    const std::vector<int>& side_list = ex.get_side_list();
    const std::vector<int>& id_list   = ex.get_id_list();
    
    for (unsigned int e=0; e<elem_list.size(); e++)
      {
	// Set any relevant node/edge maps for this element
	const ExodusII::Conversion conv =
	  em.assign_conversion(mesh.elem(elem_list[e]-1)->type());
	
	mesh.boundary_info->add_side (elem_list[e]-1,
				      conv.get_side_map(side_list[e]-1),
				      id_list[e]);
      }
  }
    
  ex.close();            // Close the exodus file, if possible

#endif
}
