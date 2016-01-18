// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_MESH_DATA_H
#define LIBMESH_MESH_DATA_H

// Local Includes
#include "libmesh/libmesh.h"
#include "libmesh/enum_xdr_mode.h"

// C++ includes
#include <cstddef>
#include <map>
#include <vector>

namespace libMesh
{

// Forward Declarations
class Node;
class Elem;
class MeshBase;
class MeshDataUnvHeader;

/**
 * The \p MeshData class handles actual data and the corresponding
 * I/O on entities (nodes, elements) of meshes.
 * The \p MeshData can be used when dealing with files
 * that contain nodal or element-oriented data, numbered in the same
 * format as a corresponding mesh file (when activated) or with
 * the \p libMesh element and node indices (when in compatibility mode).
 * To use \p MeshData, it has to be either activated or the compatibility
 * mode has to be enabled.
 *
 * \author Daniel Dreyer
 * \date 2003
 */
class MeshData
{
public:
  //----------------------------------------------------------
  // convenient typedefs
  /**
   * A const iterator over the nodal data entries of
   * \p MeshData.  Use this when a loop over all \p Node *
   * in the \p MeshData is wanted.  Note that only const versions
   * are provided.  Also these iterators should @e not be
   * confused with the \p node_iterators provided
   * for the \p Mesh classes!
   */
  typedef std::map<const Node *, std::vector<Number> >::const_iterator const_node_data_iterator;

  /**
   * A const iterator over the element-associated data entries of
   * \p MeshData.  Use this when a loop over all \p Node *
   * in the \p MeshData is wanted.  Note that only const versions
   * are provided.  Also these iterators should @e not be
   * confused with the \p node_iterators provided
   * for the \p Mesh classes!
   */
  typedef std::map<const Elem *, std::vector<Number> >::const_iterator const_elem_data_iterator;


  //----------------------------------------------------------
  /**
   * Default Constructor.  Takes const reference
   * to the mesh it belongs to.
   */
  MeshData (const MeshBase & m);

  /**
   * Destructor.
   */
  ~MeshData ();

  /**
   * When \p MeshData should be used, it has to be activated
   * first, @e prior to reading in a mesh with the \p Mesh::read()
   * methods. This will ensure that element and node ids
   * given in the mesh file, i.e. the foreign node and element
   * ids, are stored in the corresponding id maps.
   * Optionally takes a string that should help the user
   * in identifying the data later on.
   */
  void activate (const std::string & descriptor="");

  /**
   * When the \p MeshData should be used, but was @e not activated
   * prior to reading in a mesh, then the compatibility mode enables
   * to still use this object as if the \p MeshData was active.
   * The foreign node and element ids are simply assigned the
   * indices used in \p libMesh.  Note that the compatibility mode
   * should be used with caution, since the node and element
   * indices in \p libMesh may be renumbered any time.  This
   * \p MeshData always employs the current node and element ids,
   * it does @e not create an image of ids when compatibility
   * mode was activated.
   */
  void enable_compatibility_mode (const std::string & descriptor="");

  /**
   * Clears the data fields, but leaves the id maps
   * untouched.  Useful for clearing data for a new
   * data file.  Use \p slim() to delete the maps.
   */
  void clear ();

  /**
   * Once the data is properly read from file, the id
   * maps can safely be cleared.  However, if this object
   * should remain able to @e write nodal or element oriented
   * data to file, this method should better @e not be used.
   * Use the appropriate \p bool to select the id map that
   * should be cleared.  By default, both id maps are deleted.
   */
  void slim (const bool node_id_map = true,
             const bool elem_id_map = true);

  /**
   * Translates the @e nodal data contained in this object
   * to \p data_values and \p data_names.  These two
   * vectors are particularly suitable for use with
   * the \p MeshBase::write method that takes nodal
   * data.  E.g., the export method may be used for
   * inspecting boundary conditions.  A reference
   * to the mesh for which the data should be written
   * has to be provided.  Note that this mesh @e has
   * to contain the nodes for which this \p MeshData
   * holds data.  I.e., \p out_mesh may only refer to
   * the \p MeshBase itself (that this \p MeshData belongs
   * to), or its \p BoundaryMesh, cf. \p Mesh.
   */
  void translate (const MeshBase & out_mesh,
                  std::vector<Number> & data_values,
                  std::vector<std::string> & data_names) const;

  /**
   * Read mesh data from file named \p name.
   * Guess format from the file extension.  Note that
   * prior to this you have to at least either
   * \p close_node_map() or \p close_elem_map().
   */
  void read (const std::string & name);

  /**
   * Write mesh data to file named \p name.
   * Guess format from the file extension.
   */
  void write (const std::string & name);

  /**
   * @returns a string containing relevant information
   * about the mesh.
   */
  std::string get_info () const;

  /**
   * Prints relevant information about the mesh.
   */
  void print_info (std::ostream & os=libMesh::out) const;

  /**
   * Same as above, but allows you to use the stream syntax.
   */
  friend std::ostream & operator << (std::ostream & os, const MeshData & m);


  //----------------------------------------------------------
  // Node-associated data
  /**
   * @returns the \f$ i^{th} \f$ value (defaults to 0) associated
   * with node \p node.  Returns \p libMesh::zero when there
   * is no such \p node in the map.
   */
  Number operator() (const Node * node,
                     const unsigned int i=0) const;

  /**
   * @returns \p true when the node \p node has data,
   * \p false otherwise.
   */
  bool has_data (const Node * node) const;

  /**
   * @returns a const reference to the values associated with
   * the node \p node.  @e Beware: this method will crash
   * when there is no data associated with the node \p node!
   * Check existence through \p has_data() first.
   */
  const std::vector<Number> & get_data (const Node * node) const;

  /**
   * Sets all the data values associated with
   * the node \p node, overwriting any existing vector
   */
  void set_data (const Node * node, const std::vector<Number> & val);

  /**
   * @returns the number of \p Number -type data
   * (i.e., the size of the \p std::vector<Number>
   * returned through the \p operator() methods)
   * associated with a node.  Returns 0 when no
   * nodal data exists.
   */
  unsigned int n_val_per_node () const;

  /**
   * @returns the number of nodes for which this
   * \p MeshData has data stored.
   */
  dof_id_type n_node_data () const;

  /**
   * Returns the \p MeshData::const_node_data_iterator which points
   * to the beginning of the \p Node * data containers
   * used here.
   */
  const_node_data_iterator node_data_begin () const;

  /**
   * Returns the \p MeshData::const_node_data_iterator which points
   * to the end of the \p Node * data containers used here.
   */
  const_node_data_iterator node_data_end () const;

  /**
   * For the desperate user, nodal boundary conditions
   * may be inserted directly through the map \p nd.
   * It is mandatory that there does not yet exist any
   * other node data in this object, that the id maps
   * are closed, that the size of the std::vector's of
   * each map have identical length and that the Node *
   * point to nodes of the associated mesh.
   * Note that this method takes a non-const reference
   * and essentially clears the passed-in data.
   * If \p close_elem_data is \p true (default), then
   * this \p MeshData is ready for use: write to file,
   * use the operator() methods etc. If \p false, the
   * user @e has to add element-associated data, too.
   */
  void insert_node_data (std::map<const Node *,
                         std::vector<Number> > & nd,
                         const bool close_elem_data = true);


  //----------------------------------------------------------
  // Element-associated data
  /**
   * @returns the \f$ i^{th} \f$ value (defaults to 0) associated
   * with element \p elem.  Returns \p libMesh::zero when there
   * is no data for \p elem in the map.
   */
  Number operator() (const Elem * elem,
                     const unsigned int i=0) const;

  /**
   * @returns \p true when the element \p elem has data,
   * \p false otherwise.
   */
  bool has_data (const Elem * elem) const;

  /**
   * @returns a const reference to the values associated with
   * the element \p elem.  @e Beware: this method will crash
   * when there is no data associated with the element \p elem!
   * Check existence through \p has_data() first.
   */
  const std::vector<Number> & get_data (const Elem * elem) const;

  /**
   * Sets all the data values associated with
   * the element \p elem, overwriting any existing vector
   */
  void set_data (const Elem * elem, const std::vector<Number> & val);

  /**
   * @returns the number of \p Number -type data
   * (i.e., the size of the \p std::vector<Number>
   * returned through the \p operator() methods)
   * associated with an element.  Returns 0 when
   * there is no element-associated data.
   */
  unsigned int n_val_per_elem () const;

  /**
   * @returns the number of elements for which this
   * \p MeshData has data stored.
   */
  dof_id_type n_elem_data () const;

  /**
   * Returns a \p MeshData::const_elem_data_iterators which points
   * to the beginning of the \p Elem * data containers
   * used here.
   */
  const_elem_data_iterator elem_data_begin () const;

  /**
   * Returns a \p MeshData::const_elem_data_iterators which points
   * to the end of the \p Elem * data containers used here.
   */
  const_elem_data_iterator elem_data_end () const;

  /**
   * For the desperate user, element-associated boundary
   * conditions may be inserted directly through the
   * map \p ed.  Similar to the version for nodal data,
   * it is imperative that the local \p _elem_data is empty,
   * that the id maps are closed, that the size of the
   * \p std::vector's of each map have identical length
   * and that the \p Elem points to elements of the
   * associated mesh.
   * Note that this method takes a non-const reference
   * and essentially clears the passed-in data.
   * If \p close_node_data is \p true (default), then
   * this \p MeshData is ready for use: write to file,
   * use the operator() methods etc. If \p false, the
   * user @e has to add nodal data, too.
   */
  void insert_elem_data (std::map<const Elem *, std::vector<Number> > & ed,
                         const bool close_node_data = true);


  //----------------------------------------------------------
  /**
   * @returns \p true when this object is active and working.
   * Use \p activate() to bring this object alive.
   */
  bool active () const;
  /**
   * @returns \p true when this object is in compatibility
   * mode.  See \p enable_compatibility_mode() for details.
   */
  bool compatibility_mode () const;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use for @e element associated data, \p false
   * otherwise.
   */
  bool elem_initialized () const;

  /**
   * @returns \p true when this object is properly initialized
   * and ready for use for @e nodal data, \p false otherwise.
   */
  bool node_initialized () const;


  //----------------------------------------------------------
  // Methods for accessing the node and element maps.
  // Heavily used by the \p read() and \p write() methods.
  /**
   * @returns the \p Node * that this foreign id maps to.
   */
  const Node * foreign_id_to_node (const unsigned int fid) const;

  /**
   * @returns the \p Elem * that this foreign id maps to.
   */
  const Elem * foreign_id_to_elem (const unsigned int fid) const;

  /**
   * @returns the foreign id this \p Node * maps to.
   */
  unsigned int node_to_foreign_id (const Node * n) const;

  /**
   * @returns the foreign id this \p Elem * maps to.
   */
  unsigned int elem_to_foreign_id (const Elem * n) const;

  //----------------------------------------------------------
  // Methods for the header information in universal formated
  // datasets.

  /**
   * Read access to the \p MeshDataUnvHeader data structure.
   */
  const MeshDataUnvHeader & get_unv_header() const;

  /**
   * Set the \p MeshDataUnvHeader data structure that will be
   * used for output.
   */
  void set_unv_header(MeshDataUnvHeader * unv_header);


  /**
   * Assign to \p this the data from the other \p MeshData.
   * Used by \p BoundaryInfo when copying the \p MeshData
   * from the \p d dimensional mesh to the \p d-1 dimensional mesh
   * (the boundary mesh).
   */
  void assign (const MeshData & omd);


  //----------------------------------------------------------
  // Methods used by mesh importes to communicate node/element
  // labels to this \p MeshData

  /**
   * In general, \p MeshData gathers nodal data
   * from a file, but it needs to relate this data
   * with the \p Node * of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_node_id (const Node * node,
                            const unsigned int foreign_node_id);

  /**
   * In general, \p MeshData gathers element-associated
   * data from file, but it needs to relate this data
   * with the \p Elem * of the current mesh.  Mesh
   * importers simply use this method to add such
   * a map.
   */
  void add_foreign_elem_id (const Elem * elem,
                            const unsigned int foreign_elem_id);

  /**
   * Signal to this object that the mesh importer finished
   * adding node and element foreign-id maps.
   */
  void close_foreign_id_maps ();


protected:


  //----------------------------------------------------------
  // read/write Methods
  /**
   * Read nodal/element oriented data in TetGen format.
   */
  void read_tetgen (const std::string & name);

  /**
   * Read nodal/element oriented data in UNV format,
   * either from an ASCII file or from a gzip'ed ASCII
   * file, using the C++ wrapper \p gzstream to \p zlib.h.
   */
  void read_unv (const std::string & file_name);

  /**
   * Actual implementation of reading nodal/element
   * oriented data in UNV format.  This has to be
   * decoupled from \p read_unv() in order to allow
   * reading both \p .unv and \p .unv.gz files.
   */
  void read_unv_implementation (std::istream & in_file);

  /**
   * Write nodal/element oriented data in UNV format,
   * either to an ASCII file or to a gzip'ed ASCII
   * file, using the C++ wrapper \p gzstream to \p zlib.h.
   */
  void write_unv (const std::string & file_name);

  /**
   * Actual implementation of writing nodal/element
   * oriented data in UNV format.  This has to be
   * decoupled from \p write_unv() in order to allow
   * writing both \p .unv and \p .unv.gz files.
   */
  void write_unv_implementation (std::ostream & out_file);


  /**
   * Read nodal/element oriented data using the
   * \p Xdr class that enables both ASCII and
   * binary format through the same interface.
   * By default uses ASCII format, but may easily
   * be changed setting \p mode to \p DECODE.
   */
  void read_xdr (const std::string & name,
                 const XdrMODE mode = READ);

  /**
   * Write nodal data in format comparable to
   * the XDR format already known from \p Mesh.
   * By default uses ASCII format, but may easily
   * be changed setting \p mode to \p ENCODE.
   */
  void write_xdr (const std::string & name,
                  const XdrMODE mode = WRITE);


  /**
   * The mesh this object belongs to
   */
  const MeshBase & _mesh;

  /**
   * Some name the user gave to the data when this
   * object got activated
   */
  std::string _data_descriptor;


  //--------------------------------------------------
  // node associated data & maps
  /**
   * The map containing pointers to nodes in the mesh
   * and the corresponding data.
   */
  std::map<const Node *,
           std::vector<Number> > _node_data;

  /**
   * Maps node pointers to node numbers in the @e foreign
   * format.
   */
  std::map<const Node *,
           unsigned int> _node_id;

  /**
   * Maps @e foreign node ids to node pointers of the
   * current mesh.
   */
  std::map<unsigned int,
           const Node *> _id_node;



  //--------------------------------------------------
  // element associated data & maps
  /**
   * Maps element pointers to the element-associated data
   */
  std::map<const Elem *,
           std::vector<Number> > _elem_data;

  /**
   * Maps element pointers to element labels in the @e foreign
   * format.
   */
  std::map<const Elem *,
           unsigned int> _elem_id;
  /**
   * Maps @e foreign element labels to element pointers of the
   * current mesh.
   */
  std::map<unsigned int,
           const Elem *> _id_elem;



  //--------------------------------------------------------
  /**
   * \p true when the mesh importer finished adding
   * node-foreign-id maps, and the node-foreign-id maps
   * exist.  Note that these maps may be deleted through
   * \p slim() to save memory.  Then the data is
   * still accessible through the \p Node * or \p Elem *,
   * but the foreign id's are lost.
   */
  bool _node_id_map_closed;

  /**
   * \p true when the nodal data are properly initialized,
   * false otherwise.
   */
  bool _node_data_closed;


  //--------------------------------------------------------
  /**
   * \p true when the mesh importer finished adding
   * element-id maps, and the element-id maps exist.
   * Note that these maps may be deleted through
   * \p slim() to save memory.  Then the data is
   * still accessible through the \p Elem *,
   * but the foreign element id's are lost.
   */
  bool _elem_id_map_closed;

  /**
   * \p true when the element based data are properly initialized,
   * false otherwise.
   */
  bool _elem_data_closed;


  //--------------------------------------------------------
  /**
   * \p true when this object is set active (to gather data
   * during mesh import).
   */
  bool _active;

  /**
   * \p true when this object is in compatibility mode
   * (use libMesh's node and element numbers as fake
   * foreign id's)
   */
  bool _compatibility_mode;

  /**
   * The header information of universal files.
   */
  MeshDataUnvHeader * _unv_header;

  /**
   * Make the \p MeshDataUnvHeader class a friend.
   */
  friend class MeshDataUnvHeader;

};



/**
 * Class \p MeshDataUnvHeader handles the data specified at
 * the @e beginning of a dataset 2414 in a universal file.
 * This header is structured in records 1 to 13.  A typical
 * header is described here.  The text after the # are comments
 * and are @e not part of such a dataset.  The text in brackets
 * after the # are the corresponding class members names.
 *
 * \verbatim
 *
 * -1                                                                           # beginning of dataset
 * 2414                                                                         # type of dataset: data at mesh entities
 * 1                                                                            # R.  1: unique number of dataset (dataset_label)
 * STRUCTURAL MODE     1                                                        # R.  2: text describing content (dataset_name)
 * 1                                                                            # R.  3: data belongs to: nodes, elements,...
 * #        (dataset_location)
 * Default Model                                                                # R.  4: user-specified text (id_lines_1_to_5[0])
 * I-DEAS Master Series                                                         # R.  5: user-specified text (id_lines_1_to_5[1])
 * 18-AUG-2003 20:00:12    HPUX11_64     MAR2003                                # R.  6: user-specified text (id_lines_1_to_5[2])
 * MODE   1 FREQUENCY       501.25 Hz                                           # R.  7: user-specified text (id_lines_1_to_5[3])
 * STRUCTURAL MODE     1                                                        # R.  8: user-specified text (id_lines_1_to_5[4])
 * 0         2         3         8         2         6                          # R.  9: (model_type) (analysis_type)
 * #        (data_characteristic) (result_type)
 * #        (data_type) (nvaldc)
 * 0         0         0         0         0         1         0         0      # R. 10: analysis-specific data (record_10)
 * 0         0                                                                  # R. 11: analysis-specific data (record_11)
 * 0.00000E+00  0.50125E+03  0.99192E+07  0.10000E+01  0.00000E+00  0.00000E+00 # R. 12: analysis-specific data (record_12)
 * 0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00 # R. 13: analysis-specific data (record_13)
 * \endverbatim
 *
 * For more details we refer to the general description of the I-DEAS
 * universal file format.
 *
 * An instance of this class may be attached to the \p MeshData
 * of some mesh.  Then the \p read() and \p write() methods
 * of \p MeshData use this \p MeshDataUnvHeader instead of
 * some empty default.  Also files that contain multiple
 * datasets of type \p 2414 may be handled through the
 * \p which_dataset() method.
 * Note that an instance of this class has to be attached
 * to the \p MeshData @e prior to using the \p read() or
 * \p write() methods of the \p MeshData.
 */
class MeshDataUnvHeader
{
public:

  /**
   * Default Constructor.  Initializes the respective
   * data.
   */
  MeshDataUnvHeader ();

  /**
   * Destructor.
   */
  ~MeshDataUnvHeader ();

  /**
   * Universal files may contain multiple data sets of type
   * \p 2414.  These sets are identified through their
   * labels (not to be confused with the dataset label \p 2414!).
   * The user may provide a label of the dataset that she
   * wants.  Then the file is scanned for this dataset, and
   * datasets with a different label are skipped.
   *
   * When this method is @e not called, then simply the first
   * dataset in the file is used.  Note that for this method
   * to have any effect, this method has to be called prior to
   * using the \p MeshData::read() or \p MeshData::write()
   * methods.
   */
  void which_dataset (const unsigned int ds_label);

  /**
   * Assignment operator.  Simply assigns all values from
   \p omduh to \p this.
  */
  void operator = (const MeshDataUnvHeader & omduh);

  /**
   * @returns \p true when \p this and \p omduh are equal,
   * \p false otherwise.
   */
  bool operator == (const MeshDataUnvHeader & omduh) const;

  /**
   * Record 1.  User specified analysis dataset label.
   */
  unsigned int dataset_label;

  /**
   * Record 2. User specified analysis dataset name.
   */
  std::string dataset_name;

  /**
   * Record 3. The dataset location (e.g. data at nodes,
   * data on elements, etc.).
   */
  unsigned int dataset_location;

  /**
   * Record 4 trough 8 are ID lines.
   */
  std::vector<std::string> id_lines_1_to_5;

  /**
   * Record 9, first part. This record contains data specifying
   * the model type (e.g. unknown, structural, etc.),
   * the analysis type (e.g. unknown, static, transient,
   * normal mode, etc.),
   * the data characteristics (such as scalar, 3 dof global
   * translation vector, etc.),
   * the result type (e.g. stress, strain, velocity, etc.).
   */
  unsigned int model_type,
    analysis_type,
    data_characteristic,
    result_type;

  /**
   * Record 9, second part. See first part, then we have:
   * the data type (currently supported: 2,4 for \p Real,
   * and 5,6 for \p Complex. other possibilities: e.g. integer),
   */
  unsigned int data_type;

  /**
   * Record 9, third and last part. See first and second part,
   * then we have: the number of data values for the mesh data.
   */
  unsigned int nvaldc;

  /**
   * Record 10 and 11 are analysis specific data of
   * type integer.
   */
  std::vector<int> record_10,
    record_11;

  /**
   * Record 12 and 13 are analysis specific data of
   * type Real.
   */
  std::vector<Real> record_12,
    record_13;


protected:

  /**
   * @returns \p true when this dataset is the one
   * that the user wants, \p false otherwise.  When
   * no desired dataset is given, always returns
   * \p true.  Aside from this return value, this method
   * also reads the header information from the
   * stream \p in_file.
   */
  bool read (std::istream & in_file);

  /**
   * Write the header information to the stream \p out_file.
   */
  void write (std::ostream & out_file);


private:

  /**
   * the desired dataset label.  defaults to -1
   * if not given
   */
  unsigned int _desired_dataset_label;

  /**
   * @returns \p true when the string \p number
   * has a 'D' that needs to be replaced by 'e',
   * \p false otherwise.  Also actually replaces
   * the 'D' by an 'e'.
   */
  static bool need_D_to_e (std::string & number);

  /**
   * Make the \p MeshData class a friend.
   */
  friend class MeshData;

};



// ------------------------------------------------------------
// MeshData inline methods

//-------------------------------------------------------------
// element data inline methods
inline
Number MeshData::operator() (const Node * node,
                             const unsigned int i) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_node_data_closed);

  std::map<const Node *,
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

  if (pos == _node_data.end())
    return libMesh::zero;

  // we only get here when pos != _node_data.end()
  libmesh_assert_less (i, pos->second.size());
  return pos->second[i];
}



inline
bool MeshData::has_data (const Node * node) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_node_data_closed);

  std::map<const Node *,
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

  return (pos != _node_data.end());
}



inline
const std::vector<Number> & MeshData::get_data (const Node * node) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_node_data_closed);

  std::map<const Node *,
           std::vector<Number> >::const_iterator pos = _node_data.find(node);

#ifdef DEBUG
  if (pos == _node_data.end())
    libmesh_error_msg("ERROR: No data for this node.  Use has_data() first!");
#endif

  return pos->second;
}



inline
void MeshData::set_data (const Node * node,
                         const std::vector<Number> & val)
{
  this->_node_data[node] = val;
}



inline
MeshData::const_node_data_iterator MeshData::node_data_begin () const
{
  return _node_data.begin();
}



inline
MeshData::const_node_data_iterator MeshData::node_data_end () const
{
  return _node_data.end();
}



//-------------------------------------------------------------
// element data inline methods
inline
Number MeshData::operator() (const Elem * elem,
                             const unsigned int i) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_elem_data_closed);

  std::map<const Elem *,
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  if (pos == _elem_data.end())
    return libMesh::zero;

  // we only get here when pos != _elem_data.end()
  libmesh_assert_less (i, pos->second.size());
  return pos->second[i];
}



inline
bool MeshData::has_data (const Elem * elem) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_elem_data_closed);

  std::map<const Elem *,
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

  return (pos != _elem_data.end());
}



inline
const std::vector<Number> & MeshData::get_data (const Elem * elem) const
{
  libmesh_assert (_active || _compatibility_mode);
  libmesh_assert (_elem_data_closed);

  std::map<const Elem *,
           std::vector<Number> >::const_iterator pos = _elem_data.find(elem);

#ifdef DEBUG
  if (pos == _elem_data.end())
    libmesh_error_msg("ERROR: No data for this element.  Use has_data() first!");
#endif

  return pos->second;
}



inline
void MeshData::set_data (const Elem * elem,
                         const std::vector<Number> & val)
{
  this->_elem_data[elem] = val;
}



inline
MeshData::const_elem_data_iterator MeshData::elem_data_begin () const
{
  return _elem_data.begin();
}



inline
MeshData::const_elem_data_iterator MeshData::elem_data_end () const
{
  return _elem_data.end();
}



//-------------------------------------------------------------
// other inline methods
inline
bool MeshData::active() const
{
  return _active;
}



inline
bool MeshData::compatibility_mode() const
{
  return _compatibility_mode;
}



inline
bool MeshData::elem_initialized() const
{
  return (_active && _elem_data_closed);
}



inline
bool MeshData::node_initialized() const
{
  return (_active && _node_data_closed);
}



inline
void MeshData::add_foreign_node_id (const Node * node,
                                    const unsigned int foreign_node_id)
{
  if (_active)
    {
      libmesh_assert (!_node_id_map_closed);
      libmesh_assert(node);
      libmesh_assert (_node_id.find(node) == _node_id.end());
      libmesh_assert (_id_node.find(foreign_node_id) == _id_node.end());

      /*
       * _always_ insert in _id_node and _node_id.  If we would
       * use the mesh.node(unsigned int) method or the node.id()
       * to get Node * and unsigned int, respectively, we would not
       * be safe any more when the mesh gets refined or re-numbered
       * within libMesh. And we could get in big trouble that would
       * be hard to find when importing data _after_ having refined...
       */
      _node_id.insert(std::make_pair(node, foreign_node_id));
      _id_node.insert(std::make_pair(foreign_node_id, node));
    }
}



inline
void MeshData::add_foreign_elem_id (const Elem * elem,
                                    const unsigned int foreign_elem_id)
{
  if (_active)
    {
      libmesh_assert (!_elem_id_map_closed);
      libmesh_assert(elem);
      libmesh_assert (_elem_id.find(elem) == _elem_id.end());
      libmesh_assert (_id_elem.find(foreign_elem_id) == _id_elem.end());

      _elem_id.insert(std::make_pair(elem, foreign_elem_id));
      _id_elem.insert(std::make_pair(foreign_elem_id, elem));
    }
}


inline
const MeshDataUnvHeader & MeshData::get_unv_header () const
{
  libmesh_assert(this->_unv_header);
  return *this->_unv_header;
}


inline
void MeshData::set_unv_header (MeshDataUnvHeader * unv_header)
{
  libmesh_assert(unv_header);
  this->_unv_header = unv_header;
}


//-----------------------------------------------------------
// MeshDataUnvHeader inline methods


} // namespace libMesh

#endif // LIBMESH_MESH_DATA_H
