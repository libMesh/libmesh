// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __unv_io_h__
#define __unv_io_h__


// C++ inludes
#include <vector>
#include <map>
#include <string>

// Local includes
#include "mesh_input.h"
#include "mesh_output.h"

// Forward declarations
class MeshBase;
class MeshData;

/**
 * The \p UNVIO class implements the Ideas \p UNV universal
 * file format.  This class enables both reading and writing
 * \p UNV files.
 */

// ------------------------------------------------------------
// UNVIO class definition
class UNVIO : public MeshInput<MeshBase>,
	      public MeshOutput<MeshBase>
{

 public:
  
  /**
   * Constructor.  Takes a writeable reference to a mesh object.
   * This is the constructor required to read a mesh.
   */
  UNVIO (MeshBase& mesh, MeshData& mesh_data);

  /**
   * Constructor.  Takes a reference to a constant mesh object.
   * This constructor will only allow us to write the mesh.
   */
  UNVIO (const MeshBase& mesh, MeshData& mesh_data);
  
  /**
   * Destructor.
   */
  virtual ~UNVIO ();
  
  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string& );
  
  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string& );

  /**
   * Set the flag indicationg if we should be verbose.
   */
  bool & verbose ();

  
 private:
  

  /**
   * The actual implementation of the read function.
   * The public read interface simply decides which
   * type of stream to pass the implementation.
   */
  void read_implementation (std::istream& in_stream);

  /**
   * The actual implementation of the write function.
   * The public write interface simply decides which
   * type of stream to pass the implementation.
   */
  void write_implementation (std::ostream& out_stream);
  
  /**
   * Clears the data structures to a pristine
   * state.
   */
  void clear();


  //-------------------------------------------------------------
  // read support methods
  /**
   * When reading, counting the nodes first
   * helps pre-allocation.  Also determine
   * whether we need to convert from "D" to "e".
   */
  void count_nodes (std::istream& in_file);
  
  /**
   * Method reads nodes from \p in_file and stores them in
   * vector<Node*> \p nodes in the order they come in.
   * The original node labels are being stored in
   * \p _assign_nodes in order to assign the elements to
   * the correct nodes later.  In addition, provided it is 
   * active, the \p MeshData gets to know the node id from
   * the Universal file, too.
   */
  void node_in (std::istream& in_file);

  /**
   * When reading, counting the elements first
   * helps pre-allocation.
   */
  void count_elements (std::istream& in_file);

  /**
   * Method reads elements and stores them in
   * \p std::vector<Elem*> \p _elements in the same order as they
   * come in. Within \p UNVIO, element labels are
   * ignored, but \p MeshData takes care of such things
   * (if active).
   */
  void element_in (std::istream& in_file);

  /**
   * @returns \p false when error occured, \p true otherwise.
   * Adjusts the \p in_stream to the beginning of the
   * dataset \p ds_name.
   */
  bool beginning_of_dataset (std::istream& in_file, 
			     const std::string& ds_name) const;

  /**
   * Method for converting exponential notation
   * from "D" to "e", for example
   * \p 3.141592654D+00 \p --> \p 3.141592654e+00
   * in order to make it readable for C++.
   */
  Real D_to_e (std::string& number) const;


  //-------------------------------------------------------------
  // write support methods
  /**
   * Outputs nodes to the file \p out_file.
   * For this to work, the \p MeshData of the current
   * \p MeshBase has to be active.  Do not use this directly,
   * but through the proper write method.
   */
  void node_out (std::ostream& out_file);

  /**
   * Outputs the element data to the file \p out_file.
   * For this to work, the \p MeshData of the current
   * \p Mesh has to be active. Do not use this directly,
   * but through the proper write method.
   */
  void element_out (std::ostream& out_file);


  //-------------------------------------------------------------
  // local data

  /**
   * should be be verbose?
   */
  bool _verbose;

  /**
   * maps node id's from UNV to internal.  Used when reading.
   */
  std::vector<unsigned int> _assign_nodes;

  /**
   * stores positions of relevant datasets in the file, should
   * help to re-read the data faster.  Used when reading.
   */
  std::map<std::string,std::streampos> _ds_position;

  /**
   * total number of nodes, determined through \p count_nodes().
   * Primarily used when reading.
   */
  unsigned int _n_nodes;

  /**
   * total number of elements, determined through 
   * \p count_elements().  Primarily used when reading.
   */
  unsigned int _n_elements;

  /**
   * label for the node dataset
   */
  static const std::string _label_dataset_nodes;

  /**
   * label for the element dataset
   */
  static const std::string _label_dataset_elements;

  /**
   * whether we need to convert notation of exponentials.
   * Used when reading.
   */
  bool _need_D_to_e;

  /**
   * A pointer to the MeshData object you would like to use.
   * with this UNVIO object.  Can be NULL.
   */
  MeshData& _mesh_data;

};



// ------------------------------------------------------------
// MeshIO inline members
inline
UNVIO::UNVIO (MeshBase& mesh, MeshData& mesh_data) :
  MeshInput<MeshBase> (mesh),
  MeshOutput<MeshBase>(mesh),
  _verbose (false),
  _mesh_data (mesh_data)
{
}



inline
UNVIO::UNVIO (const MeshBase& mesh, MeshData& mesh_data) :
  MeshOutput<MeshBase> (mesh),
  _verbose (false),
  _mesh_data (mesh_data)
{
}



inline
UNVIO::~UNVIO ()
{
  this->clear ();
}



inline
bool & UNVIO::verbose ()
{
  return _verbose;
}



inline
bool UNVIO::beginning_of_dataset (std::istream& in_file, 
				  const std::string& ds_name) const
{
  assert (in_file.good());
  assert (!ds_name.empty());

  std::string olds, news;

  while (true)
    {
      in_file >> olds >> news;

      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while( ((olds != "-1") || (news == "-1") ) && !in_file.eof() )
	{	  
	  olds = news;
	  in_file >> news;
	}

      if (in_file.eof())
	return false;
      
      if (news == ds_name)
	return true;
    }

  // should never end up here
  libmesh_error();
  return false;
}



inline
Real UNVIO::D_to_e (std::string& number) const
{
  /* find "D" in string, start looking at 
   * 6th element, to improve speed.
   * We dont expect a "D" earlier
   */

#ifdef __HP_aCC
  // Use an int instead of an unsigned int,
  // otherwise HP aCC may crash!
  const int position = number.find("D",6);
#else
  const std::string::size_type position = number.find("D",6);
#endif

  assert (position != std::string::npos);
  number.replace(position, 1, "e"); 

  return std::atof (number.c_str());
}



#endif // #define __unv_io_h__
