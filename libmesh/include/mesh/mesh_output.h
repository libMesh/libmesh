// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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



#ifndef __mesh_output_h__
#define __mesh_output_h__


// C++ inludes
#include <string>
#include <vector>

// Local includes
#include "libmesh_common.h"
#include "mesh_base.h"

// Forward declares
class EquationSystems;


/**
 * This class defines an abstract interface for \p Mesh output.
 * Specific classes derived from this class actually implement
 * writing various mesh formats.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision$
 */

// ------------------------------------------------------------
// MeshOutput class definition
template <class MT>
class MeshOutput
{
 protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  MeshOutput (const bool is_parallel_format = false);
  
  /**
   * Constructor.  Takes a reference to a constant object.
   * This constructor will only allow us to write the object.
   */
  MeshOutput (const MT&, const bool is_parallel_format = false);

  
 public:

  /**
   * Destructor.
   */
  virtual ~MeshOutput ();
  
  /**
   * This method implements writing a mesh to a specified file.
   */
  virtual void write (const std::string&) = 0;

  /**
   * This method implements writing a mesh with data to a specified file
   * where the data is taken from the \p EquationSystems object.
   */
  virtual void write_equation_systems (const std::string&,
				       const EquationSystems&);

  /**
   * This method implements writing a mesh with nodal data to a
   * specified file where the nodal data and variable names are provided.
   */
  virtual void write_nodal_data (const std::string&,
				 const std::vector<Number>&,
				 const std::vector<std::string>&)
  { libmesh_error(); }

  
 protected:


  /**
   * Returns the object as a read-only reference.
   */
  const MT& mesh() const;

  
 private:
  

  /**
   *  A pointer to a constant object.
   * This allows us to write the object to file.
   */
  const MT* const _obj;

  /**
   * Flag specifying whether this format is parallel-capable.
   * If this is false (default) I/O is only permitted when the mesh
   * has been serialized.
   */
  const bool _is_parallel_format;

  /**
   * A helper function which allows us to fill temporary
   * name and solution vectors with an EquationSystems object
   */
  void _build_variable_names_and_solution_vector(const EquationSystems& es,
						 std::vector<Number>& soln,
						 std::vector<std::string>& names);
};



// ------------------------------------------------------------
// MeshOutput inline members
template <class MT>
inline
MeshOutput<MT>::MeshOutput (const bool is_parallel_format) :
  _obj(NULL),
  _is_parallel_format(is_parallel_format)
{}



template <class MT>
inline
MeshOutput<MT>::MeshOutput (const MT& obj, const bool is_parallel_format) :
  _obj (&obj),
  _is_parallel_format(is_parallel_format)
{
  if (!_is_parallel_format && !this->mesh().is_serial())
    {
      if (libMesh::processor_id() == 0)
	{
          std::cerr << "WARNING:  This I/O operation may only be supported for meshes which have been serialized!"
		    << std::endl;
          here();
        }
//      libmesh_error();
    }
}



template <class MT>
inline
MeshOutput<MT>::~MeshOutput ()
{
}



template <class MT>
inline
void MeshOutput<MT>::write_equation_systems (const std::string& fname,
					     const EquationSystems& es)
{
  // Build the nodal solution values & get the variable
  // names from the EquationSystems object
  std::vector<Number>      soln;
  std::vector<std::string> names;

  this->_build_variable_names_and_solution_vector(es, soln, names);
  //es.build_variable_names  (names);
  //es.build_solution_vector (soln);

  this->write_nodal_data (fname, soln, names);  
}



template <class MT>
inline
const MT& MeshOutput<MT>::mesh () const
{
  libmesh_assert (_obj != NULL);
  return *_obj;
}


#endif // #define __mesh_io_h__
