// $Id: data_map.h,v 1.1.2.1 2003-05-09 12:55:39 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __data_map_h__
#define __data_map_h__


// C++ includes
#include <map>
#include <string>
#include <typeinfo>

// Local includes
#include "mesh_common.h"




//----------------------------------------
// DataObjectBase
  
/**
 * The base class for a generic data object
 */
class DataObjectBase
{
public:

  /**
   * Destructor. 
   */
  virtual ~DataObjectBase () {}
  
  /**
   * Returns the name of the class used in construction.
   */
  const std::string & class_name () const { return _class_name; }

  
protected:

  
  /**
   * Constructor.  Takes the name of the class pointed to.
   */
  DataObjectBase (const std::string& class_name) : _class_name(class_name) {}
  
  /**
   * Build an object of the requested type.
   * Note that this COPIES the object \p data,
   * so you should pass POINTERS to large objects,
   * not the objects themselves!
   */
  template <typename T>
  static DataObjectBase* build (const T& data);
  
  /**
   * Build an object of the requested type.
   */
  template <typename T>
  static DataObjectBase* build (const T* data_ptr);


private:
  
  /**
   * The class name.
   */
  const std::string _class_name;
    
  /**
   * Friends 
   */
  friend class DataMap;
};





//-----------------------------------------------------------
// DataObject

/**
 * Derived class.  Implements storage of arbtrary data types.
 * This class copies the data, so it is safe to use with literals.
 */
template <typename T>
class DataObject : public DataObjectBase
{
private:
  
  /**
   * Construct from an object.  NOTE THAT THIS COPIES THE OBJECT!
   */
  DataObject (const T& data) :
    DataObjectBase(typeid(T).name()),
    _data (data)
  {}
  
  /**
   * @returns a constant reference to the data.
   */
  const T& get_data () const { return _data; }

  /**
   * The actual data
   */
  const T _data;

  /**
   * Friends 
   */
  friend class DataObjectBase;
  friend class DataMap;
};



/**
 * Derived class.  This class is multiply inherited based on input type,
 * so it is only safe to use with objects.  Furthermore, this class
 * must be able to inherit from \p T, and T needs an accessible default
 * constructor.
 */
template <class T>
class DataObjectPtr : public DataObjectBase, public T
{ 
private:

  /**
   * Constructor.
   */
  DataObjectPtr (const T* data_ptr) :
    DataObjectBase (typeid(T).name()),
    _data_ptr (data_ptr)
  {}

  /**
   * Pointer to the object.
   */
  const T* _data_ptr;

  /**
   * Friends 
   */
  friend class DataObjectBase;
};



//-------------------------------------------------------------
// DataMap

/**
 * Class that maps a string to an arbitrary data type with
 * type safety.
 */
class DataMap
{
  
public:
  
  /**
   * Constructor.
   */
  DataMap () {}

  /**
   * Destructor.
   */
  ~DataMap () { this->clear(); } 

  /**
   * Clears the object & returns to a pristine state.
   */
  void clear ();

  /**
   * Add the data \p data with the associated \p name;
   */
  template <typename T>
  void add_data (const std::string& name, const T& data);
  
  /**
   * Add the data \p data with the associated \p name;
   */
  template <typename T>
  void add_data_ptr (const std::string& name, const T* data_ptr);
  
  /**
   * Get the data specified by \p name.  The return type is
   * specified by the template parameter.
   */
  template <typename T>
  const T & get_data (const std::string& name) const;
  
  /**
   * Get the data specified by \p name.  The return type is
   * specified by the template parameter.
   */
  template <typename T>
  const T * get_data_ptr (const std::string& name) const;
   
private:

  /**
   * Data container.
   */
  std::map<std::string, DataObjectBase*> _data_map;
};



//----------------------------------------------------------------
// DataMap inline members
template <typename T>
inline
DataObjectBase* DataObjectBase::build (const T& data)
{
  return new DataObject<T> (data);
}



template <typename T>
inline
DataObjectBase* DataObjectBase::build (const T* data_ptr)
{
  assert (data_ptr != NULL);
  return new DataObjectPtr<T> (data_ptr);
}



template <typename T>
inline
void DataMap::add_data (const std::string& name, const T& data)
{
  // Check to see if we have somethint with this name already
  std::map<std::string, DataObjectBase*>::iterator
    pos = _data_map.find (name);


  // Make sure it wasn't already there
  if (pos != _data_map.end())
    {
      std::cerr << "ERROR:  Data already exists for \""
		<< name
		<< "\", aborting!"
		<< std::endl;

      error();
    }

  
  // Didn't find it
  else
    {
      _data_map[name] = DataObjectBase::build (data);
    }
}



template <typename T>
inline
void DataMap::add_data_ptr (const std::string& name, const T* data_ptr)
{
  // Check to see if we have somethint with this name already
  std::map<std::string, DataObjectBase*>::iterator
    pos = _data_map.find (name);


  // Make sure it wasn't already there
  if (pos != _data_map.end())
    {
      std::cerr << "ERROR:  Data already exists for \""
		<< name
		<< "\", aborting!"
		<< std::endl;

      error();
    }
  
  // Didn't find it
  else
    {
      _data_map[name] = DataObjectBase::build (data_ptr);
    }
}



template <typename T>
inline
const T& DataMap::get_data (const std::string& name) const
{
  std::map<std::string, DataObjectBase*>::const_iterator
    pos = _data_map.find(name);
  
  if (pos == _data_map.end())
    {
      std::cerr << "ERROR:  No data associated with the name \""
		<< name << "\"!"
		<< std::endl;

      error();
    }
  
  assert (pos->second != NULL);

  
  // the object must be there.  Look above
  DataObject<T>* obj = dynamic_cast<DataObject<T>* >(pos->second);

  // Check for failed cast
  if (obj == NULL)
    {      
      std::cerr << "ERROR: Cast failed, invalid type requested!"
		<< std::endl
		<< "The known type is:     \""         << pos->second->class_name() << "\""
		<< std::endl
		<< "The requested type is: \""         << typeid(T).name() << "\""
		<< std::endl;
      error();
    }


  return obj->get_data ();
}



template <typename T>
inline
const T* DataMap::get_data_ptr (const std::string& name) const
{
  std::map<std::string, DataObjectBase*>::const_iterator
    pos = _data_map.find(name);
  
  if (pos == _data_map.end())
    {
      std::cerr << "ERROR:  No data associated with the name \""
		<< name << "\"!"
		<< std::endl;

      error();
    }


  // Can't happen...
  assert (pos->second != NULL);
  
  // the object must be there.  Look above
  T* obj = dynamic_cast<T*>(pos->second);

  // Check for failed cast
  if (obj == NULL)
    {      
      std::cerr << "ERROR: Cast failed, invalid type requested!"
		<< std::endl
		<< "The known type is:     \""  << pos->second->class_name() << "\""
		<< std::endl
		<< "The requested type is: \""  << typeid(T).name() << "\""
		<< std::endl
		<< "Cannot cast from a \""      << pos->second->class_name()
		<< "\" to a \"" << typeid(T).name() << "\", aborting!"
		<< std::endl;
      error();
    }

  return obj;
}


#endif // #define __data_map_h__
