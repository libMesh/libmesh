// $Id: data_map.h,v 1.3 2003-09-02 18:02:36 benkirk Exp $

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



#ifndef __data_map_h__
#define __data_map_h__


// C++ includes
#include <map>
#include <string>
#include <typeinfo>

// Local includes
#include "mesh_common.h"
#include "data_object.h"




//-----------------------------------------------------------
// DataObjectCpy

/**
 * Derived class.  Implements storage of arbtrary data types.
 * This class copies the data, so it is safe to use with literals.
 * Since it copies the data be careful not to call this with
 * heavy objects.
 */
template <typename T>
class DataObjectCpy : public DataObject
{
private:
  
  /**
   * Construct from an object.  NOTE THAT THIS COPIES THE OBJECT!
   * Tells the \p DataObject class what type \p T is for
   * possible debugging.
   */
  DataObjectCpy (const T& data) :
    _data (data)
  {}
    
  /**
   * Build an object of the requested type.
   * Note that this COPIES the object \p data,
   * so you should pass POINTERS to large objects,
   * not the objects themselves!
   */
  static DataObject* build (const T& data);
  
  /**
   * @returns a constant reference to the data.
   */
  const T& get_data () const { return _data; }

  /**
   * The actual data
   */
  const T _data;

  /**
   * Friends.  These are the only classes that should
   * ever be able to use one of these objects since the
   * class is completely private.
   */
  friend class DataMap;
};



// /**
//  * Derived class.  This class is multiply inherited based on input type,
//  * so it is only safe to use with objects.  Furthermore, this class
//  * must be able to inherit from \p T, and T needs an accessible default
//  * constructor.  By inheriting from \p T a \p DataObjectPtr<T> can
//  * be successfully cast to a \p T or any related class via
//  * \p dynamic_cast<>
//  */
// template <class T>
// class DataObjectPtr : public T, public DataObject
// { 
// private:

//   /**
//    * Constructor.
//    */
//   DataObjectPtr (const T* data_ptr) :
//     _data_ptr (data_ptr)
//   {}
  
//   /**
//    * Build an object of the requested type, taking a pointer
//    * to the type. For this to work properly we must be able
//    * to derive from \p T, and furthermore \p T needs an accessible
//    * default constructor.
//    */
//   static DataObject* build (const T* data_ptr);

//   /**
//    * Pointer to the object.
//    */
//   const T* _data_ptr;

//   /**
//    * Friend. This is the only class that should
//    * ever be able to use one of these objects since
//    * the class is completely private. 
//    */
//   friend class DataMap;
// };



//-------------------------------------------------------------
// DataMap

/**
 * Class that maps a string to an arbitrary data type with
 * type safety.  A \p DataMap can map a name to a copy of
 * data or a name to a pointer to an object.  Which method
 * is right for you depends on the data type.
 *
 * @author Michael Anderson, Bill Barth, Benjamin Kirk,
 *  and John Peterson, 2003. 
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
   * Add the data \p data with the associated \p name
   * when \p data has been derived from a \p DataObject
   */
  void add_data_object (const std::string& name, DataObject& data);
  
  /**
   * Add the data \p data with the associated \p name.
   * This method will create a copy of the data.
   */
  template <typename T>
  void add_data (const std::string& name, const T& data);
  
//   /**
//    * Add the data \p data with the associated \p name.
//    * We must be able to derive from a \p T, and \p T
//    * must have an accessible default constructor.
//    */
//   template <typename T>
//   void add_data_ptr (const std::string& name, const T* data_ptr);
  
  /**
   * Get the data specified by \p name.  The return type is
   * specified by the template parameter.  Note that \p T
   * must be the same type that was used in the \p add_data()
   * member.
   */
  template <typename T>
  const T & get_data (const std::string& name) const;
  
  /**
   * Get the data specified by \p name.  The return type is
   * specified by the template parameter.  Note that \p T
   * must be the same type that was used in the \p add_data()
   * member.
   */
  template <class T>
  const T & get_data_object (const std::string& name) const;
  
//   /**
//    * Get the data specified by \p name.  The return type is
//    * specified by the template parameter.  Note that \p T
//    * may be any type that is related to the type used
//    * in the \p add_data() member.  
//    */
//   template <typename T>
//   const T * get_data_ptr (const std::string& name) const;

  
private:

  
  /**
   * Data container.
   */
  std::map<std::string, DataObject*> _data_map;
};




//----------------------------------------------------------------
// DataMap inline members
template <typename T>
inline
DataObject* DataObjectCpy<T>::build (const T& data)
{
  return new DataObjectCpy<T> (data);
}



// template <typename T>
// inline
// DataObject* DataObjectPtr<T>::build (const T* data_ptr)
// {
//   assert (data_ptr != NULL);
//   return new DataObjectPtr<T> (data_ptr);
// }



inline
void DataMap::add_data_object (const std::string& name, DataObject& data)
{
  here();
  
  // Check to see if we have somethint with this name already
  std::map<std::string, DataObject*>::iterator
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
      _data_map[name] = &data;
    }
}



template <typename T>
inline
void DataMap::add_data (const std::string& name, const T& data)
{
  // Check to see if we have somethint with this name already
  std::map<std::string, DataObject*>::iterator
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
      _data_map[name] = DataObjectCpy<T>::build (data);
    }
}



// template <typename T>
// inline
// void DataMap::add_data_ptr (const std::string& name, const T* data_ptr)
// {
//   // Check to see if we have somethint with this name already
//   std::map<std::string, DataObject*>::iterator
//     pos = _data_map.find (name);


//   // Make sure it wasn't already there
//   if (pos != _data_map.end())
//     {
//       std::cerr << "ERROR:  Data already exists for \""
// 		<< name
// 		<< "\", aborting!"
// 		<< std::endl;

//       error();
//     }
  
//   // Didn't find it
//   else
//     {
//       _data_map[name] = DataObjectPtr<T>::build (data_ptr);
//     }
// }



template <typename T>
inline
const T& DataMap::get_data (const std::string& name) const
{
  std::map<std::string, DataObject*>::const_iterator
    pos = _data_map.find(name);
  
  if (pos == _data_map.end())
    {
      std::cerr << "ERROR:  No data associated with the name \""
		<< name << "\"!"
		<< std::endl;

      error();
    }
  
  assert (pos->second != NULL);


  // Maybe it was a copy
  DataObjectCpy<T>* obj = dynamic_cast<DataObjectCpy<T>* >(pos->second);

  // Check for failed cast
  if (obj == NULL)
    {      
      std::cerr << "ERROR: Cast failed, invalid type requested!"
		<< "The requested type is: \"" << typeid(T).name() << "\""
		<< std::endl;
      error();
    }


  return obj->get_data ();
}



template <class T>
inline
const T& DataMap::get_data_object (const std::string& name) const
{
  std::map<std::string, DataObject*>::const_iterator
    pos = _data_map.find(name);
  
  if (pos == _data_map.end())
    {
      std::cerr << "ERROR:  No data associated with the name \""
		<< name << "\"!"
		<< std::endl;

      error();
    }
  
  assert (pos->second != NULL);


  // Try this first
  T* obj = dynamic_cast<T*>(pos->second);

  // Check for failed cast
  if (obj == NULL)
    {      
      std::cerr << "ERROR: Cast failed, invalid type requested!"
		<< "The requested type is: \"" << typeid(T).name() << "\""
		<< std::endl;
      error();
    }

  return *obj;
}



// template <typename T>
// inline
// const T* DataMap::get_data_ptr (const std::string& name) const
// {
//   std::map<std::string, DataObject*>::const_iterator
//     pos = _data_map.find(name);
  
//   if (pos == _data_map.end())
//     {
//       std::cerr << "ERROR:  No data associated with the name \""
// 		<< name << "\"!"
// 		<< std::endl;

//       error();
//     }


//   // Can't happen...
//   assert (pos->second != NULL);
  
//   // the object must be there.  Look above
//   T* obj = dynamic_cast<T*>(pos->second);

//   // Check for failed cast
//   if (obj == NULL)
//     {      
//       std::cerr << "ERROR: Cast failed, invalid type requested!"
// 		<< "The requested type is: \"" << typeid(T).name() << "\""
// 		<< std::endl;
//       error();
//     }

//   return obj;
// }


#endif // #define __data_map_h__
