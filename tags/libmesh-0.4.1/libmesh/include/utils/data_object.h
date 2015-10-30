// $Id: data_object.h,v 1.1 2003-11-05 22:26:44 benkirk Exp $

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



#ifndef __data_object_h__
#define __data_object_h__



//----------------------------------------
// DataObject
  
/**
 * The base class for a generic data object.
 */
class DataObject
{
public:

  /**
   * Destructor. This class is pure-virtual.
   */
  virtual ~DataObject () = 0;
};


inline
DataObject::~DataObject () {}


#endif // #define __data_map_h__
