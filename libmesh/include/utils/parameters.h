// $Id: parameters.h,v 1.1 2004-11-19 13:57:52 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __parameters_h__
#define __parameters_h__

// C++ includes
#include <string>
#include <map>

// Local includes
#include "libmesh_common.h"



/**
 * This class provides the ability to map between
 * arbitrary, user-defined strings and several data
 * types.  This can be used to provide arbitrary
 * user-specified options.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision: 1.1 $
 */

// ------------------------------------------------------------
// Parameters class definition
class Parameters 
{
public:

  template <typename T>
  bool have_parameter (const std::string&) const;

  template <typename T>
  const T& get (const std::string&) const;
  
  template <typename T>
  T& set (const std::string&);

  void clear () { error(); }
  
private:

  mutable std::map<std::string, int>         _int;
  mutable std::map<std::string, Real>        _real;
  mutable std::map<std::string, Complex>     _complex;
  mutable std::map<std::string, std::string> _string;
};



// ------------------------------------------------------------
// Parameters class inline methods
template <>
inline
bool Parameters::have_parameter<int>(const std::string& name) const
{
  return _int.count(name);
}



template <>
inline
const int& Parameters::get (const std::string& name) const
{
  if (!this->have_parameter<int>(name))
    {
      std::cerr << "ERROR: no integer parameter named "
		<< name << std::endl;
      error();
    }
  
  return _int[name];
}



template <>
inline
int& Parameters::set (const std::string& name)
{
  return _int[name];
}



#endif // #define __parameters_h__
