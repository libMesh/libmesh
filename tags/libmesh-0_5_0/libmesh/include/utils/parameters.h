// $Id: parameters.h,v 1.5 2005-02-22 22:17:35 jwpeterson Exp $

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



#ifndef __parameters_h__
#define __parameters_h__

// C++ includes
#include <typeinfo>
#include <string>
#include <map>

// Local includes
#include "libmesh_common.h"
#include "reference_counted_object.h"



/**
 * This class provides the ability to map between
 * arbitrary, user-defined strings and several data
 * types.  This can be used to provide arbitrary
 * user-specified options.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \version $Revision: 1.5 $
 */

// ------------------------------------------------------------
// Parameters class definition
class Parameters 
{
public:

  /**
   * Default constructor.  Does nothing.
   */
  Parameters () {};

  /**
   * Copy constructor.
   */
  Parameters (const Parameters&);

  /**
   * Destructor.  Clears any allocated memory.
   */
  ~Parameters ();

  /**
   * Assignment operator.
   */
  Parameters& operator= (const Parameters&);
  
  /**
   * @returns \p true if a parameter of type \p T
   * with a specified name exists, \p false otherwise.
   */
  template <typename T>
  bool have_parameter (const std::string&) const;

  /**
   * @returns a constant reference to the specified parameter
   * value.  Requires, of course, that the parameter exists.
   */
  template <typename T>
  const T& get (const std::string&) const;

  /**
   * @returns a writeable reference to the specified parameter.
   * This method will create the parameter if it does not exist,
   * so it can be used to define parameters which will later be
   * accessed with the \p get() member.
   */
  template <typename T>
  T& set (const std::string&);

  /**
   * Removes the specified parameter from the list, if it exists.
   */
  void remove (const std::string&);

  /**
   * @returns the total number of parameters.
   */
  unsigned int n_parameters () const { return _values.size(); }

  /**
   * @returns the number of parameters of the requested type.
   */
  template <typename T>
  unsigned int n_parameters () const;
  
  /**
   * Clears internal data structures & frees any allocated memory.
   */
  void clear ();

  /**
   * Prints the contents to the specified stream.
   */
  void print (std::ostream& os=std::cout) const;

  
private:

  /**
   * Abstract definition of a parameter value.
   */
  class Value : public ReferenceCountedObject<Value>
  {
  public:

    /**
     * Destructor.
     */
    virtual ~Value() {};

    /**
     * String identifying the type of parameter stored.
     * Must be reimplemented in derived classes.
     */
    virtual std::string type () const = 0;

    /**
     * Prints the parameter value to the specified stream.
     * Must be reimplemented in derived classes.
     */
    virtual void print(std::ostream&) const = 0;

    /**
     * Clone this value.  Useful in copy-construction.
     * Must be reimplemented in derived classes.
     */
    virtual Value* clone () const = 0;
  };


  /**
   * Concrete definition of a parameter value
   * for a specified type.
   */
  template <typename T>
  class Parameter : public Value
  {
  public:

    /**
     * @returns a read-only reference to the parameter value.
     */
    const T& get () const { return _value; }
    
    /**
     * @returns a writeable reference to the parameter value.
     */
    T& set () { return _value; }

    /**
     * String identifying the type of parameter stored.
     */
    virtual std::string type () const;

    /**
     * Prints the parameter value to the specified stream.
     */
    virtual void print(std::ostream&) const;

    /**
     * Clone this value.  Useful in copy-construction.
     */
    virtual Value* clone () const;
    
  private:

    /**
     * Stored parameter value.
     */
    T _value;
  };

  /**
   * Data structure to map names with values.
   */
  std::map<std::string, Value*> _values;

  /**
   * Parameter map iterator.
   */
  typedef std::map<std::string, Value*>::iterator iterator;

  /**
   * Constant parameter map iterator.
   */
  typedef std::map<std::string, Value*>::const_iterator const_iterator;
};



// ------------------------------------------------------------
// Parameters::Parameter<> class inline methods
template <typename T>
inline
std::string Parameters::Parameter<T>::type () const
{
  return typeid(T).name();
}



template <typename T>
inline
void Parameters::Parameter<T>::print (std::ostream& os) const
{
  os << _value;
}



template <typename T>
inline
Parameters::Value* Parameters::Parameter<T>::clone () const
{
  Parameters::Parameter<T>
    *copy = new Parameters::Parameter<T>;

  assert (copy != NULL);
  
  copy->_value = _value;

  return copy;
}


// ------------------------------------------------------------
// Parameters class inline methods
inline
void Parameters::clear () // since this is inline we must define it
{                         // before its first use (for some compilers)
  while (!_values.empty())
    {
      Parameters::iterator it = _values.begin();
      
      delete it->second;      
      it->second = NULL;
      
      _values.erase(it);
    }
}



inline
Parameters& Parameters::operator= (const Parameters& rhs)
{
  this->clear();
  
  for (Parameters::const_iterator it = rhs._values.begin();
       it != rhs._values.end(); ++it)
    _values[it->first] = it->second->clone();
  
  return *this;
}



inline
Parameters::Parameters (const Parameters& p)
{
  *this = p;
}



inline
Parameters::~Parameters ()
{
  this->clear ();
}



inline
void Parameters::print (std::ostream& os) const
{
  Parameters::const_iterator it = _values.begin();

  os << "Name\t Type\t Value\n"
     << "---------------------\n";
  while (it != _values.end())
    {
      os << " "   << it->first
	 << "\t " << it->second->type()
	 << "\t ";   it->second->print(os);      
      os << '\n';

      ++it;
    }
}



template <typename T>
inline
bool Parameters::have_parameter (const std::string& name) const
{
  Parameters::const_iterator it = _values.find(name);

  if (it != _values.end())
    if (dynamic_cast<const Parameter<T>*>(it->second) != NULL)
      return true;

  return false;
}



template <typename T>
inline
const T& Parameters::get (const std::string& name) const
{
  if (!this->have_parameter<T>(name))
    {
      std::cerr << "ERROR: no "
		<< typeid(T).name()
		<< " parameter named "
		<< name << ":" << std::endl
		<< *this;
      
      error();
    }

  Parameters::const_iterator it = _values.find(name);

  assert (it != _values.end());
  assert (it->second != NULL);
  assert (dynamic_cast<const Parameter<T>*>(it->second) != NULL);
  
  return dynamic_cast<Parameter<T>*>(it->second)->get();
}



template <typename T>
inline
T& Parameters::set (const std::string& name)
{
  Parameter<T>* param = NULL;
  
  if (!this->have_parameter<T>(name))
    _values[name] = new Parameter<T>;

  param = dynamic_cast<Parameter<T>*>(_values[name]);

  assert (param != NULL);
  
  return param->set();
}



inline
void Parameters::remove (const std::string& name)
{
  Parameters::iterator it = _values.find(name);

  if (it != _values.end())
    {
      delete it->second;
      it->second = NULL;

      _values.erase(it);
    }
}



template <typename T>
inline
unsigned int Parameters::n_parameters () const
{
  unsigned int cnt = 0;  
  
  Parameters::const_iterator       it  = _values.begin();
  const Parameters::const_iterator end = _values.end();
  
  for (; it != end; ++it)
    if (dynamic_cast<Parameter<T>*>(it->second) != NULL)
      cnt++;

  return cnt;	 
}


inline
std::ostream& operator << (std::ostream& os, const Parameters& p)
{
  p.print(os);
  return os;
}

#endif // #define __parameters_h__
