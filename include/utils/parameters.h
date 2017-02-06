// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARAMETERS_H
#define LIBMESH_PARAMETERS_H

// C++ includes
#include <typeinfo>
#include <string>
#include <map>

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/print_trace.h"

// C++ includes
#include <cstddef>
#include <map>
#include <string>
#include <typeinfo>
#include <sstream>

namespace libMesh
{
/**
 * Helper functions for printing scalar, vector and vector<vector> types.  Called from Parameters::Parameter<T>::print(...).
 */
template<typename P>
void print_helper(std::ostream & os, const P * param);

template<typename P>
void print_helper(std::ostream & os, const std::vector<P> * param);

template<typename P>
void print_helper(std::ostream & os, const std::vector<std::vector<P> > * param);

/**
 * This class provides the ability to map between
 * arbitrary, user-defined strings and several data
 * types.  This can be used to provide arbitrary
 * user-specified options.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 */
class Parameters
{
public:

  /**
   * Default constructor.  Does nothing.
   */
  Parameters () {}

  /**
   * Copy constructor.
   */
  Parameters (const Parameters &);

  /**
   * Destructor.  Clears any allocated memory.
   */
  virtual ~Parameters ();

  /**
   * Assignment operator.  Removes all parameters in \p this
   * and inserts copies of all parameters from \p source
   */
  virtual Parameters & operator= (const Parameters & source);

  /**
   * Addition/Assignment operator.  Inserts copies of all parameters
   * from \p source.  Any parameters of the same name already in \p
   * this are replaced.
   */
  virtual Parameters & operator+= (const Parameters & source);

  /**
   * @returns \p true if a parameter of type \p T
   * with a specified name exists, \p false otherwise.
   *
   * If RTTI has been disabled then we return \p true
   * if a parameter of specified name exists regardless of its type.
   */
  template <typename T>
  bool have_parameter (const std::string &) const;

  /**
   * @returns a constant reference to the specified parameter
   * value.  Requires, of course, that the parameter exists.
   */
  template <typename T>
  const T & get (const std::string &) const;

  /**
   * Inserts a new Parameter into the object but does not return
   * a writable reference.  The value of the newly inserted
   * parameter may not be valid.
   */
  template <typename T>
  void insert (const std::string &);

  /**
   * @returns a writeable reference to the specified parameter.
   * This method will create the parameter if it does not exist,
   * so it can be used to define parameters which will later be
   * accessed with the \p get() member.
   */
  template <typename T>
  T & set (const std::string &);

  /**
   * Overridable function to set any extended attributes for
   * classes inheriting from this class.
   */
  virtual void set_attributes(const std::string &, bool /*inserted_only*/) {}

  /**
   * Removes the specified parameter from the list, if it exists.
   */
  void remove (const std::string &);

  /**
   * @returns the total number of parameters.
   */
  std::size_t n_parameters () const { return _values.size(); }

#ifdef LIBMESH_HAVE_RTTI
  /**
   * @returns the number of parameters of the requested type.
   */
  template <typename T>
  unsigned int n_parameters () const;
#endif // LIBMESH_HAVE_RTTI

  /**
   * Clears internal data structures & frees any allocated memory.
   */
  virtual void clear ();

  /**
   * Prints the contents, by default to libMesh::out.
   */
  void print (std::ostream & os=libMesh::out) const;

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
    virtual ~Value() {}

#ifdef LIBMESH_HAVE_RTTI
    /**
     * String identifying the type of parameter stored.
     * Must be reimplemented in derived classes.
     */
    virtual std::string type () const = 0;
#endif // LIBMESH_HAVE_RTTI

    /**
     * Prints the parameter value to the specified stream.
     * Must be reimplemented in derived classes.
     */
    virtual void print(std::ostream &) const = 0;

    /**
     * Clone this value.  Useful in copy-construction.
     * Must be reimplemented in derived classes.
     */
    virtual Value * clone () const = 0;
  };

public:

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
    const T & get () const { return _value; }

    /**
     * @returns a writeable reference to the parameter value.
     */
    T & set () { return _value; }

#ifdef LIBMESH_HAVE_RTTI
    /**
     * String identifying the type of parameter stored.
     */
    virtual std::string type () const;
#endif // LIBMESH_HAVE_RTTI

    /**
     * Prints the parameter value to the specified stream.
     */
    virtual void print(std::ostream &) const;

    /**
     * Clone this value.  Useful in copy-construction.
     */
    virtual Value * clone () const;

  private:
    /**
     * Stored parameter value.
     */
    T _value;
  };

  /**
   * Parameter map iterator.
   */
  typedef std::map<std::string, Value *>::iterator iterator;

  /**
   * Constant parameter map iterator.
   */
  typedef std::map<std::string, Value *>::const_iterator const_iterator;

  /**
   * Iterator pointing to the beginning of the set of parameters.
   */
  iterator begin();

  /**
   * Iterator pointing to the beginning of the set of parameters.
   */
  const_iterator begin() const;

  /**
   * Iterator pointing to the end of the set of parameters
   */
  iterator end();

  /**
   * Iterator pointing to the end of the set of parameters
   */
  const_iterator end() const;

protected:

  /**
   * Data structure to map names with values.
   */
  std::map<std::string, Value *> _values;

};

// ------------------------------------------------------------
// Parameters::Parameter<> class inline methods

// This only works with Run-Time Type Information, even though
// typeid(T) *should* be determinable at compile time regardless...
#ifdef LIBMESH_HAVE_RTTI
template <typename T>
inline
std::string Parameters::Parameter<T>::type () const
{
  return demangle(typeid(T).name());
}
#endif

template <typename T>
inline
void Parameters::Parameter<T>::print (std::ostream & os) const
{
  // Call helper function overloaded for basic scalar and vector types
  print_helper(os, static_cast<const T *>(&_value));
}

template <typename T>
inline
Parameters::Value * Parameters::Parameter<T>::clone () const
{
  Parameter<T> * copy = new Parameter<T>;

  libmesh_assert(copy);

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
      it->second = libmesh_nullptr;

      _values.erase(it);
    }
}



inline
Parameters & Parameters::operator= (const Parameters & source)
{
  this->clear();
  *this += source;

  return *this;
}

inline
Parameters & Parameters::operator+= (const Parameters & source)
{
  for (Parameters::const_iterator it = source._values.begin();
       it != source._values.end(); ++it)
    {
      if (_values.find(it->first) != _values.end())
        delete _values[it->first];
      _values[it->first] = it->second->clone();
    }

  return *this;
}

inline
Parameters::Parameters (const Parameters & p)
{
  *this = p;
}



inline
Parameters::~Parameters ()
{
  this->clear ();
}



inline
void Parameters::print (std::ostream & os) const
{
  Parameters::const_iterator it = _values.begin();

  os << "Name\t Type\t Value\n"
     << "---------------------\n";
  while (it != _values.end())
    {
      os << " "   << it->first
#ifdef LIBMESH_HAVE_RTTI
         << "\t " << it->second->type()
#endif // LIBMESH_HAVE_RTTI
         << "\t ";   it->second->print(os);
      os << '\n';

      ++it;
    }
}



// Declare this now that Paramers::print() is defined.
// By declaring this early we can use it in subsequent
// methods.  Required for gcc-4.0.2 -- 11/30/2005, BSK
inline
std::ostream & operator << (std::ostream & os, const Parameters & p)
{
  p.print(os);
  return os;
}



template <typename T>
inline
bool Parameters::have_parameter (const std::string & name) const
{
  Parameters::const_iterator it = _values.find(name);

  if (it != _values.end())
#ifdef LIBMESH_HAVE_RTTI
    if (dynamic_cast<const Parameter<T> *>(it->second) != libmesh_nullptr)
#else // LIBMESH_HAVE_RTTI
      if (cast_ptr<const Parameter<T> *>(it->second) != libmesh_nullptr)
#endif // LIBMESH_HAVE_RTTI
        return true;

  return false;
}



template <typename T>
inline
const T & Parameters::get (const std::string & name) const
{
  if (!this->have_parameter<T>(name))
    {
      std::ostringstream oss;

      oss << "ERROR: no";
#ifdef LIBMESH_HAVE_RTTI
      oss << ' ' << demangle(typeid(T).name());
#endif
      oss << " parameter named \""
          << name << "\" found.\n\n"
          << "Known parameters:\n"
          << *this;

      libmesh_error_msg(oss.str());
    }

  Parameters::const_iterator it = _values.find(name);

  libmesh_assert(it != _values.end());
  libmesh_assert(it->second);

  return cast_ptr<Parameter<T> *>(it->second)->get();
}

template <typename T>
inline
void Parameters::insert (const std::string & name)
{
  if (!this->have_parameter<T>(name))
    _values[name] = new Parameter<T>;

  set_attributes(name, true);
}


template <typename T>
inline
T & Parameters::set (const std::string & name)
{
  if (!this->have_parameter<T>(name))
    _values[name] = new Parameter<T>;

  set_attributes(name, false);

  return cast_ptr<Parameter<T> *>(_values[name])->set();
}

inline
void Parameters::remove (const std::string & name)
{
  Parameters::iterator it = _values.find(name);

  if (it != _values.end())
    {
      delete it->second;
      it->second = libmesh_nullptr;

      _values.erase(it);
    }
}



#ifdef LIBMESH_HAVE_RTTI
template <typename T>
inline
unsigned int Parameters::n_parameters () const
{
  unsigned int cnt = 0;

  Parameters::const_iterator       it  = _values.begin();
  const Parameters::const_iterator vals_end = _values.end();

  for (; it != vals_end; ++it)
    if (dynamic_cast<Parameter<T> *>(it->second) != libmesh_nullptr)
      cnt++;

  return cnt;
}
#endif

inline
Parameters::iterator Parameters::begin()
{
  return _values.begin();
}

inline
Parameters::const_iterator Parameters::begin() const
{
  return _values.begin();
}

inline
Parameters::iterator Parameters::end()
{
  return _values.end();
}

inline
Parameters::const_iterator Parameters::end() const
{
  return _values.end();
}

//non-member scalar print function
template<typename P>
void print_helper(std::ostream & os, const P * param)
{
  os << *param;
}

template<>
inline
void print_helper(std::ostream & os, const char * param)
{
  // Specialization so that we don't print out unprintable characters
  os << static_cast<int>(*param);
}

template<>
inline
void print_helper(std::ostream & os, const unsigned char * param)
{
  // Specialization so that we don't print out unprintable characters
  os << static_cast<int>(*param);
}

//non-member vector print function
template<typename P>
void print_helper(std::ostream & os, const std::vector<P> * param)
{
  for (std::size_t i=0; i<param->size(); ++i)
    os << (*param)[i] << " ";
}

//non-member vector<vector> print function
template<typename P>
void print_helper(std::ostream & os, const std::vector<std::vector<P> > * param)
{
  for (std::size_t i=0; i<param->size(); ++i)
    for (std::size_t j=0; j<(*param)[i].size(); ++j)
      os << (*param)[i][j] << " ";
}

} // namespace libMesh

#endif // LIBMESH_PARAMETERS_H
