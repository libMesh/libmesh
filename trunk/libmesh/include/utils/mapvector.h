// $Id: mapvector.h,v 1.1 2007-10-22 19:57:57 roystgnr Exp $

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



#ifndef __mapvector_h__
#define __mapvector_h__

// C++ Includes   -----------------------------------
#include <map>

template <typename Val>
class mapvector : public std::map<unsigned int, Val>
{
public:
  typedef std::map<unsigned int, Val> maptype;

  Val& operator[] (const unsigned int &k)
  {
    return maptype::operator[](k);
  }
  Val operator[] (const unsigned int &k) const 
  {
    typename maptype::const_iterator it = this->find(k);
      return it == this->end().it? Val() : it->second;
  }

  class veclike_iterator
  {
  public:
    veclike_iterator(const typename maptype::iterator &i)
      : it(i) {}

    veclike_iterator(const veclike_iterator &i)
      : it(i.it) {}

    Val& operator*() const { return it->second; }

    veclike_iterator& operator++() { ++it; return *this; }

    veclike_iterator operator++(int) {
      veclike_iterator i = *this;
      ++(*this);
      return i;
    }

    bool operator==(const veclike_iterator &other) const {
      return it == other.it;
    }

    bool operator!=(const veclike_iterator &other) const {
      return it != other.it;
    }

    typename maptype::iterator it;
  };

  class const_veclike_iterator
  {
  public:
    const_veclike_iterator(const typename maptype::const_iterator &i)
      : it(i) {}

    const_veclike_iterator(const const_veclike_iterator &i)
      : it(i.it) {}

    const_veclike_iterator(const veclike_iterator &i)
      : it(i.it) {}

    const Val& operator*() const { return it->second; }

    const_veclike_iterator& operator++() { ++it; return *this; }

    const_veclike_iterator operator++(int) {
      veclike_iterator i = *this;
      ++(*this);
      return i;
    }

    bool operator==(const const_veclike_iterator &other) const {
      return it == other.it;
    }

    bool operator!=(const const_veclike_iterator &other) const {
      return it != other.it;
    }

    typename maptype::const_iterator it;
  };

  void erase(unsigned int i) {
      maptype::erase(i);
  }

  void erase(const veclike_iterator &pos) {
      maptype::erase(pos.it);
  }

  veclike_iterator begin() {
    return veclike_iterator(maptype::begin());
  }

  const_veclike_iterator begin() const {
    return const_veclike_iterator(maptype::begin());
  }

  veclike_iterator end() {
    return veclike_iterator(maptype::end());
  }

  const_veclike_iterator end() const {
    return const_veclike_iterator(maptype::end());
  }
};

#endif // __mapvector_h__
