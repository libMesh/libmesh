// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_HASHWORD_H
#define LIBMESH_HASHWORD_H

// All functions in this file by Bob Jenkins, May 2006, Public Domain.

#include <stdint.h> // uint32_t

// ::size_t is defined in the backward compatibility header stddef.h.
// It's been part of ANSI/ISO C and ISO C++ since their very
// beginning. Every C++ implementation has to ship with stddef.h
// (compatibility) and cstddef where only the latter defines
// std::size_t and not necessarily ::size_t. See Annex D of the C++
// Standard.
//
// http://stackoverflow.com/questions/237370/does-stdsize-t-make-sense-in-c
#include <stddef.h>


// Anonymous namespace for utility functions used locally
namespace
{
// Rotate x by k bits
// @author Bob Jenkins, May 2006, Public Domain.
// http://burtleburtle.net/bob/hash/index.html
inline
uint32_t rot(uint32_t x, uint32_t k)
{
  return (x<<k) | (x>>(32-k));
}



// mix 3 32-bit values reversibly
// @author Bob Jenkins, May 2006, Public Domain.
// http://burtleburtle.net/bob/hash/index.html
inline
void mix(uint32_t& a, uint32_t& b, uint32_t& c)
{
  a -= c;  a ^= rot(c, 4);  c += b;
  b -= a;  b ^= rot(a, 6);  a += c;
  c -= b;  c ^= rot(b, 8);  b += a;
  a -= c;  a ^= rot(c,16);  c += b;
  b -= a;  b ^= rot(a,19);  a += c;
  c -= b;  c ^= rot(b, 4);  b += a;
}


// 'final' mixing of 3 32-bit numbers, result is stored in c.
// @author Bob Jenkins, May 2006, Public Domain.
// http://burtleburtle.net/bob/hash/index.html
inline
void final(uint32_t& a, uint32_t& b, uint32_t& c)
{
  c ^= b; c -= rot(b,14);
  a ^= c; a -= rot(c,11);
  b ^= a; b -= rot(a,25);
  c ^= b; c -= rot(b,16);
  a ^= c; a -= rot(c,4);
  b ^= a; b -= rot(a,14);
  c ^= b; c -= rot(b,24);
}
} // end anonymous namespace



namespace libMesh
{
namespace Utility
{
// The hashword function takes an array of uint32_t's of length 'length'
// and computes a single key from it.
// @author Bob Jenkins, May 2006, Public Domain.
// http://burtleburtle.net/bob/hash/index.html
inline
uint32_t hashword(const uint32_t *k, size_t length, uint32_t initval=0)
{
  uint32_t a,b,c;

  // Set up the internal state
  a = b = c = 0xdeadbeef + ((static_cast<uint32_t>(length))<<2) + initval;

  //------------------------------------------------- handle most of the key
  while (length > 3)
    {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      length -= 3;
      k += 3;
    }

  //------------------------------------------- handle the last 3 uint32_t's
  switch(length)                     // all the case statements fall through
    {
    case 3 : c+=k[2];
    case 2 : b+=k[1];
    case 1 : a+=k[0];
      final(a,b,c);
    default:     // case 0: nothing left to add
      break;
    }

  //------------------------------------------------------ report the result
  return c;
}


// This is a hard-coded version of hashword for hashing exactly 2 numbers
// @author Bob Jenkins, May 2006, Public Domain.
// http://burtleburtle.net/bob/hash/index.html
inline
uint32_t hashword2(const uint32_t& first, const uint32_t& second, uint32_t initval=0)
{
  uint32_t a,b,c;

  // Set up the internal state
  a = b = c = 0xdeadbeef + 8 + initval;

  b+=second;
  a+=first;
  final(a,b,c);

  return c;
}

// Homegrown implementations if we don't have 32-bit unsigned ints
inline
uint64_t hashword2(const uint64_t first, const uint64_t second)
{
  // big prime number
  const unsigned int bp = 65449;

  return (first%bp + (second<<5)%bp);
}

inline
uint16_t hashword2(const uint16_t first, const uint16_t second)
{
  // "big" prime number
  const uint16_t bp = 257;

  return static_cast<uint16_t>(first%bp + (second<<3)%bp);
}

// Another homegrown implementation for 64-bit unsigned ints
inline
uint64_t hashword(const uint64_t *k, size_t length)
{
  // big prime number
  const unsigned int bp = 65449;

  uint64_t c = 0;
  unsigned int shift=0;
  for (size_t i=0; i != length; ++i)
    {
      c += (k[i] << shift) % bp;
      shift += 5;
    }

  return c;
}

// Another homegrown implementation for 16-bit unsigned ints
inline
uint16_t hashword(const uint16_t *k, size_t length)
{
  // "big" prime number
  const uint16_t bp = 257;

  uint16_t c = 0;
  uint16_t shift=0;
  for (size_t i=0; i != length; ++i)
    {
      c = static_cast<uint16_t>
        (c + static_cast<uint16_t>(k[i] << shift) % bp);
      shift = static_cast<uint16_t>(shift+3);
    }

  return c;
}

} // end Utility namespace
} // end libMesh namespace

#endif // LIBMESH_HASHWORD_H
