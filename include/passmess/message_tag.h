// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef PASSMESS_MESSAGE_TAG_H
#define PASSMESS_MESSAGE_TAG_H

// libMesh Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_MPI
#  include "mpi.h"
#endif // LIBMESH_HAVE_MPI

// C++ includes
#include <climits> // INT_MIN

namespace PassMess
{
//-------------------------------------------------------------------
/**
 * Forward declarations of classes we will define later.
 */
class Communicator;

//-------------------------------------------------------------------
/**
 * Encapsulates the MPI tag integers.
 */
class MessageTag
{
public:

  /**
   * Invalid tag, to allow for default construction.
   */
  static const int invalid_tag = INT_MIN;

  /**
   * Explicit constructor, to discourage using "magic numbers"
   * as tags.  Communicator::get_unique_tag is recommended instead.
   */
  explicit MessageTag(int tagvalue = invalid_tag)
    : _tagvalue(tagvalue), _comm(nullptr) {}

  /**
   * Copy constructor.  Helps Communicator do reference counting on
   * unique tags
   */
  MessageTag(const MessageTag & other);

  /**
   * Move constructor.  Helps Communicator do reference counting on
   * unique tags
   */
  MessageTag(MessageTag && other);

  /**
   * Copy assignment operator.  Helps Communicator do reference
   * counting on unique tags
   */
  MessageTag & operator = (const MessageTag & other);

  /**
   * Move assignment operator.  Helps Communicator do reference
   * counting on unique tags
   */
  MessageTag & operator = (MessageTag && other);

  /**
   * Destructor.  Helps Communicator do reference counting on unique
   * tags
   */
  ~MessageTag();

  int value() const {
    return _tagvalue;
  }

private:
  int _tagvalue;
  const Communicator * _comm;

  // Constructor for reference-counted unique tags
  MessageTag(int tagvalue, const Communicator * comm)
    : _tagvalue(tagvalue), _comm(comm) {}

  // Let Communicator handle the reference counting
  friend class Communicator;
};


//-------------------------------------------------------------------
/**
 * Default message tag ids
 */
#ifdef LIBMESH_HAVE_MPI
const MessageTag any_tag = MessageTag(MPI_ANY_TAG);
#else
const MessageTag any_tag = MessageTag(-1);
#endif

const MessageTag no_tag = MessageTag(0);

} // namespace PassMess

#endif // PASSMESS_MESSAGE_TAG_H
