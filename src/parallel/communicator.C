// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Parallel includes
#include "libmesh/communicator.h"

// libMesh includes
#include "libmesh/libmesh_logging.h"
#include "libmesh/print_trace.h"

namespace libMesh
{

namespace Parallel
{


// ------------------------------------------------------------
// Simple Communicator member functions

void Communicator::reference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  libmesh_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]++;
}


void Communicator::dereference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  libmesh_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]--;
  // If we don't have any more outstanding references, we
  // don't even need to keep this tag in our "used" set.
  if (!used_tag_values[tagvalue])
    used_tag_values.erase(tagvalue);
}


Communicator::Communicator () :
#ifdef LIBMESH_HAVE_MPI
  _communicator(MPI_COMM_SELF),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _next_tag(0),
  _max_tag(std::numeric_limits<int>::max()),
  _I_duped_it(false)
{
}


Communicator::Communicator (const communicator & comm) :
#ifdef LIBMESH_HAVE_MPI
  _communicator(MPI_COMM_SELF),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _next_tag(0),
  _max_tag(std::numeric_limits<int>::max()),
  _I_duped_it(false)
{
  this->assign(comm);
}


Communicator::~Communicator ()
{
  this->clear();
}


#ifdef LIBMESH_HAVE_MPI
void Communicator::split(int color, int key, Communicator & target) const
{
  target.clear();
  MPI_Comm newcomm;
  libmesh_call_mpi
    (MPI_Comm_split(this->get(), color, key, &newcomm));

  target.assign(newcomm);
  target._I_duped_it = (color != MPI_UNDEFINED);
  target.send_mode(this->send_mode());
}
#else
void Communicator::split(int, int, Communicator & target) const
{
  target.assign(this->get());
}
#endif


void Communicator::duplicate(const Communicator & comm)
{
  this->duplicate(comm._communicator);
  this->send_mode(comm.send_mode());
}


#ifdef LIBMESH_HAVE_MPI
void Communicator::duplicate(const communicator & comm)
{
  if (_communicator != MPI_COMM_NULL)
    {
      libmesh_call_mpi
        (MPI_Comm_dup(comm, &_communicator));

      _I_duped_it = true;
    }
  this->assign(_communicator);
}
#else
void Communicator::duplicate(const communicator &) { }
#endif


void Communicator::clear() {
#ifdef LIBMESH_HAVE_MPI
  if (_I_duped_it)
    {
      libmesh_assert (_communicator != MPI_COMM_NULL);
      libmesh_call_mpi
        (MPI_Comm_free(&_communicator));

      _communicator = MPI_COMM_NULL;
    }
  _I_duped_it = false;
#endif
}


Communicator & Communicator::operator= (const communicator & comm)
{
  this->clear();
  this->assign(comm);
  return *this;
}


void Communicator::assign(const communicator & comm)
{
  _communicator = comm;
#ifdef LIBMESH_HAVE_MPI
  if (_communicator != MPI_COMM_NULL)
    {
      int i;
      libmesh_call_mpi
        (MPI_Comm_size(_communicator, &i));

      libmesh_assert_greater_equal (i, 0);
      _size = cast_int<processor_id_type>(i);

      libmesh_call_mpi
        (MPI_Comm_rank(_communicator, &i));

      libmesh_assert_greater_equal (i, 0);
      _rank = cast_int<processor_id_type>(i);

      // Get the maximum tag value
      {
        MPI_Aint * maxTag;
        int flag;
        libmesh_call_mpi(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &maxTag, &flag));

        _max_tag = *maxTag;
      }
    }
  else
    {
      _rank = 0;
      _size = 1;
    }
#endif
  _send_mode = DEFAULT;
}


/**
 * Pause execution until all processors reach a certain point.
 */
#ifdef LIBMESH_HAVE_MPI
void Communicator::barrier () const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("barrier()", "Parallel");
      libmesh_call_mpi(MPI_Barrier (this->get()));
    }
}
#else
void Communicator::barrier () const {}
#endif

#ifdef LIBMESH_HAVE_MPI
void Communicator::nonblocking_barrier (Request & req) const
{
  if (this->size() > 1)
    {
      LOG_SCOPE("nonblocking_barrier()", "Parallel");
      libmesh_call_mpi(MPI_Ibarrier (this->get(), req.get()));
    }
}
#else
void Communicator::nonblocking_barrier (Request & /*req*/) const {}
#endif


MessageTag Communicator::get_unique_tag(int /*tagvalue*/) const
{
  // Don't give out tag values greater than what MPI can handle!
  if (_next_tag == _max_tag)
    _next_tag = 0;

  auto new_tag = _next_tag++;

  // We don't want to simply get the largest one
  // because it might be near MPI_TAG_UB and then
  // we can end up in a sticky situation
  // Most-likely we will amost never hit a
  // tag that is in use because we have
  // ~2 billion of them
  while(used_tag_values.count(new_tag))
    new_tag++;

  used_tag_values[new_tag] = 1;

  return MessageTag(new_tag, this);
}


} // namespace Parallel

} // namespace libMesh
