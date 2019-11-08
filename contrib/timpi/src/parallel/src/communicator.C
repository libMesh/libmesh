// The TIMPI Message-Passing Parallelism Library.
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


// Local includes
#include "timpi/communicator.h"

// TIMPI includes
#include "timpi/parallel_implementation.h" // for inline max(int)
#include "timpi/timpi_assert.h"

// Disable libMesh logging until we decide how to port it best
// #include "libmesh/libmesh_logging.h"
#define TIMPI_LOG_SCOPE(f,c)

namespace TIMPI
{


// ------------------------------------------------------------
// Simple Communicator member functions

void Communicator::reference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  timpi_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]++;
}


void Communicator::dereference_unique_tag(int tagvalue) const
{
  // This had better be an already-acquired tag.
  timpi_assert(used_tag_values.count(tagvalue));

  used_tag_values[tagvalue]--;
  // If we don't have any more outstanding references, we
  // don't even need to keep this tag in our "used" set.
  if (!used_tag_values[tagvalue])
    used_tag_values.erase(tagvalue);
}


Communicator::Communicator () :
#ifdef TIMPI_HAVE_MPI
  _communicator(MPI_COMM_SELF),
#endif
  _rank(0),
  _size(1),
  _send_mode(DEFAULT),
  used_tag_values(),
  _next_tag(0),
  _max_tag(std::numeric_limits<int>::max()),
  _I_duped_it(false) {}


Communicator::Communicator (const communicator & comm) :
#ifdef TIMPI_HAVE_MPI
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


#ifdef TIMPI_HAVE_MPI
void Communicator::split(int color, int key, Communicator & target) const
{
  target.clear();
  MPI_Comm newcomm;
  timpi_call_mpi
    (MPI_Comm_split(this->get(), color, key, &newcomm));

  target.assign(newcomm);
  target._I_duped_it = (color != MPI_UNDEFINED);
  target.send_mode(this->send_mode());
}


void Communicator::split_by_type(int split_type, int key, info i, Communicator & target) const
{
  target.clear();
  MPI_Comm newcomm;
  timpi_call_mpi
    (MPI_Comm_split_type(this->get(), split_type, key, i, &newcomm));

  target.assign(newcomm);
  target._I_duped_it = (split_type != MPI_UNDEFINED);
  target.send_mode(this->send_mode());
}

#else
void Communicator::split(int, int, Communicator & target) const
{
  target.assign(this->get());
}

void Communicator::split_by_type(int, int, info, Communicator & target) const
{
  target.assign(this->get());
}
#endif


void Communicator::duplicate(const Communicator & comm)
{
  this->duplicate(comm._communicator);
  this->send_mode(comm.send_mode());
}


#ifdef TIMPI_HAVE_MPI
void Communicator::duplicate(const communicator & comm)
{
  if (_communicator != MPI_COMM_NULL)
    {
      timpi_call_mpi
        (MPI_Comm_dup(comm, &_communicator));

      _I_duped_it = true;
    }
  this->assign(_communicator);
}
#else
void Communicator::duplicate(const communicator &) { }
#endif


void Communicator::clear() {
#ifdef TIMPI_HAVE_MPI
  if (_I_duped_it)
    {
      timpi_assert (_communicator != MPI_COMM_NULL);
      timpi_call_mpi
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
#ifdef TIMPI_HAVE_MPI
  if (_communicator != MPI_COMM_NULL)
    {
      int i;
      timpi_call_mpi
        (MPI_Comm_size(_communicator, &i));

      timpi_assert_greater_equal (i, 0);
      _size = cast_int<processor_id_type>(i);

      timpi_call_mpi
        (MPI_Comm_rank(_communicator, &i));

      timpi_assert_greater_equal (i, 0);
      _rank = cast_int<processor_id_type>(i);

      int * maxTag;
      int flag = false;
      timpi_call_mpi(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &maxTag, &flag));
      timpi_assert(flag);
      _max_tag = *maxTag;
    }
  else
    {
      _rank = 0;
      _size = 1;
      _max_tag = std::numeric_limits<int>::max();
    }
  _next_tag = _max_tag / 2;
#endif
  _send_mode = DEFAULT;
}


/**
 * Pause execution until all processors reach a certain point.
 */
#ifdef TIMPI_HAVE_MPI
void Communicator::barrier () const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("barrier()", "Communicator");
      timpi_call_mpi(MPI_Barrier (this->get()));
    }
}
#else
void Communicator::barrier () const {}
#endif

#ifdef TIMPI_HAVE_MPI
void Communicator::nonblocking_barrier (Request & req) const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("nonblocking_barrier()", "Communicator");
      timpi_call_mpi(MPI_Ibarrier (this->get(), req.get()));
    }
}
#else
void Communicator::nonblocking_barrier (Request & /*req*/) const {}
#endif


MessageTag Communicator::get_unique_tag(int tagvalue) const
{
  if (tagvalue == MessageTag::invalid_tag)
    {
#ifndef NDEBUG
      // Automatic tag values have to be requested in sync
      int maxval = _next_tag;
      this->max(maxval);
      timpi_assert_equal_to(_next_tag, maxval);
#endif
      tagvalue = _next_tag++;
    }

  if (used_tag_values.count(tagvalue))
    {
      // Get the largest value in the used values, and pick one
      // larger
      tagvalue = used_tag_values.rbegin()->first+1;
      timpi_assert(!used_tag_values.count(tagvalue));
    }
  if (tagvalue >= _next_tag)
    _next_tag = tagvalue+1;

  if (_next_tag >= _max_tag)
    _next_tag = _max_tag/2;

  used_tag_values[tagvalue] = 1;

  return MessageTag(tagvalue, this);
}


status Communicator::probe (const unsigned int src_processor_id,
                            const MessageTag & tag) const
{
  TIMPI_LOG_SCOPE("probe()", "Communicator");

#ifndef TIMPI_HAVE_MPI
  timpi_not_implemented();
  ignore(src_processor_id, tag);
#endif

  status stat;

  timpi_assert(src_processor_id < this->size() ||
                  src_processor_id == any_source);

  timpi_call_mpi
    (MPI_Probe (src_processor_id, tag.value(), this->get(), &stat));

  return stat;
}


bool Communicator::verify(const bool & r) const
{
  const unsigned char rnew = r;
  return this->verify(rnew);
}


bool Communicator::semiverify(const bool * r) const
{
  if (r)
    {
      const unsigned char rnew = *r;
      return this->semiverify(&rnew);
    }

  const unsigned char * rptr = nullptr;
  return this->semiverify(rptr);
}


bool Communicator::verify(const std::string & r) const
{
  if (this->size() > 1)
    {
      // Cannot use <char> since MPI_MIN is not
      // strictly defined for chars!
      std::vector<short int> temp; temp.reserve(r.size());
      for (std::size_t i=0; i != r.size(); ++i)
        temp.push_back(r[i]);
      return this->verify(temp);
    }
  return true;
}


bool Communicator::semiverify(const std::string * r) const
{
  if (this->size() > 1)
    {
      std::size_t rsize = r ? r->size() : 0;
      std::size_t * psize = r ? &rsize : nullptr;

      if (!this->semiverify(psize))
        return false;

      this->max(rsize);

      // Cannot use <char> since MPI_MIN is not
      // strictly defined for chars!
      std::vector<short int> temp (rsize);
      if (r)
        {
          temp.reserve(rsize);
          for (std::size_t i=0; i != rsize; ++i)
            temp.push_back((*r)[i]);
        }

      std::vector<short int> * ptemp = r ? &temp: nullptr;

      return this->semiverify(ptemp);
    }
  return true;
}


void Communicator::min(bool & r) const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("min(bool)", "Communicator");

      unsigned int temp = r;
      timpi_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &temp, 1,
                        StandardType<unsigned int>(),
                        OpFunction<unsigned int>::min(),
                        this->get()));
      r = temp;
    }
}


void Communicator::minloc(bool & r,
                          unsigned int & min_id) const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("minloc(bool)", "Communicator");

      DataPlusInt<int> data_in;
      ignore(data_in); // unused ifndef TIMPI_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out = data_in;

      timpi_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type_acquire<int>().first,
                        OpFunction<int>::min_location(), this->get()));
      r = data_out.val;
      min_id = data_out.rank;
    }
  else
    min_id = this->rank();
}


void Communicator::max(bool & r) const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("max(bool)", "Communicator");

      unsigned int temp = r;
      timpi_call_mpi
        (MPI_Allreduce (MPI_IN_PLACE, &temp, 1,
                        StandardType<unsigned int>(),
                        OpFunction<unsigned int>::max(),
                        this->get()));
      r = temp;
    }
}


void Communicator::maxloc(bool & r,
                          unsigned int & max_id) const
{
  if (this->size() > 1)
    {
      TIMPI_LOG_SCOPE("maxloc(bool)", "Communicator");

      DataPlusInt<int> data_in;
      ignore(data_in); // unused ifndef TIMPI_HAVE_MPI
      data_in.val = r;
      data_in.rank = this->rank();
      DataPlusInt<int> data_out = data_in;

      timpi_call_mpi
        (MPI_Allreduce (&data_in, &data_out, 1,
                        dataplusint_type_acquire<int>().first,
                        OpFunction<int>::max_location(),
                        this->get()));
      r = data_out.val;
      max_id = data_out.rank;
    }
  else
    max_id = this->rank();
}


} // namespace TIMPI
