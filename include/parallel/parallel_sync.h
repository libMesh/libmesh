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



#ifndef LIBMESH_PARALLEL_SYNC_H
#define LIBMESH_PARALLEL_SYNC_H

// Local Includes
#include "libmesh/parallel.h"

// C++ includes
#include <map>
#include <type_traits>
#include <vector>


namespace libMesh
{



//--------------------------------------------------------------------------
namespace Parallel {

//------------------------------------------------------------------------
/**
 * Send and receive and act on vectors of data.
 *
 * The \p data map is indexed by processor ids as keys, and for each
 * processor id in the map there should be a vector of data to send.
 *
 * Data which is received from other processors will be operated on by
 * act_on_data(processor_id_type pid, const std::vector<datum> & data)
 *
 * No guarantee about operation ordering is made - this function will
 * attempt to act on data in the order in which it is received.
 *
 * All receives and actions are completed before this function
 * returns.
 *
 * Not all sends may have yet completed.  The supplied container of
 * Request objects, \p req, has more requests inserted, one for each
 * of the data sends.  These requests must be waited on before the \p
 * data map is deleted.
 */
template <typename MapToVectors,
          typename RequestContainer,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & data,
                               RequestContainer & reqs,
                               ActionFunctor & act_on_data);



/**
 * Send and receive and act on vectors of data.
 *
 * The \p data map is indexed by processor ids as keys, and for each
 * processor id in the map there should be a vector of data to send.
 *
 * Data which is received from other processors will be operated on by
 * act_on_data(processor_id_type pid, const std::vector<datum> & data);
 *
 * No guarantee about operation ordering is made - this function will
 * attempt to act on data in the order in which it is received.
 *
 * All communication and actions are complete when this function
 * returns.
 */
template <typename MapToVectors,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & data,
                               ActionFunctor & act_on_data);


/**
 * Send query vectors, receive and answer them with vectors of data,
 * then act on those answers.
 *
 * The \p data map is indexed by processor ids as keys, and for each
 * processor id in the map there should be a vector of query ids to send.
 *
 * Query data which is received from other processors will be operated
 * on by
 * gather_data(processor_id_type pid, const std::vector<id> & ids,
 *             std::vector<datum> & data)
 *
 * Answer data which is received from other processors will be operated on by
 * act_on_data(processor_id_type pid, const std::vector<id> & ids,
 *             const std::vector<datum> & data);
 *
 * The example pointer may be null; it merely needs to be of the
 * correct type.  It's just here because function overloading in C++
 * is easy, whereas SFINAE is hard and partial template specialization
 * of functions is impossible.
 *
 * No guarantee about operation ordering is made - this function will
 * attempt to act on data in the order in which it is received.
 *
 * All receives and actions are completed before this function
 * returns.
 *
 * Not all sends may have yet completed.  The supplied container of
 * Request objects, \p req, has more requests inserted, one for each
 * of the data sends.  These requests must be waited on before the \p
 * data map is deleted.
 */
template <typename datum,
          typename MapToVectors,
          typename RequestContainer,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               RequestContainer & reqs,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const datum * example);

/**
 * Send query vectors, receive and answer them with vectors of data,
 * then act on those answers.
 *
 * The \p data map is indexed by processor ids as keys, and for each
 * processor id in the map there should be a vector of query ids to send.
 *
 * Query data which is received from other processors will be operated
 * on by
 * gather_data(processor_id_type pid, const std::vector<id> & ids,
 *             std::vector<datum> & data)
 *
 * Answer data which is received from other processors will be operated on by
 * act_on_data(processor_id_type pid, const std::vector<id> & ids,
 *             const std::vector<datum> & data);
 *
 * The example pointer may be null; it merely needs to be of the
 * correct type.  It's just here because function overloading in C++
 * is easy, whereas SFINAE is hard and partial template specialization
 * of functions is impossible.
 *
 * No guarantee about operation ordering is made - this function will
 * attempt to act on data in the order in which it is received.
 *
 * All communication and actions are complete when this function
 * returns.
 */
template <typename datum,
          typename MapToVectors,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const datum * example);

//------------------------------------------------------------------------
// Parallel function overloads
//

/*
 * A specialization for types that are harder to non-blocking receive.
 */
template <template <typename, typename, typename ...> class MapType,
          typename KeyType,
          typename ValueType,
          typename A1,
          typename A2,
          typename ... ExtraTypes,
          typename RequestContainer,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapType<processor_id_type, std::vector<std::vector<ValueType,A1>,A2>, ExtraTypes...> & data,
                               RequestContainer & reqs,
                               ActionFunctor & act_on_data);



/*
 * A specialization for types that are harder to non-blocking receive.
 */
template <template <typename, typename, typename ...> class MapType,
          typename KeyType,
          typename ValueType,
          typename A1,
          typename A2,
          typename ... ExtraTypes,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapType<processor_id_type, std::vector<std::vector<ValueType,A1>,A2>, ExtraTypes...> & data,
                               ActionFunctor & act_on_data);

/*
 * A specialization for types that are harder to non-blocking receive.
 */
template <typename datum,
          typename A,
          typename MapToVectors,
          typename RequestContainer,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               RequestContainer & reqs,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const std::vector<datum,A> * example);







//------------------------------------------------------------------------
// Parallel members
//

template <typename MapToVectors,
          typename RequestContainer,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & data,
                               RequestContainer & reqs,
                               ActionFunctor & act_on_data)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  processor_id_type num_procs = comm.size();

  // Size of vectors to send to each procesor
  std::vector<std::size_t> will_send_to(num_procs, 0);
  processor_id_type num_sends = 0;
  for (auto & datapair : data)
    {
      // Don't try to send anywhere that doesn't exist
      libmesh_assert_less(datapair.first, num_procs);

      // Don't give us empty vectors to send
      libmesh_assert_greater(datapair.second.size(), 0);

      will_send_to[datapair.first] = datapair.second.size();
      num_sends++;
    }

  // Tell everyone about where everyone will send to
  comm.alltoall(will_send_to);

  // will_send_to now represents who we'll receive from
  // give it a good name
  auto & will_receive_from = will_send_to;

  // This function only works for "flat" data that we can pre-size
  // receive buffers for: a map to vectors-of-standard-types, not e.g.
  // vectors-of-vectors.
  //
  // Trying to instantiate a StandardType<T> gives us a compiler error
  // where otherwise we would have had a runtime error.
  //
  // Creating a StandardType<T> manually also saves our APIs from
  // having to do a bunch of automatic creations later.
  //
  // This object will be free'd before all non-blocking communications
  // complete, but the MPI standard for MPI_Type_free specifies "Any
  // communication that is currently using this datatype will
  // complete normally." so we're cool.
  typedef decltype(data.begin()->second.front()) ref_type;
  typedef typename std::remove_reference<ref_type>::type nonref_type;
  StandardType<typename std::remove_const<nonref_type>::type> datatype;

  // We'll grab a tag so we can overlap request sends and receives
  // without confusing one for the other
  MessageTag tag = comm.get_unique_tag(1225);

  MapToVectors received_data;

  // Post all of the sends, non-blocking
  for (auto & datapair : data)
    {
      processor_id_type destid = datapair.first;
      libmesh_assert_less(destid, num_procs);
      auto & datum = datapair.second;

      // Just act on data if the user requested a send-to-self
      if (destid == comm.rank())
        act_on_data(destid, datum);
      else
        {
          Request sendreq;
          comm.send(destid, datum, datatype, sendreq, tag);
          reqs.insert(reqs.end(), sendreq);
        }
    }

  // Post all of the receives, non-blocking
  std::vector<Request> receive_reqs;
  std::vector<processor_id_type> receive_procids;
  for (processor_id_type proc_id = 0; proc_id < num_procs; proc_id++)
    if (will_receive_from[proc_id] && proc_id != comm.rank())
      {
        Request req;
        auto & incoming_data = received_data[proc_id];
        incoming_data.resize(will_receive_from[proc_id]);
        comm.receive(proc_id, incoming_data, datatype, req, tag);
        receive_reqs.push_back(req);
        receive_procids.push_back(proc_id);
      }

  while(receive_reqs.size())
    {
      std::size_t completed = waitany(receive_reqs);
      processor_id_type proc_id = receive_procids[completed];
      receive_reqs.erase(receive_reqs.begin() + completed);
      receive_procids.erase(receive_procids.begin() + completed);

      act_on_data(proc_id, received_data[proc_id]);
      received_data.erase(proc_id);
    }
}



template <template <typename, typename, typename ...> class MapType,
          typename ValueType,
          typename A1,
          typename A2,
          typename ... ExtraTypes,
          typename RequestContainer,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapType<processor_id_type, std::vector<std::vector<ValueType,A1>,A2>, ExtraTypes...> & data,
                               RequestContainer & reqs,
                               ActionFunctor & act_on_data)
{
  // This function must be run on all processors at once
  libmesh_parallel_only(comm);

  processor_id_type num_procs = comm.size();

  // Size of vectors to send to each procesor
  std::vector<std::size_t> will_send_to(num_procs, 0);
  processor_id_type num_sends = 0;
  for (auto & datapair : data)
    {
      // Don't try to send anywhere that doesn't exist
      libmesh_assert_less(datapair.first, num_procs);

      // Don't give us empty vectors to send
      libmesh_assert_greater(datapair.second.size(), 0);

      will_send_to[datapair.first] = datapair.second.size();
      num_sends++;
    }

  // Tell everyone about where everyone will send to
  comm.alltoall(will_send_to);

  // will_send_to now represents who we'll receive from
  // give it a good name
  auto & will_receive_from = will_send_to;

  processor_id_type n_receives = 0;
  for (processor_id_type proc_id = 0; proc_id < num_procs; proc_id++)
    if (will_receive_from[proc_id])
      n_receives++;

  // We'll construct a datatype once for repeated use
  StandardType<ValueType> datatype;

  // We'll grab a tag so we can overlap request sends and receives
  // without confusing one for the other
  MessageTag tag = comm.get_unique_tag(1225);

  // Post all of the sends, non-blocking
  for (auto & datapair : data)
    {
      processor_id_type destid = datapair.first;
      libmesh_assert_less(destid, num_procs);
      auto & datum = datapair.second;

      // Just act on data if the user requested a send-to-self
      if (destid == comm.rank())
        {
          act_on_data(destid, datum);
          n_receives--;
        }
      else
        {
          Request sendreq;
          comm.send(destid, datum, datatype, sendreq, tag);
          reqs.insert(reqs.end(), sendreq);
        }
    }

  // Post all of the receives.
  //
  // Use blocking API here since we can't use the pre-sized
  // non-blocking APIs with this data type.
  //
  // FIXME - implement Derek's API from #1684, switch to that!
  for (processor_id_type i = 0; i != n_receives; ++i)
    {
      Status stat(comm.probe(any_source, tag));
      const processor_id_type
        proc_id = cast_int<processor_id_type>(stat.source());

      std::vector<std::vector<ValueType,A1>,A2> received_data;
      comm.receive(proc_id, received_data, datatype, tag);
      act_on_data(proc_id, received_data);
    }
}



template <typename MapToVectors,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & data,
                               ActionFunctor & act_on_data)
{
  std::vector<Request> requests;

  push_parallel_vector_data(comm, data, requests, act_on_data);

  wait(requests);
}



template <template <typename, typename, typename ...> class MapType,
          typename ValueType,
          typename A1,
          typename A2,
          typename ... ExtraTypes,
          typename ActionFunctor>
void push_parallel_vector_data(const Communicator & comm,
                               const MapType<processor_id_type, std::vector<std::vector<ValueType,A1>,A2>, ExtraTypes...> & data,
                               ActionFunctor & act_on_data)
{
  std::vector<Request> requests;

  push_parallel_vector_data(comm, data, requests, act_on_data);

  wait(requests);
}



template <typename datum,
          typename MapToVectors,
          typename RequestContainer,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               RequestContainer & reqs,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const datum *)
{
  typedef typename MapToVectors::mapped_type query_type;

  std::map<processor_id_type, std::vector<datum> >
    response_data, received_data;
  std::vector<Request> response_reqs;

  StandardType<datum> datatype;

  // We'll grab a tag so we can overlap request sends and receives
  // without confusing one for the other
  MessageTag tag = comm.get_unique_tag(105);

  auto gather_functor =
    [&comm, &gather_data, &response_data, &response_reqs, &datatype, &tag]
    (processor_id_type pid, query_type query)
    {
      gather_data(pid, query, response_data[pid]);
      libmesh_assert_equal_to(query.size(), response_data[pid].size());

      // Just act on data later if the user requested a send-to-self
      if (pid != comm.rank())
        {
          Request sendreq;
          comm.send(pid, response_data[pid], datatype, sendreq, tag);
          response_reqs.push_back(sendreq);
        }
    };

  push_parallel_vector_data (comm, queries, reqs, gather_functor);

  // Every outgoing query should now have an incoming response.
  // Post all of the receives, non-blocking
  std::vector<Request> receive_reqs;
  std::vector<processor_id_type> receive_procids;
  for (auto & querypair : queries)
    {
      processor_id_type proc_id = querypair.first;
      libmesh_assert_less(proc_id, comm.size());

      if (proc_id == comm.rank())
        {
          libmesh_assert(queries.count(proc_id));
          libmesh_assert_equal_to(queries.at(proc_id).size(),
                                  response_data.at(proc_id).size());
          act_on_data(proc_id, queries.at(proc_id), response_data.at(proc_id));
        }
      else
        {
          auto & querydata = querypair.second;
          Request req;
          auto & incoming_data = received_data[proc_id];
          incoming_data.resize(querydata.size());
          comm.receive(proc_id, incoming_data, datatype, req, tag);
          receive_reqs.push_back(req);
          receive_procids.push_back(proc_id);
        }
    }

  while(receive_reqs.size())
    {
      std::size_t completed = waitany(receive_reqs);
      processor_id_type proc_id = receive_procids[completed];
      receive_reqs.erase(receive_reqs.begin() + completed);
      receive_procids.erase(receive_procids.begin() + completed);

      libmesh_assert(queries.count(proc_id));
      libmesh_assert_equal_to(queries.at(proc_id).size(),
                              received_data[proc_id].size());
      act_on_data(proc_id, queries.at(proc_id), received_data[proc_id]);
      received_data.erase(proc_id);
    }

  wait(response_reqs);
}


template <typename datum,
          typename MapToVectors,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const datum * example)
{
  std::vector<Request> requests;

  pull_parallel_vector_data(comm, queries, requests, gather_data,
                            act_on_data, example);

  wait(requests);
}


template <typename datum,
          typename A,
          typename MapToVectors,
          typename RequestContainer,
          typename GatherFunctor,
          typename ActionFunctor>
void pull_parallel_vector_data(const Communicator & comm,
                               const MapToVectors & queries,
                               RequestContainer & reqs,
                               GatherFunctor & gather_data,
                               ActionFunctor & act_on_data,
                               const std::vector<datum,A> *)
{
  typedef typename MapToVectors::mapped_type query_type;

  std::map<processor_id_type, std::vector<std::vector<datum,A>>>
    response_data;
  std::vector<Request> response_reqs;

  // We'll grab a tag so we can overlap request sends and receives
  // without confusing one for the other
  MessageTag tag = comm.get_unique_tag(105);

  auto gather_functor =
    [&comm, &gather_data, &act_on_data,
     &response_data, &response_reqs, &tag]
    (processor_id_type pid, query_type query)
    {
      gather_data(pid, query, response_data[pid]);
      libmesh_assert_equal_to(query.size(),
                              response_data[pid].size());

      // Just act on data if the user requested a send-to-self
      if (pid == comm.rank())
        {
          act_on_data(pid, query, response_data[pid]);
        }
      else
        {
          Request sendreq;
          comm.send(pid, response_data[pid], sendreq, tag);
          response_reqs.push_back(sendreq);
        }
    };

  push_parallel_vector_data (comm, queries, reqs, gather_functor);

  // Every outgoing query should now have an incoming response.
  //
  // Post all of the receives.
  //
  // Use blocking API here since we can't use the pre-sized
  // non-blocking APIs with this data type.
  //
  // FIXME - implement Derek's API from #1684, switch to that!
  std::vector<Request> receive_reqs;
  std::vector<processor_id_type> receive_procids;
  for (std::size_t i = 0,
       n_queries = queries.size() - queries.count(comm.rank());
       i != n_queries; ++i)
    {
      Status stat(comm.probe(any_source, tag));
      const processor_id_type
        proc_id = cast_int<processor_id_type>(stat.source());

      std::vector<std::vector<datum,A>> received_data;
      comm.receive(proc_id, received_data, tag);

      libmesh_assert(queries.count(proc_id));
      auto & querydata = queries.at(proc_id);
      libmesh_assert_equal_to(querydata.size(), received_data.size());
      act_on_data(proc_id, querydata, received_data);
    }

  wait(response_reqs);
}


} // namespace Parallel


} // namespace libMesh

#endif // LIBMESH_PARALLEL_SYNC_H
