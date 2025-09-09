// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/thread_buffered_syncbuf.h"
#include "libmesh/threads.h"

namespace libMesh
{
ThreadBufferedSyncbuf::ThreadBufferedSyncbuf(std::streambuf & sink, bool flush_on_newline)
  : _sink(sink), _flush_on_newline(flush_on_newline), _mu(std::make_unique<Threads::spin_mutex>())
{
}

ThreadBufferedSyncbuf::~ThreadBufferedSyncbuf() = default;

ThreadBufferedSyncbuf::int_type ThreadBufferedSyncbuf::overflow(int_type ch)
{
  if (traits_type::eq_int_type(ch, traits_type::eof()))
    return traits_type::not_eof(ch);

  auto & t = this->thread_local_buffer();
  t.buf.push_back(traits_type::to_char_type(ch));

  if (_flush_on_newline && ch == '\n')
    this->emit_from_thread_local_buffer(t.buf, /*force_flush=*/false);

  return ch;
}

std::streamsize ThreadBufferedSyncbuf::xsputn(const char * s, std::streamsize n)
{
  auto & t = this->thread_local_buffer();
  t.buf.append(s, static_cast<size_t>(n));

  if (_flush_on_newline && !t.buf.empty() && t.buf.back() == '\n')
    this->emit_from_thread_local_buffer(t.buf, /*force_flush=*/false);

  return n;
}

int ThreadBufferedSyncbuf::sync()
{
  this->emit_from_thread_local_buffer(this->thread_local_buffer().buf, /*force_flush=*/true);
  return 0;
}

ThreadBufferedSyncbuf::ThreadLocalBuffer & ThreadBufferedSyncbuf::thread_local_buffer()
{
  static thread_local ThreadLocalBuffer t(*this);
  return t;
}

void ThreadBufferedSyncbuf::emit_from_thread_local_buffer(std::string & b, bool force_flush)
{
  if (b.empty())
    {
      if (force_flush)
        {
          Threads::spin_mutex::scoped_lock lk(*_mu);
          _sink.pubsync();
        }
      return;
    }

  {
    Threads::spin_mutex::scoped_lock lk(*_mu);
    _sink.sputn(b.data(), static_cast<std::streamsize>(b.size()));
    if (force_flush)
      _sink.pubsync();
  }
  b.clear();
}

ThreadBufferedSyncbuf::ThreadLocalBuffer::~ThreadLocalBuffer()
{
  // Runs at thread exit: commit any leftovers for this thread.
  if (!buf.empty())
    _owner.emit_from_thread_local_buffer(buf, /*force_flush=*/false);
}

}
