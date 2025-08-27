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

#ifndef LIBMESH_THREAD_BUFFERED_SYNCBUF_H
#define LIBMESH_THREAD_BUFFERED_SYNCBUF_H

#include <streambuf>
#include <string>
#include <memory>

namespace libMesh
{
namespace Threads
{
class spin_mutex;
}

class ThreadBufferedSyncbuf : public std::streambuf
{
public:
  explicit ThreadBufferedSyncbuf(std::streambuf & sink, bool flush_on_newline = true);

  /**
   * Defaulted destructor defined in implementation so we don't need to include threads.h
   */
  ~ThreadBufferedSyncbuf();

protected:
  // Single-char path
  int_type overflow(int_type ch) override;

  // Bulk write path
  std::streamsize xsputn(const char * s, std::streamsize n) override;

  // Flush path (std::flush / std::endl / ostream::flush)
  int sync() override;

private:
  /**
   * A class that wraps a thread-local string. We make sure when we're destructing and our buffer is
   * non-empty to write to our main output sink
   */
  class ThreadLocalBuffer
  {
  public:
    explicit ThreadLocalBuffer(ThreadBufferedSyncbuf & owner) : _owner(owner) {}

    /**
     * Ensures we write to our sink upon destruction if our buf is non-empty
     */
    ~ThreadLocalBuffer();

    /// per-thread buffer
    std::string buf;

  private:
    /// owning syncing stream buffer
    ThreadBufferedSyncbuf & _owner;
  };

  /**
   * One ThreadLocalBuffer instance per thread, constructed on first use with *this.
   */
  ThreadLocalBuffer & thread_local_buffer();

  /**
   * Emit from the thread local buffer to our wrapped \p _sink
   */
  void emit_from_thread_local_buffer(std::string & b, bool force_flush);

  /// Wrapped output sink
  std::streambuf & _sink;
  /// Whether to flush our sink to terminal/file on terminating new-line characters (\n)
  const bool _flush_on_newline;

  /// Serialization for commits to the shared sink
  std::unique_ptr<Threads::spin_mutex> _mu;
};

}
#endif // LIBMESH_THREAD_BUFFERED_SYNCBUF_H
