#include <libmesh/libmesh.h>

/**
 * This class uses RAII to control redirecting the libMesh::err stream
 * to NULL and restoring it around some operation where we do not want
 * to see output to the screen.
 */
class StreamRedirector
{
public:

  /**
   * Constructor; saves the original libMesh::err streambuf.
   */
  StreamRedirector()
    : _errbuf(libMesh::err.rdbuf())
  {
    libMesh::err.rdbuf(nullptr);
  }

  /**
   * Destructor: restores the stream.
   */
  ~StreamRedirector()
  {
    libMesh::err.rdbuf(_errbuf);
  }

private:
  std::streambuf * _errbuf;
};
