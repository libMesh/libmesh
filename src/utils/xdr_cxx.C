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


// C/C++ includes
#include <cstring>
#include <limits>
#include <iomanip>
#include <memory>
#include <sstream>
#include <fstream>

// Local includes
#include "libmesh/xdr_cxx.h"
#include "libmesh/libmesh_logging.h"
#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h" // For reading/writing compressed streams
#endif
#include "libmesh/utility.h" // unzip_file

#ifdef LIBMESH_HAVE_UNISTD_H
#include <unistd.h> // for getpid() on Unix
#endif
#ifdef LIBMESH_HAVE_PROCESS_H
#include <process.h> // for getpid() on Windows
#endif

// Anonymous namespace for implementation details.
namespace {

// Nasty hacks for reading/writing zipped files
void bzip_file (std::string_view unzipped_name)
{
#ifdef LIBMESH_HAVE_BZIP
  LOG_SCOPE("system(bzip2)", "XdrIO");

  std::string system_string = "bzip2 -f ";
  system_string += unzipped_name;
  if (std::system(system_string.c_str()))
    libmesh_file_error(system_string);
#else
  libmesh_error_msg("ERROR: need bzip2/bunzip2 to create " << unzipped_name << ".bz2");
#endif
}

void xzip_file (std::string_view unzipped_name)
{
#ifdef LIBMESH_HAVE_XZ
  LOG_SCOPE("system(xz)", "XdrIO");

  std::string system_string = "xz -f ";
  system_string += unzipped_name;
  if (std::system(system_string.c_str()))
    libmesh_file_error(system_string);
#else
  libmesh_error_msg("ERROR: need xz to create " << unzipped_name << ".xz");
#endif
}


// remove an unzipped file
void remove_unzipped_file (std::string_view name)
{
  std::ostringstream pid_suffix;
  pid_suffix << '_' << getpid();

  // If we temporarily decompressed a file, remove the
  // uncompressed version
  if (libMesh::Utility::ends_with(name, ".bz2"))
    {
      std::string new_name(name.begin(), name.end()-4);
      new_name += pid_suffix.str();
      std::remove(new_name.c_str());
    }
  if (libMesh::Utility::ends_with(name, ".xz"))
    {
      std::string new_name(name.begin(), name.end()-3);
      new_name += pid_suffix.str();
      std::remove(new_name.c_str());
    }
}
}

namespace libMesh
{

//-------------------------------------------------------------
// Xdr class implementation
Xdr::Xdr (std::string name,
          const XdrMODE m) :
  mode(m),
  file_name(std::move(name)),
#ifdef LIBMESH_HAVE_XDR
  fp(nullptr),
#endif
  in(),
  out(),
  comm_len(xdr_MAX_STRING_LENGTH),
  gzipped_file(false),
  bzipped_file(false),
  xzipped_file(false),
  version_number(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION, LIBMESH_MINOR_VERSION, LIBMESH_MICRO_VERSION))
{
  this->open(file_name);
}



Xdr::Xdr (std::ostream & stream) :
  mode(WRITE),
  file_name(),
#ifdef LIBMESH_HAVE_XDR
  fp(nullptr),
#endif
  in(),
  out(std::make_unique<std::ostream>(stream.rdbuf())),
  comm_len(xdr_MAX_STRING_LENGTH),
  gzipped_file(false),
  bzipped_file(false),
  xzipped_file(false),
  version_number(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION, LIBMESH_MINOR_VERSION, LIBMESH_MICRO_VERSION))
{
}



Xdr::Xdr (std::istream & stream) :
  mode(READ),
  file_name(),
#ifdef LIBMESH_HAVE_XDR
  fp(nullptr),
#endif
  in(std::make_unique<std::istream>(stream.rdbuf())),
  out(),
  comm_len(xdr_MAX_STRING_LENGTH),
  gzipped_file(false),
  bzipped_file(false),
  xzipped_file(false),
  version_number(LIBMESH_VERSION_ID(LIBMESH_MAJOR_VERSION, LIBMESH_MINOR_VERSION, LIBMESH_MICRO_VERSION))
{
}



Xdr::~Xdr()
{
  this->close();
}



void Xdr::open (std::string name)
{
  file_name = std::move(name);

  if (file_name == "")
    return;

  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        fp = fopen(file_name.c_str(), (mode == ENCODE) ? "w" : "r");
        if (!fp)
          libmesh_file_error(file_name.c_str());
        xdrs = std::make_unique<XDR>();
        xdrstdio_create (xdrs.get(), fp, (mode == ENCODE) ? XDR_ENCODE : XDR_DECODE);
#else

        libmesh_error_msg("ERROR: Functionality is not available.\n" \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;

      }

    case READ:
      {
        gzipped_file = Utility::ends_with(file_name, ".gz");
        bzipped_file = Utility::ends_with(file_name, ".bz2");
        xzipped_file = Utility::ends_with(file_name, ".xz");

        if (gzipped_file)
          {
#ifdef LIBMESH_HAVE_GZSTREAM
            auto inf = std::make_unique<igzstream>();
            libmesh_assert(inf);
            inf->open(file_name.c_str(), std::ios::in);
            in = std::move(inf);
#else
            libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
          }
        else
          {
            auto inf = std::make_unique<std::ifstream>();
            libmesh_assert(inf);

            std::string new_name = Utility::unzip_file(file_name);

            inf->open(new_name.c_str(), std::ios::in);
            in = std::move(inf);
          }

        libmesh_assert(in.get());

        if (!in->good())
          libmesh_file_error(file_name);
        return;
      }

    case WRITE:
      {
        gzipped_file = (file_name.rfind(".gz") == file_name.size() - 3);
        bzipped_file = (file_name.rfind(".bz2") == file_name.size() - 4);
        xzipped_file = (file_name.rfind(".xz") == file_name.size() - 3);

        if (gzipped_file)
          {
#ifdef LIBMESH_HAVE_GZSTREAM
            auto outf = std::make_unique<ogzstream>();
            libmesh_assert(outf);
            outf->open(file_name.c_str(), std::ios::out);
            out = std::move(outf);
#else
            libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
          }
        else
          {
            auto outf = std::make_unique<std::ofstream>();
            libmesh_assert(outf);

            std::string new_name = file_name;

            if (bzipped_file)
              new_name.erase(new_name.end() - 4, new_name.end());

            if (xzipped_file)
              new_name.erase(new_name.end() - 3, new_name.end());

            outf->open(new_name.c_str(), std::ios::out);
            out = std::move(outf);
          }

        libmesh_assert(out.get());
        libmesh_assert (out->good());
        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}



void Xdr::close ()
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        if (xdrs)
          {
            xdr_destroy (xdrs.get());
            xdrs.reset();
          }

        if (fp)
          {
            fflush(fp);
            fclose(fp);
            fp = nullptr;
          }
#else

        libmesh_error_msg("ERROR: Functionality is not available.\n" \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        file_name = "";
        return;
      }

    case READ:
      {
        if (in.get() != nullptr)
          {
            in.reset();

            if (bzipped_file || xzipped_file)
              remove_unzipped_file(file_name);
          }
        file_name = "";
        return;
      }

    case WRITE:
      {
        if (out.get() != nullptr)
          {
            out.reset();

            if (bzipped_file)
              bzip_file(std::string(file_name.begin(), file_name.end()-4));

            else if (xzipped_file)
              xzip_file(std::string(file_name.begin(), file_name.end()-3));
          }
        file_name = "";
        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}



bool Xdr::is_open() const
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        if (fp)
          if (xdrs)
            return true;

        return false;

#else

        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

        return false;

#endif

      }

    case READ:
      {
        if (in.get() != nullptr)
          return in->good();
        return false;
      }

    case WRITE:
      {
        if (out.get() != nullptr)
          return out->good();
        return false;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }

  return false;
}



bool Xdr::is_eof()
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR
        libmesh_assert(fp);

        // Are we already at eof?
        if (feof(fp))
          return true;

        // Or about to reach eof?
        int next = fgetc(fp);
        if (next == EOF)
          {
            // We should *only* be at EOF, not otherwise broken
            libmesh_assert(feof(fp));
            libmesh_assert(!ferror(fp));

            // Reset the EOF indicator
            clearerr(fp);
            libmesh_assert(!ferror(fp));

            // We saw EOF
            return true;
          }

        // We didn't see EOF; restore whatever we did see.
        ungetc(next, fp);
        break;
#else

        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

        return false;

#endif

      }
    case READ:
      {
        libmesh_assert(in.get());

        // Are we already at eof?
        if (in->eof())
          return true;

        // Or about to reach eof?
        int next = in->peek();
        if (next == EOF)
          {
            // We should *only* be at EOF, not otherwise broken
            libmesh_assert(in->eof());
            libmesh_assert(!in->fail());

            // Reset the EOF indicator
            in->clear();
            libmesh_assert(in->good());

            // We saw EOF
            return true;
          }
        break;
      }
    default:
      libmesh_error();
    }

  return false;
}



#ifdef LIBMESH_HAVE_XDR

// Anonymous namespace for Xdr::data helper functions
namespace
{

template <typename T>
xdrproc_t xdr_translator();

template <typename T>
bool xdr_translate(XDR * x, T & a)
{
  return (xdr_translator<T>())(x, &a, 0);
}

template <>
bool xdr_translate(XDR * x,
                   std::string & s)
{
  char sptr[xdr_MAX_STRING_LENGTH+1];
  std::copy(s.begin(), s.end(), sptr);
  sptr[s.size()] = 0;
  unsigned int length = xdr_MAX_STRING_LENGTH;

  // Get a pointer to the beginning of the buffer.  We need to pass
  // its address to the xdr API.
  char * begin = sptr;
  bool b = xdr_string(x, &begin, length);

  // This is necessary when reading, but inefficient when writing...
  length = cast_int<unsigned int>(std::strlen(sptr));
  s.resize(length);
  std::copy(sptr, sptr+length, s.begin());

  return b;
}

template <typename T>
bool xdr_translate(XDR * x, std::complex<T> & a)
{
  T r = a.real(), i = a.imag();
  bool b1 = xdr_translate(x, r);
  bool b2 = xdr_translate(x, i);
  a = std::complex<T>(r,i);
  return (b1 && b2);
}

template <typename T>
bool xdr_translate(XDR * x, std::vector<T> & a)
{
  unsigned int length = cast_int<unsigned int>(a.size());
  xdr_u_int(x, &length);
  if (length > 0)
    {
      a.resize(length);
      return xdr_vector(x, reinterpret_cast<char *>(a.data()), length, sizeof(T),
                        xdr_translator<T>());
    }
  else
    return true;
}

template <typename T>
bool xdr_translate(XDR * x, std::vector<std::complex<T>> & a)
{
  unsigned int length = cast_int<unsigned int>(a.size());
  bool b = xdr_u_int(x, &length);
  a.resize(length);
  for (auto & val : a)
    if (!xdr_translate(x, val))
      b = false;
  return b;
}

template <>
bool xdr_translate(XDR * x, std::vector<std::string> & s)
{
  unsigned int length = cast_int<unsigned int>(s.size());
  bool b = xdr_u_int(x, &length);
  s.resize(length);
  for (auto & val : s)
    if (!xdr_translate(x, val))
      b = false;
  return b;
}


// https://lists.libguestfs.org/archives/list/guestfs@lists.libguestfs.org/message/JNRVVA5VOIDBBQHTF7QYRSWZ2AQG4P7U/
//
//  With Mac OS X 10.9, xdrproc_t is no longer defined as:
//  typedef bool_t (*xdrproc_t)(XDR *, ...);
//  but instead as:
//  typdef bool_t (*xdrproc_t)(XDR *, void *, unsigned int);
//  For reference, Linux systems typically define it as:
//  typedef bool_t (*xdrproc_t)(XDR *, void *, ...);
//
// I'm not sure why it took over a decade later, but we're starting to
// get "converts to incompatible function type" errors on some OSX
// builds when trying to pass XDR's own functions to xdrproc_t
// arguments - apparently because the xdr_foo functions now only take
// two arguments?  We'll add some shims to convert them.
#ifdef __APPLE__
  #define libmesh_define_xdr(foo) \
    bool_t libmesh_xdr_##foo(XDR * x, void * v, unsigned int ui) { return xdr_##foo(x, v); }
#else
  #define libmesh_define_xdr(foo) \
    const xdrproc_t libmesh_xdr_##foo = (xdrproc_t)xdr_##foo;
#endif

libmesh_define_xdr(char)
libmesh_define_xdr(short)
libmesh_define_xdr(int)
libmesh_define_xdr(long)
libmesh_define_xdr(longlong_t)
libmesh_define_xdr(u_char)
libmesh_define_xdr(u_short)
libmesh_define_xdr(u_int)
libmesh_define_xdr(u_long)
libmesh_define_xdr(u_longlong_t)
libmesh_define_xdr(float)
libmesh_define_xdr(double)


template <>
xdrproc_t xdr_translator<int>()
{
  // Don't let XDR truncate anything on systems where int is 64-bit;
  // xdr_int is hard-coded to write 32 bits.
  if (sizeof(int) <= 4)
    return (xdrproc_t)(libmesh_xdr_int);
  else if (sizeof(int) == sizeof(long long))
    return (xdrproc_t)(libmesh_xdr_longlong_t);
  else if (sizeof(int) == sizeof(long))
    return (xdrproc_t)(libmesh_xdr_long);
  else
    libmesh_error();
}

template <>
xdrproc_t xdr_translator<unsigned int>()
{
  // Don't let XDR truncate anything on systems where int is 64-bit;
  // xdr_u_int is hard-coded to write 32 bits.
  if (sizeof(unsigned int) <= 4)
    return (xdrproc_t)(libmesh_xdr_u_int);
  else if (sizeof(unsigned int) == sizeof(unsigned long))
    return (xdrproc_t)(libmesh_xdr_u_long);
  else if (sizeof(unsigned int) == sizeof(unsigned long long))
    return (xdrproc_t)(libmesh_xdr_u_longlong_t);
  else
    libmesh_error();
}

template <>
xdrproc_t xdr_translator<long int>()
{
  // Don't let XDR truncate anything on systems where long is 64-bit;
  // xdr_long is hard-coded to write 32 bits.
  if (sizeof(long int) <= 4)
    return (xdrproc_t)(libmesh_xdr_long);
  else if (sizeof(long int) == sizeof(long long))
    return (xdrproc_t)(libmesh_xdr_longlong_t);
  else
    libmesh_error();
}

template <>
xdrproc_t xdr_translator<unsigned long int>()
{
  // Don't let XDR truncate anything on systems where long is 64-bit;
  // xdr_u_long is hard-coded to write 32 bits.  This bit us HARD.
  if (sizeof(unsigned long int) <= 4)
    return (xdrproc_t)(libmesh_xdr_u_long);
  else if (sizeof(unsigned long int) == sizeof(unsigned long long))
    return (xdrproc_t)(libmesh_xdr_u_longlong_t);
  else
    libmesh_error();
}

// All the other XDR routines should be safe, unless
// 1. You're on a system with sizeof(short)==8 and you want to store
// n>2^32 in a short; this will never happen since we assume
// sizeof(short) may be as small as 2 bytes and we use at least int
// for anything larger.
// 2. You're on a system with sizeof(long long) > 8, and you want to
// store n>2^64 in one.  Welcome, 22nd century programmers; sorry we
// couldn't accommodate you.
template <>
xdrproc_t xdr_translator<long long>() { return (xdrproc_t)(libmesh_xdr_longlong_t); }

template <>
xdrproc_t xdr_translator<unsigned long long>() { return (xdrproc_t)(libmesh_xdr_u_longlong_t); }

template <>
xdrproc_t xdr_translator<short int>() { return (xdrproc_t)(libmesh_xdr_short); }

template <>
xdrproc_t xdr_translator<unsigned short int>() { return (xdrproc_t)(libmesh_xdr_u_short); }

template <>
xdrproc_t xdr_translator<char>() { return (xdrproc_t)(libmesh_xdr_char); }

template <>
xdrproc_t xdr_translator<signed char>() { return (xdrproc_t)(libmesh_xdr_char); }

template <>
xdrproc_t xdr_translator<unsigned char>() { return (xdrproc_t)(libmesh_xdr_u_char); }

template <>
xdrproc_t xdr_translator<float>() { return (xdrproc_t)(libmesh_xdr_float); }

template <>
xdrproc_t xdr_translator<double>() { return (xdrproc_t)(libmesh_xdr_double); }

// FIXME - no XDR love for quadruple precision; not even for long double?
template <>
xdrproc_t xdr_translator<long double>() { return (xdrproc_t)(libmesh_xdr_double); }

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
template <>
xdrproc_t xdr_translator<Real>() { return (xdrproc_t)(libmesh_xdr_double); }
#endif

} // end anonymous namespace

#endif

template <typename T>
void Xdr::do_read(T & a)
{
  *in >> a;
  in->getline(comm, comm_len);
}

template <typename T>
void Xdr::do_read(std::complex<T> & a)
{
  T r, i;
  *in >> r >> i;
  a = std::complex<T>(r,i);
  in->getline(comm, comm_len);
}

template <>
void Xdr::do_read(std::string & a)
{
  in->getline(comm, comm_len);

  a = "";

  for (unsigned int c=0, sl=std::strlen(comm); c!=sl; c++)
    {
      if (comm[c] == '\t')
        break;
      a.push_back(comm[c]);
    }
}

template <typename T>
void Xdr::do_read(std::vector<T> & a)
{
  unsigned int length=0;
  data(length, "# vector length");
  a.resize(length);

  for (T & a_i : a)
    {
      libmesh_assert(in.get());
      libmesh_assert (in->good());
      *in >> a_i;
    }
  in->getline(comm, comm_len);
}

template <typename T>
void Xdr::do_read(std::vector<std::complex<T>> & a)
{
  unsigned int length=0;
  data(length, "# vector length x 2 (complex)");
  a.resize(length);

  for (std::complex<T> & a_i : a)
    {
      T r, im;
      libmesh_assert(in.get());
      libmesh_assert (in->good());
      *in >> r >> im;
      a_i = std::complex<T>(r,im);
    }
  in->getline(comm, comm_len);
}

template <typename T>
void Xdr::do_write(T & a) { *out << a; }

template <typename T>
void Xdr::do_write(std::complex<T> & a)
{
  *out << a.real() << "\t " << a.imag();
}

template <typename T>
void Xdr::do_write(std::vector<T> & a)
{
  std::size_t length = a.size();
  data(length, "# vector length");

  // Use scientific precision with lots of digits for the original type T.
  *out << std::scientific
       << std::setprecision(std::numeric_limits<T>::max_digits10);

  for (T & a_i : a)
    {
      libmesh_assert(out.get());
      libmesh_assert (out->good());
      this->do_write(a_i);
      *out << "\t ";
    }
}

template <typename T>
void Xdr::do_write(std::vector<std::complex<T>> & a)
{
  std::size_t length=a.size();
  data(length, "# vector length x 2 (complex)");

  // Use scientific precision with lots of digits for the original type T.
  *out << std::scientific
       << std::setprecision(std::numeric_limits<T>::max_digits10);

  for (std::complex<T> & a_i : a)
    {
      libmesh_assert(out.get());
      libmesh_assert (out->good());
      this->do_write(a_i);
      *out << "\t ";
    }
}



template <typename T>
void Xdr::data (T & a, std::string_view comment_in)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (is_open());

        xdr_translate(xdrs.get(), a);

#else

        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;
      }

    case READ:
      {
        libmesh_assert(in.get());
        libmesh_assert (in->good());

        this->do_read(a);

        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());

        // We will use scientific notation sufficient to exactly
        // represent our floating point precision in the following
        // output.  The desired precision and format will
        // automatically determine the width.
        *out << std::scientific
             << std::setprecision(std::numeric_limits<T>::max_digits10);

        this->do_write(a);

        // If there's a comment provided, write a tab character and
        // then the comment.
        if (comment_in != "")
          *out << "\t " << comment_in;

        // Go to the next line.
        *out << '\n';

        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}


template <typename T>
void Xdr::data_stream (T * val, const unsigned int len, const unsigned int line_break)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (this->is_open());

        unsigned int size_of_type = cast_int<unsigned int>(sizeof(T));

        xdr_vector(xdrs.get(),
                   (char *) val,
                   len,
                   size_of_type,
                   xdr_translator<T>());
#else
        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (this->is_open());

        unsigned int size_of_type = cast_int<unsigned int>(sizeof(T));

        if (len > 0)
          xdr_vector(xdrs.get(),
                     (char *) val,
                     len,
                     size_of_type,
                     xdr_translator<T>());
#else
        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;
      }

    case READ:
      {
        libmesh_assert(in.get());
        libmesh_assert (in->good());

        for (unsigned int i=0; i<len; i++)
          {
            libmesh_assert(in.get());
            libmesh_assert (in->good());
            *in >> val[i];
          }

        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());

        // We will use scientific notation sufficient to exactly
        // represent our floating point precision in the following
        // output.  The desired precision and format will
        // automatically determine the width.
        *out << std::scientific
             << std::setprecision(std::numeric_limits<T>::max_digits10);

        if (line_break == libMesh::invalid_uint)
          for (unsigned int i=0; i<len; i++)
            {
              libmesh_assert(out.get());
              libmesh_assert (out->good());
              *out << val[i] << " ";
            }
        else
          {
            const unsigned imax = std::min(line_break, len);
            unsigned int cnt=0;
            while (cnt < len)
              {
                for (unsigned int i=0; (i<imax && cnt<len); i++)
                  {
                    libmesh_assert(out.get());
                    libmesh_assert (out->good());
                    *out << val[cnt++];

                    // Write a space unless this is the last character on the current line.
                    if (i+1 != imax)
                      *out << " ";
                  }
                libmesh_assert(out.get());
                libmesh_assert (out->good());
                *out << '\n';
              }
          }

        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}



template <>
void Xdr::data_stream (double * val, const unsigned int len, const unsigned int line_break)
{
  this->_xfp_data_stream
    (val, len,
#ifdef LIBMESH_HAVE_XDR
     (xdrproc_t)libmesh_xdr_double,
#else
     nullptr,
#endif
     line_break, std::numeric_limits<double>::max_digits10);
}



template <>
void Xdr::data_stream (float * val, const unsigned int len, const unsigned int line_break)
{
  this->_xfp_data_stream
    (val, len,
#ifdef LIBMESH_HAVE_XDR
     (xdrproc_t)libmesh_xdr_float,
#else
     nullptr,
#endif
     line_break, std::numeric_limits<float>::max_digits10);
}



template <>
void Xdr::data_stream (long double * val, const unsigned int len, const unsigned int line_break)
{
  this->_xfp_data_stream
    (val, len, nullptr, line_break,
     std::numeric_limits<long double>::max_digits10);
}


#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
template <>
void Xdr::data_stream (Real * val, const unsigned int len, const unsigned int line_break)
{
  this->_xfp_data_stream(val, len, nullptr, line_break, 36);
}
#endif // LIBMESH_DEFAULT_QUADRUPLE_PRECISION



template <typename XFP>
void Xdr::_xfp_data_stream (XFP * val, const unsigned int len,
#ifdef LIBMESH_HAVE_XDR
                            xdrproc_t xdr_proc,
#else
                            void *,
#endif
                            const unsigned int line_break,
                            const int n_digits)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR
        libmesh_assert (this->is_open());

        if (len > 0)
          {
            if (xdr_proc)
              {
                xdr_vector(xdrs.get(),
                           (char *) val,
                           len,
                           sizeof(XFP),
                           xdr_proc);
                return;
              }

            // FIXME[JWP]: How to implement this for long double?  Mac
            // OS X defines 'xdr_quadruple' but AFAICT, it does not
            // exist for Linux... for now, reading/writing XDR files
            // with long doubles drops back to double precision, but
            // you can still write long double ASCII files of course.

            // FIXME[RHS]: 128 bit FP has the same problem as long
            // double, only much worse since even _Quad/__float128
            // aren't standard either.

            std::vector<double> io_buffer (len);

            // Fill io_buffer if we are writing.
            if (mode == ENCODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                io_buffer[cnt++] = double(val[i]);

            xdr_vector(xdrs.get(),
                       reinterpret_cast<char *>(io_buffer.data()),
                       len,
                       sizeof(double),
                       (xdrproc_t) libmesh_xdr_double);

            // Fill val array if we are reading.
            if (mode == DECODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  val[i] = io_buffer[cnt++];
                }
          }

#else

        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;
      }

    case READ:
      {
        libmesh_assert(in.get());
        libmesh_assert (in->good());

        for (unsigned int i=0; i<len; i++)
          {
            libmesh_assert(in.get());
            libmesh_assert (in->good());
            *in >> val[i];
          }

        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());

        // Save stream flags
        std::ios_base::fmtflags out_flags = out->flags();

        // We will use scientific notation with specified digit
        // count in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(n_digits);

        if (line_break == libMesh::invalid_uint)
          for (unsigned int i=0; i<len; i++)
            {
              libmesh_assert(out.get());
              libmesh_assert (out->good());
              *out << val[i] << ' ';
            }
        else
          {
            const unsigned imax = std::min(line_break, len);
            unsigned int cnt=0;
            while (cnt < len)
              {
                for (unsigned int i=0; (i<imax && cnt<len); i++)
                  {
                    libmesh_assert(out.get());
                    libmesh_assert (out->good());
                    *out << val[cnt++];

                    // Write a space unless this is the last character on the current line.
                    if (i+1 != imax)
                      *out << " ";
                  }
                libmesh_assert(out.get());
                libmesh_assert (out->good());
                *out << '\n';
              }
          }

        // Restore stream flags
        out->flags(out_flags);

        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}



template <>
void Xdr::data_stream (std::complex<double> * val, const unsigned int len, const unsigned int line_break)
{
  this->_complex_data_stream(val, len, line_break);
}



template <>
void Xdr::data_stream (std::complex<long double> * val, const unsigned int len, const unsigned int line_break)
{
  this->_complex_data_stream(val, len, line_break);
}



template <typename T>
void Xdr::_complex_data_stream (std::complex<T> * val, const unsigned int len, const unsigned int line_break)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (this->is_open());

        // FIXME[JWP]: How to implement this for long double?  Mac OS
        // X defines 'xdr_quadruple' but AFAICT, it does not exist for
        // Linux... for now, reading/writing XDR files with long
        // doubles drops back to double precision, but you can still
        // write long double ASCII files of course.

        if (len > 0)
          {
            std::vector<double> io_buffer (2*len);

            // Fill io_buffer if we are writing.
            if (mode == ENCODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  io_buffer[cnt++] = val[i].real();
                  io_buffer[cnt++] = val[i].imag();
                }

            xdr_vector(xdrs.get(),
                       reinterpret_cast<char *>(io_buffer.data()),
                       2*len,
                       sizeof(double),
                       (xdrproc_t) libmesh_xdr_double);

            // Fill val array if we are reading.
            if (mode == DECODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  double re = io_buffer[cnt++];
                  double im = io_buffer[cnt++];
                  val[i] = std::complex<T>(re, im);
                }
          }
#else

        libmesh_error_msg("ERROR: Functionality is not available.\n"    \
                          << "Make sure LIBMESH_HAVE_XDR is defined at build time\n" \
                          << "The XDR interface is not available in this installation");

#endif
        return;
      }

    case READ:
      {
        libmesh_assert(in.get());
        libmesh_assert (in->good());

        for (unsigned int i=0; i<len; i++)
          {
            libmesh_assert(in.get());
            libmesh_assert (in->good());
            T re, im;
            *in >> re >> im;
            val[i] = std::complex<T>(re,im);
          }

        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());


        // Save stream flags
        std::ios_base::fmtflags out_flags = out->flags();

        // We will use scientific notation with a precision of
        // 'max_digits10' digits in the following output.  The desired
        // precision and format will automatically determine the
        // width.  Note: digit10 is the number of digits (in decimal
        // base) that can be represented without change.  Equivalent
        // to FLT_DIG, DBL_DIG or LDBL_DIG for floating types.
        *out << std::scientific
             << std::setprecision(std::numeric_limits<T>::max_digits10);

        if (line_break == libMesh::invalid_uint)
          for (unsigned int i=0; i<len; i++)
            {
              libmesh_assert(out.get());
              libmesh_assert (out->good());
              *out << val[i].real() << ' ' << val[i].imag() << ' ';
            }
        else
          {
            const unsigned imax = std::min(line_break, len);
            unsigned int cnt=0;
            while (cnt < len)
              {
                for (unsigned int i=0; (i<imax && cnt<len); i++)
                  {
                    libmesh_assert(out.get());
                    libmesh_assert (out->good());
                    *out << val[cnt].real() << ' ' << val[cnt].imag();
                    cnt++;

                    // Write a space unless this is the last character on the current line.
                    if (i+1 != imax)
                      *out << " ";
                  }
                libmesh_assert(out.get());
                libmesh_assert (out->good());
                *out << '\n';
              }
          }

        // Restore stream flags
        out->flags(out_flags);

        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}



void Xdr::comment (std::string & comment_in)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
        return;
      }

    case READ:
      {
        libmesh_assert(in.get());
        libmesh_assert (in->good());
        in->getline(comm, comm_len);
        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());
        *out << "\t " << comment_in << '\n';
        return;
      }

    default:
      libmesh_error_msg("Invalid mode = " << mode);
    }
}


#undef xdr_REAL


//
template LIBMESH_EXPORT void Xdr::data<int>                              (int &,                             std::string_view);
template LIBMESH_EXPORT void Xdr::data<unsigned int>                     (unsigned int &,                    std::string_view);
template LIBMESH_EXPORT void Xdr::data<unsigned short int>               (unsigned short int &,              std::string_view);
template LIBMESH_EXPORT void Xdr::data<short int>                        (short int &,                       std::string_view);
template LIBMESH_EXPORT void Xdr::data<unsigned long int>                (unsigned long int &,               std::string_view);
template LIBMESH_EXPORT void Xdr::data<unsigned long long>               (unsigned long long &,              std::string_view);
template LIBMESH_EXPORT void Xdr::data<long int>                         (long int &,                        std::string_view);
template LIBMESH_EXPORT void Xdr::data<long long>                        (long long &,                       std::string_view);
template LIBMESH_EXPORT void Xdr::data<char>                             (char &,                            std::string_view);
template LIBMESH_EXPORT void Xdr::data<signed char>                      (signed char &,                     std::string_view);
template LIBMESH_EXPORT void Xdr::data<unsigned char>                    (unsigned char &,                   std::string_view);
template LIBMESH_EXPORT void Xdr::data<float>                            (float &,                           std::string_view);
template LIBMESH_EXPORT void Xdr::data<double>                           (double &,                          std::string_view);
template LIBMESH_EXPORT void Xdr::data<long double>                      (long double &,                     std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::complex<float>>              (std::complex<float> &,             std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::complex<double>>             (std::complex<double> &,            std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::complex<long double>>        (std::complex<long double> &,       std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::string>                      (std::string &,                     std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<int>>                 (std::vector<int> &,                std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<unsigned int>>        (std::vector<unsigned int> &,       std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<short int>>           (std::vector<short int> &,          std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<unsigned short int>>  (std::vector<unsigned short int> &, std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<long int>>            (std::vector<long int> &,           std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<long long>>           (std::vector<long long> &,          std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<unsigned long int>>   (std::vector<unsigned long int> &,  std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<unsigned long long>>  (std::vector<unsigned long long> &, std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<char>>                (std::vector<char> &,               std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<signed char>>         (std::vector<signed char> &,        std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<unsigned char>>       (std::vector<unsigned char> &,      std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<float>>               (std::vector<float> &,              std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<double>>              (std::vector<double> &,             std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<long double>>         (std::vector<long double> &,        std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<std::complex<float>>>  (std::vector<std::complex<float>> &,  std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<std::complex<double>>> (std::vector<std::complex<double>> &, std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<std::complex<long double>>> (std::vector<std::complex<long double>> &, std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<std::string>>        (std::vector<std::string> &,        std::string_view);
template LIBMESH_EXPORT void Xdr::data_stream<unsigned char>      (unsigned char * val,      const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<short int>          (short int * val,          const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<int>                (int * val,                const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<long long>          (long long * val,          const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<unsigned short int> (unsigned short int * val, const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<unsigned int>       (unsigned int * val,       const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<unsigned long int>  (unsigned long int * val,  const unsigned int len, const unsigned int line_break);
template LIBMESH_EXPORT void Xdr::data_stream<unsigned long long> (unsigned long long * val, const unsigned int len, const unsigned int line_break);

#ifdef LIBMESH_DEFAULT_QUADRUPLE_PRECISION
template LIBMESH_EXPORT void Xdr::data<Real>                             (Real &,                            std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::complex<Real>>               (std::complex<Real> &,              std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<Real>>                (std::vector<Real> &,               std::string_view);
template LIBMESH_EXPORT void Xdr::data<std::vector<std::complex<Real>>>  (std::vector<std::complex<Real>> &, std::string_view);
#endif

} // namespace libMesh
