// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <sstream>
#include <fstream>

#include <unistd.h> // for getpid()

// Local includes
#include "libmesh/xdr_cxx.h"
#include "libmesh/libmesh_logging.h"
#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h"
#endif


// Anonymous namespace for implementation details.
namespace {

// Nasty hacks for reading/writing zipped files
void bzip_file (const std::string & unzipped_name)
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

std::string unzip_file (const std::string & name)
{
  std::ostringstream pid_suffix;
  pid_suffix << '_' << getpid();

  std::string new_name = name;
  if (name.size() - name.rfind(".bz2") == 4)
    {
#ifdef LIBMESH_HAVE_BZIP
      new_name.erase(new_name.end() - 4, new_name.end());
      new_name += pid_suffix.str();
      LOG_SCOPE("system(bunzip2)", "XdrIO");
      std::string system_string = "bunzip2 -f -k -c ";
      system_string += name + " > " + new_name;
      if (std::system(system_string.c_str()))
        libmesh_file_error(system_string);
#else
      libmesh_error_msg("ERROR: need bzip2/bunzip2 to open .bz2 file " << name);
#endif
    }
  else if (name.size() - name.rfind(".xz") == 3)
    {
#ifdef LIBMESH_HAVE_XZ
      new_name.erase(new_name.end() - 3, new_name.end());
      new_name += pid_suffix.str();
      LOG_SCOPE("system(xz -d)", "XdrIO");
      std::string system_string = "xz -f -d -k -c ";
      system_string += name + " > " + new_name;
      if (std::system(system_string.c_str()))
        libmesh_file_error(system_string);
#else
      libmesh_error_msg("ERROR: need xz to open .xz file " << name);
#endif
    }
  return new_name;
}

void xzip_file (const std::string & unzipped_name)
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
void remove_unzipped_file (const std::string & name)
{
  std::ostringstream pid_suffix;
  pid_suffix << '_' << getpid();

  // If we temporarily decompressed a file, remove the
  // uncompressed version
  if (name.size() - name.rfind(".bz2") == 4)
    {
      std::string new_name(name.begin(), name.end()-4);
      new_name += pid_suffix.str();
      std::remove(new_name.c_str());
    }
  if (name.size() - name.rfind(".xz") == 3)
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
Xdr::Xdr (const std::string & name,
          const XdrMODE m) :
  mode(m),
  file_name(name),
#ifdef LIBMESH_HAVE_XDR
  fp(libmesh_nullptr),
#endif
  in(),
  out(),
  comm_len(xdr_MAX_STRING_LENGTH),
  gzipped_file(false),
  bzipped_file(false),
  xzipped_file(false)
{
  this->open(name);
}



Xdr::~Xdr()
{
  this->close();
}



void Xdr::open (const std::string & name)
{
  file_name = name;

  if (name == "")
    return;

  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        fp = fopen(name.c_str(), (mode == ENCODE) ? "w" : "r");
        if (!fp)
          libmesh_file_error(name.c_str());
        xdrs.reset(new XDR);
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
        gzipped_file = (name.size() - name.rfind(".gz")  == 3);
        bzipped_file = (name.size() - name.rfind(".bz2") == 4);
        xzipped_file = (name.size() - name.rfind(".xz") == 3);

        if (gzipped_file)
          {
#ifdef LIBMESH_HAVE_GZSTREAM
            igzstream * inf = new igzstream;
            libmesh_assert(inf);
            in.reset(inf);
            inf->open(name.c_str(), std::ios::in);
#else
            libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
          }
        else
          {
            std::ifstream * inf = new std::ifstream;
            libmesh_assert(inf);
            in.reset(inf);

            std::string new_name = unzip_file(name);

            inf->open(new_name.c_str(), std::ios::in);
          }

        libmesh_assert(in.get());

        if (!in->good())
          libmesh_file_error(name);
        return;
      }

    case WRITE:
      {
        gzipped_file = (name.size() - name.rfind(".gz")  == 3);
        bzipped_file = (name.size() - name.rfind(".bz2") == 4);
        xzipped_file = (name.size() - name.rfind(".xz")  == 3);

        if (gzipped_file)
          {
#ifdef LIBMESH_HAVE_GZSTREAM
            ogzstream * outf = new ogzstream;
            libmesh_assert(outf);
            out.reset(outf);
            outf->open(name.c_str(), std::ios::out);
#else
            libmesh_error_msg("ERROR: need gzstream to handle .gz files!!!");
#endif
          }
        else
          {
            std::ofstream * outf = new std::ofstream;
            libmesh_assert(outf);
            out.reset(outf);

            std::string new_name = name;

            if (bzipped_file)
              new_name.erase(new_name.end() - 4, new_name.end());

            if (xzipped_file)
              new_name.erase(new_name.end() - 3, new_name.end());

            outf->open(new_name.c_str(), std::ios::out);
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
            fp = libmesh_nullptr;
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
        if (in.get() != libmesh_nullptr)
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
        if (out.get() != libmesh_nullptr)
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
        if (in.get() != libmesh_nullptr)
          return in->good();
        return false;
      }

    case WRITE:
      {
        if (out.get() != libmesh_nullptr)
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
      return xdr_vector(x, (char *) &a[0], length, sizeof(T),
                        xdr_translator<T>());
    }
  else
    return true;
}

template <typename T>
bool xdr_translate(XDR * x, std::vector<std::complex<T> > & a)
{
  unsigned int length = cast_int<unsigned int>(a.size());
  bool b = xdr_u_int(x, &length);
  a.resize(length);
  typename std::vector<std::complex<T> >::iterator iter = a.begin();
  for (; iter != a.end(); ++iter)
    if (!xdr_translate(x, *iter))
      b = false;
  return b;
}

template <>
bool xdr_translate(XDR * x, std::vector<std::string> & s)
{
  unsigned int length = cast_int<unsigned int>(s.size());
  bool b = xdr_u_int(x, &length);
  s.resize(length);
  std::vector<std::string>::iterator iter = s.begin();
  for (; iter != s.end(); ++iter)
    if (!xdr_translate(x, *iter))
      b = false;
  return b;
}

template <>
xdrproc_t xdr_translator<int>() { return (xdrproc_t)(xdr_int); }

template <>
xdrproc_t xdr_translator<unsigned int>() { return (xdrproc_t)(xdr_u_int); }

template <>
xdrproc_t xdr_translator<long int>() { return (xdrproc_t)(xdr_long); }

template <>
xdrproc_t xdr_translator<unsigned long int>() { return (xdrproc_t)(xdr_u_long); }

template <>
xdrproc_t xdr_translator<unsigned long long>() { return (xdrproc_t)(xdr_u_longlong_t); }

template <>
xdrproc_t xdr_translator<short int>() { return (xdrproc_t)(xdr_short); }

template <>
xdrproc_t xdr_translator<unsigned short int>() { return (xdrproc_t)(xdr_u_short); }

template <>
xdrproc_t xdr_translator<char>() { return (xdrproc_t)(xdr_char); }

template <>
xdrproc_t xdr_translator<signed char>() { return (xdrproc_t)(xdr_char); }

template <>
xdrproc_t xdr_translator<unsigned char>() { return (xdrproc_t)(xdr_u_char); }

template <>
xdrproc_t xdr_translator<float>() { return (xdrproc_t)(xdr_float); }

template <>
xdrproc_t xdr_translator<double>() { return (xdrproc_t)(xdr_double); }

// FIXME - no XDR love for long doubles?
template <>
xdrproc_t xdr_translator<long double>() { return (xdrproc_t)(xdr_double); }

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

  for (unsigned int c=0; c<std::strlen(comm); c++)
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

  for (std::size_t i=0; i<a.size(); i++)
    {
      libmesh_assert(in.get());
      libmesh_assert (in->good());
      *in >> a[i];
    }
  in->getline(comm, comm_len);
}

template <typename T>
void Xdr::do_read(std::vector<std::complex<T> > & a)
{
  unsigned int length=0;
  data(length, "# vector length x 2 (complex)");
  a.resize(length);

  for (std::size_t i=0; i<a.size(); i++)
    {
      T r, im;
      libmesh_assert(in.get());
      libmesh_assert (in->good());
      *in >> r >> im;
      a[i] = std::complex<T>(r,im);
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

  for (std::size_t i=0; i<a.size(); i++)
    {
      libmesh_assert(out.get());
      libmesh_assert (out->good());
      this->do_write(a[i]);
      *out << "\t ";
    }
}

template <typename T>
void Xdr::do_write(std::vector<std::complex<T> > & a)
{
  std::size_t length=a.size();
  data(length, "# vector length x 2 (complex)");

  for (std::size_t i=0; i<a.size(); i++)
    {
      libmesh_assert(out.get());
      libmesh_assert (out->good());
      this->do_write(a[i]);
      *out << "\t ";
    }
}



template <typename T>
void Xdr::data (T & a, const char * comment_in)
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

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

        this->do_write(a);

        // If there's a comment provided, write a tab character and
        // then the comment.
        if (std::string(comment_in) != "")
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

        if (size_of_type <= 4) // 32-bit types
          {
            xdr_vector(xdrs.get(),
                       (char *) val,
                       len,
                       size_of_type,
                       (xdrproc_t) xdr_u_int);
          }
        else // 64-bit types
          {
            xdr_vector(xdrs.get(),
                       (char *) val,
                       len,
                       size_of_type,
                       (xdrproc_t) xdr_u_hyper);
          }

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

        if (size_of_type <= 4) // 32-bit types
          {
            if (len > 0)
              xdr_vector(xdrs.get(),
                         (char *) val,
                         len,
                         size_of_type,
                         (xdrproc_t) xdr_u_int);
          }
        else // 64-bit types
          {
            if (len > 0)
              xdr_vector(xdrs.get(),
                         (char *) val,
                         len,
                         size_of_type,
                         (xdrproc_t) xdr_u_hyper);

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

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

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
                for (unsigned int i=0; i<imax; i++)
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
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (this->is_open());

        if (len > 0)
          xdr_vector(xdrs.get(),
                     (char *) val,
                     len,
                     sizeof(double),
                     (xdrproc_t) xdr_double);

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

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

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
                for (unsigned int i=0; i<imax; i++)
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
void Xdr::data_stream (float * val, const unsigned int len, const unsigned int line_break)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

        libmesh_assert (this->is_open());

        if (len > 0)
          xdr_vector(xdrs.get(),
                     (char *) val,
                     len,
                     sizeof(float),
                     (xdrproc_t) xdr_float);

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

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

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
                for (unsigned int i=0; i<imax; i++)
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
void Xdr::data_stream (long double * val, const unsigned int len, const unsigned int line_break)
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
        // if (len > 0)
        //   xdr_vector(xdrs.get(),
        //      (char *) val,
        //      len,
        //      sizeof(double),
        //      (xdrproc_t) xdr_quadruple);

        if (len > 0)
          {
            std::vector<double> io_buffer (len);

            // Fill io_buffer if we are writing.
            if (mode == ENCODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                io_buffer[cnt++] = val[i];

            xdr_vector(xdrs.get(),
                       (char *) &io_buffer[0],
                       len,
                       sizeof(double),
                       (xdrproc_t) xdr_double);

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

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

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
                for (unsigned int i=0; i<imax; i++)
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


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
void Xdr::data_stream (std::complex<double> * val, const unsigned int len, const unsigned int line_break)
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
            std::vector<double> io_buffer (2*len);

            // Fill io_buffer if we are writing.
            if (mode == ENCODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  io_buffer[cnt++] = val[i].real();
                  io_buffer[cnt++] = val[i].imag();
                }

            xdr_vector(xdrs.get(),
                       (char *) &io_buffer[0],
                       2*len,
                       sizeof(double),
                       (xdrproc_t) xdr_double);

            // Fill val array if we are reading.
            if (mode == DECODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  double re = io_buffer[cnt++];
                  double im = io_buffer[cnt++];
                  val[i] = std::complex<double>(re,im);
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
            double re, im;
            *in >> re >> im;
            val[i] = std::complex<double>(re,im);
          }

        return;
      }

    case WRITE:
      {
        libmesh_assert(out.get());
        libmesh_assert (out->good());

        // Save stream flags
        std::ios_base::fmtflags out_flags = out->flags();

        // We will use scientific notation with a precision of 16
        // digits in the following output.  The desired precision and
        // format will automatically determine the width.
        *out << std::scientific
             << std::setprecision(16);

        if (line_break == libMesh::invalid_uint)
          for (unsigned int i=0; i<len; i++)
            {
              libmesh_assert(out.get());
              libmesh_assert (out->good());
              *out << val[i].real() << ' ';
              *out << val[i].imag() << ' ';
            }
        else
          {
            const unsigned imax = std::min(line_break, len);
            unsigned int cnt=0;
            while (cnt < len)
              {
                for (unsigned int i=0; i<imax; i++)
                  {
                    libmesh_assert(out.get());
                    libmesh_assert (out->good());
                    *out << val[cnt].real() << ' ';
                    *out << val[cnt].imag();
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

template <>
void Xdr::data_stream (std::complex<long double> * val, const unsigned int len, const unsigned int line_break)
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
                       (char *) &io_buffer[0],
                       2*len,
                       sizeof(double),
                       (xdrproc_t) xdr_double);

            // Fill val array if we are reading.
            if (mode == DECODE)
              for (unsigned int i=0, cnt=0; i<len; i++)
                {
                  double re = io_buffer[cnt++];
                  double im = io_buffer[cnt++];
                  val[i] = std::complex<long double>(re, im);
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
            long double re, im;
            *in >> re >> im;
            val[i] = std::complex<long double>(re,im);
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
        // 'digits10' digits in the following output.  The desired
        // precision and format will automatically determine the
        // width.  Note: digit10 is the number of digits (in decimal
        // base) that can be represented without change.  Equivalent
        // to FLT_DIG, DBL_DIG or LDBL_DIG for floating types.
        *out << std::scientific
             << std::setprecision(std::numeric_limits<long double>::digits10);

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
                for (unsigned int i=0; i<imax; i++)
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
#endif // # LIBMESH_USE_COMPLEX_NUMBERS

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
template void Xdr::data<int>                              (int &,                             const char *);
template void Xdr::data<unsigned int>                     (unsigned int &,                    const char *);
template void Xdr::data<unsigned short int>               (unsigned short int &,              const char *);
template void Xdr::data<short int>                        (short int &,                       const char *);
template void Xdr::data<unsigned long int>                (unsigned long int &,               const char *);
template void Xdr::data<unsigned long long>               (unsigned long long &,              const char *);
template void Xdr::data<long int>                         (long int &,                        const char *);
template void Xdr::data<char>                             (char &,                            const char *);
template void Xdr::data<signed char>                      (signed char &,                     const char *);
template void Xdr::data<unsigned char>                    (unsigned char &,                   const char *);
template void Xdr::data<float>                            (float &,                           const char *);
template void Xdr::data<double>                           (double &,                          const char *);
template void Xdr::data<long double>                      (long double &,                     const char *);
template void Xdr::data<std::complex<float> >             (std::complex<float> &,             const char *);
template void Xdr::data<std::complex<double> >            (std::complex<double> &,            const char *);
template void Xdr::data<std::complex<long double> >       (std::complex<long double> &,       const char *);
template void Xdr::data<std::string>                      (std::string &,                     const char *);
template void Xdr::data<std::vector<int> >                (std::vector<int> &,                const char *);
template void Xdr::data<std::vector<unsigned int> >       (std::vector<unsigned int> &,       const char *);
template void Xdr::data<std::vector<short int> >          (std::vector<short int> &,          const char *);
template void Xdr::data<std::vector<unsigned short int> > (std::vector<unsigned short int> &, const char *);
template void Xdr::data<std::vector<long int> >           (std::vector<long int> &,           const char *);
template void Xdr::data<std::vector<unsigned long int> >  (std::vector<unsigned long int> &,  const char *);
template void Xdr::data<std::vector<unsigned long long> > (std::vector<unsigned long long> &, const char *);
template void Xdr::data<std::vector<char> >               (std::vector<char> &,               const char *);
template void Xdr::data<std::vector<signed char> >        (std::vector<signed char> &,        const char *);
template void Xdr::data<std::vector<unsigned char> >      (std::vector<unsigned char> &,      const char *);
template void Xdr::data<std::vector<float> >              (std::vector<float> &,              const char *);
template void Xdr::data<std::vector<double> >             (std::vector<double> &,             const char *);
template void Xdr::data<std::vector<long double> >        (std::vector<long double> &,        const char *);
template void Xdr::data<std::vector<std::complex<float> > >  (std::vector<std::complex<float> > &,  const char *);
template void Xdr::data<std::vector<std::complex<double> > > (std::vector<std::complex<double> > &, const char *);
template void Xdr::data<std::vector<std::complex<long double> > > (std::vector<std::complex<long double> > &, const char *);
template void Xdr::data<std::vector<std::string> >        (std::vector<std::string> &,        const char *);
template void Xdr::data_stream<int>                (int * val,                const unsigned int len, const unsigned int line_break);
template void Xdr::data_stream<unsigned short int> (unsigned short int * val, const unsigned int len, const unsigned int line_break);
template void Xdr::data_stream<unsigned int>       (unsigned int * val,       const unsigned int len, const unsigned int line_break);
template void Xdr::data_stream<unsigned long int>  (unsigned long int * val,  const unsigned int len, const unsigned int line_break);
template void Xdr::data_stream<unsigned long long> (unsigned long long * val, const unsigned int len, const unsigned int line_break);

} // namespace libMesh
