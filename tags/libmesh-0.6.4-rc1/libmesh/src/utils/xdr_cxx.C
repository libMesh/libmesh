// "$Id$\n"

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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

// Local includes
#include "xdr_cxx.h"
#include "libmesh_logging.h"
#include "parallel.h"
#include "o_f_stream.h"
#ifdef LIBMESH_HAVE_GZSTREAM
# include "gzstream.h"
#endif


// Anonymous namespace for implementation details.
namespace {

  // Nasty hacks for reading/writing zipped files
  void zip_file (const std::string &unzipped_name)
  {
    // There's no parallel bzip2 for us to call
    libmesh_assert(libMesh::processor_id() == 0);

#ifdef LIBMESH_HAVE_BZIP
    START_LOG("system(bzip2)", "XdrIO");

    std::string system_string = "bzip2 -f ";
    system_string += unzipped_name;
    if (std::system(system_string.c_str()))
      libmesh_file_error(system_string);
     
    STOP_LOG("system(bzip2)", "XdrIO");
#else
    std::cerr << "ERROR: need bzip2/bunzip2 to handle .bz2 files!!!"
	      << std::endl;
    libmesh_error();
#endif
  }

  std::string unzip_file (const std::string &name)
  {
    // There's no parallel bunzip2 for us to call
    libmesh_assert(libMesh::processor_id() == 0);

    std::string new_name = name;
    if (name.size() - name.rfind(".bz2") == 4)
      {
	new_name.erase(new_name.end() - 4, new_name.end());
#ifdef LIBMESH_HAVE_BZIP
	START_LOG("system(bunzip2)", "XdrIO");
	std::string system_string = "bunzip2 -f -k ";
	system_string += name;
	if (std::system(system_string.c_str()))
          libmesh_file_error(system_string);
	STOP_LOG("system(bunzip2)", "XdrIO");
#else
	std::cerr << "ERROR: need bzip2/bunzip2 to handle .bz2 files!!!"
		  << std::endl;
	libmesh_error();
#endif
      }
    return new_name;
  }

  // remove an unzipped file
  void remove_unzipped_file (const std::string &name)
  {
    // If we temporarily decompressed a .bz2 file, remove the
    // uncompressed version
    if (name.size() - name.rfind(".bz2") == 4)
      {
	const std::string new_name(name.begin(), name.end()-4);
	std::remove(new_name.c_str());
      }
  }
}

//-------------------------------------------------------------
// Xdr class implementation
Xdr::Xdr (const std::string& name, const XdrMODE m) :
  mode(m),
  file_name(name),
#ifdef LIBMESH_HAVE_XDR
  xdrs(NULL),
  fp(NULL),
#endif
  in(NULL),
  out(NULL),
  comm_len(xdr_MAX_STRING_LENGTH),
  gzipped_file(false),
  bzipped_file(false)
{
  this->open(name);
}



Xdr::~Xdr()
{
  this->close();
}



void Xdr::open (const std::string& name)
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
	libmesh_assert (fp);
	xdrs = new XDR;
	xdrstdio_create (xdrs, fp, (mode == ENCODE) ? XDR_ENCODE : XDR_DECODE);
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();
	
#endif
	return;

      }

    case READ:
      {
	gzipped_file = (name.size() - name.rfind(".gz")  == 3);
	bzipped_file = (name.size() - name.rfind(".bz2") == 4);

	if (gzipped_file)
	  {
#ifdef LIBMESH_HAVE_GZSTREAM
	    igzstream *inf = new igzstream;
	    libmesh_assert (inf != NULL);
	    in.reset(inf);
	    inf->open(name.c_str(), std::ios::in);
#else
	    std::cerr << "ERROR: need gzstream to handle .gz files!!!"
		      << std::endl;
	    libmesh_error();
#endif
	  }
	else
	  {
	    std::ifstream *inf = new std::ifstream;
	    libmesh_assert (inf != NULL);
	    in.reset(inf);
	    
	    std::string new_name(bzipped_file ? unzip_file(name) : name);

	    inf->open(new_name.c_str(), std::ios::in);
	  }

	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	return;
      }

    case WRITE:
      {
	gzipped_file = (name.size() - name.rfind(".gz")  == 3);
	bzipped_file = (name.size() - name.rfind(".bz2") == 4);

	if (gzipped_file)
	  {
#ifdef LIBMESH_HAVE_GZSTREAM
	    ogzstream *outf = new ogzstream;
	    libmesh_assert (outf != NULL);
	    out.reset(outf);
	    outf->open(name.c_str(), std::ios::out);
#else
	    std::cerr << "ERROR: need gzstream to handle .gz files!!!"
		      << std::endl;
	    libmesh_error();
#endif
	  }
	else
	  {
	    std::ofstream *outf = new std::ofstream;
	    libmesh_assert (outf != NULL);
	    out.reset(outf);

	    std::string new_name = name;

	    if (bzipped_file)
	      new_name.erase(new_name.end() - 4, new_name.end());

	    outf->open(new_name.c_str(), std::ios::out);
	  }
	
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	return;
      }
      
    default:
      libmesh_error();
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
	    xdr_destroy (xdrs);
	    delete xdrs;
	    xdrs = NULL;
	  }
	
	if (fp)
	  {
	    fflush(fp);
	    fclose(fp);
	    fp = NULL;
	  }
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();
	
#endif
	file_name = "";
	return;
      }
      
    case READ:
      {
	if (in.get() != NULL)
	  {
	    in.reset();    
	    
	    if (bzipped_file)
	      remove_unzipped_file(file_name);
	  }
	file_name = "";
	return;
      }

    case WRITE:
      {
	if (out.get() != NULL)
	  {
	    out.reset();      

	    if (bzipped_file)
	      zip_file(std::string(file_name.begin(), file_name.end()-4));
	  }
	file_name = "";
	return;
      }

    default:
      libmesh_error();
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
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

	return false;	

#endif

      }
      
    case READ:
      {
	if (in.get() != NULL)
	  return in->good();
	return false;
      }

    case WRITE:
      {
	if (out.get() != NULL)
	  return out->good();
	return false;
      }

    default:
      libmesh_error();
    }

  return false;
}



void Xdr::data (int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	xdr_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	
	*out << a << "\t " << comment << '\n';

	return;
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (unsigned int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());

	xdr_u_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	xdr_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (unsigned short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	xdr_u_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (float& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	xdr_float(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (double& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	xdr_double(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	*in >> a; in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}



#ifdef LIBMESH_USE_COMPLEX_NUMBERS

void Xdr::data (std::complex<double>& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());
	double
	  _r=a.real(),
	  _i=a.imag();
	xdr_double(xdrs, &_r);
	xdr_double(xdrs, &_i);
	a = std::complex<double>(_r,_i);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	double _r, _i;
	*in >> _r;
	*in >> _i;
	a = std::complex<double>(_r,_i);
        in->getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << a.real() << "\t " 
	    << a.imag() << "\t "
	    << comment << '\n';
	
	return;
      }

    default:
      libmesh_error();
    }
}

#endif // LIBMESH_USE_COMPLEX_NUMBERS



void Xdr::data (std::vector<int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    *out << v[i] << " ";
	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (std::vector<unsigned int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    *out << v[i] << " ";
	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (std::vector<short int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(short int),
		   (xdrproc_t) xdr_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(short int),
		   (xdrproc_t) xdr_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    *out << v[i] << " ";
	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (std::vector<unsigned short int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned short int),
		   (xdrproc_t) xdr_u_short);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned short int),
		   (xdrproc_t) xdr_u_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    *out << v[i] << " ";
	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (std::vector<float>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(float),
		   (xdrproc_t) xdr_float);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(float),
		   (xdrproc_t) xdr_float);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    OFSRealscientific(*out,17,v[i]) << " ";
 	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



void Xdr::data (std::vector<double>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());

	unsigned int length = v.size();

	this->data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(double),
		   (xdrproc_t) xdr_double);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());

	unsigned int length=0;

	this->data(length, "# vector length");

	v.resize(length);

	// Note: GCC 3.4.1 will crash in debug mode here if length
	// is zero and you attempt to access the zeroth index of v.
	if (length > 0)
	  xdr_vector(xdrs, 
		     (char*) &v[0],
		     length,
		     sizeof(double),
		     (xdrproc_t) xdr_double);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	this->data(length, "# vector length");
	
	// If you were expecting to read in a vector at this
	// point, it's not going to happen if length == 0!
	// libmesh_assert (length != 0);
	
	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> v[i];
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    OFSRealscientific(*out,17,v[i]) << " ";
 	  }



	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}




#ifdef LIBMESH_USE_COMPLEX_NUMBERS

void Xdr::data (std::vector< std::complex<double> >& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length x 2 (complex)");

	std::vector< std::complex<double> >::iterator iter = v.begin();
	
	for (; iter != v.end(); ++iter)
	    data(*iter, "");

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	unsigned int length=0;

	data(length, "# vector length x 2 (complex)");

	v.resize(length);

	std::vector< std::complex<double> >::iterator iter = v.begin();
	
	for (; iter != v.end(); ++iter)
	    data(*iter, "");
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	unsigned int length=0;

	data(length, "# vector length x 2 (complex)");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	
	    double _r, _i;
	    *in >> _r;
	    *in >> _i;
	    v[i] = std::complex<double>(_r,_i);
	  }

	in->getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	unsigned int length=v.size();

	data(length, "# vector length x 2 (complex)");

 	for (unsigned int i=0; i<v.size(); i++)
 	  {
 	    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	    OFSNumberscientific(*out,17,v[i]) << " ";
 	  }

	*out << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}

#endif // ifdef LIBMESH_USE_COMPLEX_NUMBERS




void Xdr::data (std::string& s, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	{
	  char* sptr = new char[s.size()+1];

	  for (unsigned int c=0; c<s.size(); c++)
	    sptr[c] = s[c];
	
	  sptr[s.size()] = '\0';
	  
	  xdr_string(xdrs,
		     &sptr,
		     std::strlen(sptr));

	  delete [] sptr;
	}

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (is_open());

	{
	  char* sptr = new char[xdr_MAX_STRING_LENGTH];
	  
	  xdr_string(xdrs,
		     &sptr,
		     xdr_MAX_STRING_LENGTH);

	  s.resize(std::strlen(sptr));

	  for (unsigned int c=0; c<s.size(); c++)
	    s[c] = sptr[c];
	  
	  delete [] sptr;  
	}

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	in->getline(comm, comm_len);

//#ifndef LIBMESH_BROKEN_IOSTREAM
//	s.clear();
//#else
	s = "";
//#endif

	for (unsigned int c=0; c<std::strlen(comm); c++)
	  {
	    if (comm[c] == '\t') 
	      break;
	    
	    s.push_back(comm[c]);
	  }

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	*out << s << "\t " << comment << '\n';

	return;	
      }

    default:
      libmesh_error();
    }
}



template <typename T>
void Xdr::data_stream (T *val, const unsigned int len, const unsigned int line_break)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());


	xdr_vector(xdrs, 
		   (char*) val,
		   len,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());

	if (len > 0)
	  xdr_vector(xdrs, 
		     (char*) val,
		     len,
		     sizeof(unsigned int),
		     (xdrproc_t) xdr_u_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	for (unsigned int i=0; i<len; i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> val[i];
	  }

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	if (line_break == libMesh::invalid_uint)
	  for (unsigned int i=0; i<len; i++)
	    {
	      libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	      *out << val[i] << " ";
	    }
	else
	  {
	    unsigned int cnt=0;
	    while (cnt < len)
	      {
		for (unsigned int i=0; i<std::min(line_break,len); i++)
		  {
		    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		    *out << val[cnt++] << " ";
		  }
		libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		*out << '\n';
	      }
	  }

	return;	
      }

    default:
      libmesh_error();
    }
}



template <>
void Xdr::data_stream (double *val, const unsigned int len, const unsigned int line_break)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef LIBMESH_HAVE_XDR

	libmesh_assert (this->is_open());

	if (len > 0)
	  xdr_vector(xdrs, 
		     (char*) val,
		     len,
		     sizeof(double),
		     (xdrproc_t) xdr_double);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	for (unsigned int i=0; i<len; i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> val[i];
	  }

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	if (line_break == libMesh::invalid_uint)
	  for (unsigned int i=0; i<len; i++)
	    {
	      libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	      OFSRealscientific(*out,17,val[i]) << " ";
	    }
	else
	  {
	    unsigned int cnt=0;
	    while (cnt < len)
	      {
		for (unsigned int i=0; i<std::min(line_break,len); i++)
		  {
		    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		    OFSRealscientific(*out,17,val[cnt++]) << " ";
		  }
		libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		*out << '\n';
	      }
	  }

	return;	
      }

    default:
      libmesh_error();
    }
}


template <>
void Xdr::data_stream (long double *val, const unsigned int len, const unsigned int line_break)
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
	// doubles is disabled, but you can still write long double
	// ASCII files of course.
	// if (len > 0)
	//   xdr_vector(xdrs, 
	// 	     (char*) val,
	// 	     len,
	// 	     sizeof(double),
	// 	     (xdrproc_t) xdr_quadruple);

	std::cerr << "Writing binary XDR files with long double's is not\n"
		  << "currently supported on all platforms." << std::endl;

	libmesh_error();
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	for (unsigned int i=0; i<len; i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> val[i];
	  }

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	if (line_break == libMesh::invalid_uint)
	  for (unsigned int i=0; i<len; i++)
	    {
	      libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	      OFSRealscientific(*out,17,val[i]) << " ";
	    }
	else
	  {
	    unsigned int cnt=0;
	    while (cnt < len)
	      {
		for (unsigned int i=0; i<std::min(line_break,len); i++)
		  {
		    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		    OFSRealscientific(*out,17,val[cnt++]) << " ";
		  }
		libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		*out << '\n';
	      }
	  }

	return;	
      }

    default:
      libmesh_error();
    }
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
void Xdr::data_stream (std::complex<double> *val, const unsigned int len, const unsigned int line_break)
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
	      
	    xdr_vector(xdrs, 
		       (char*) &io_buffer[0],
		       2*len,
		       sizeof(double),
		       (xdrproc_t) xdr_double);
	    
	    // Fill val array if we are reading.
	    if (mode == DECODE)
	      for (unsigned int i=0, cnt=0; i<len; i++)
		{
		  val[i].real() = io_buffer[cnt++];
		  val[i].imag() = io_buffer[cnt++];
		} 
	  }
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure LIBMESH_HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	libmesh_error();

#endif
	return;
      }

    case READ:
      {
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());

	for (unsigned int i=0; i<len; i++)
	  {
	    libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	    *in >> val[i].real() >> val[i].imag();
	  }

	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());

	if (line_break == libMesh::invalid_uint)
	  for (unsigned int i=0; i<len; i++)
	    {
	      libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	      OFSRealscientific(*out,17,val[i].real()) << " ";
	      OFSRealscientific(*out,17,val[i].imag()) << " ";
	    }
	else
	  {
	    unsigned int cnt=0;
	    while (cnt < len)
	      {
		for (unsigned int i=0; i<std::min(line_break,len); i++)
		  {
		    libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		    OFSRealscientific(*out,17,val[cnt].real()) << " ";
		    OFSRealscientific(*out,17,val[cnt].imag()) << " ";
		    cnt++;
		  }
		libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
		*out << '\n';
	      }
	  }

	return;	
      }

    default:
      libmesh_error();
    }
}
#endif // # LIBMESH_USE_COMPLEX_NUMBERS

void Xdr::comment (std::string &comment)
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
	libmesh_assert (in.get() != NULL); libmesh_assert (in->good());
	in->getline(comm, comm_len);
	return;	
      }

    case WRITE:
      {
	libmesh_assert (out.get() != NULL); libmesh_assert (out->good());
	*out << "\t " << comment << '\n';
	return;	
      }

    default:
      libmesh_error();
    }
}


#undef xdr_REAL


//
template void Xdr::data_stream<int>          (int *val,          const unsigned int len, const unsigned int line_break);
template void Xdr::data_stream<unsigned int> (unsigned int *val, const unsigned int len, const unsigned int line_break);
