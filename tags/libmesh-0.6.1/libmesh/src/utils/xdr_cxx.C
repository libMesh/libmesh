// "$Id: xdr_cxx.C,v 1.26 2007-10-21 20:48:54 benkirk Exp $\n"

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <iostream>
#include <cstring>

// Local includes
#include "xdr_cxx.h"



//-------------------------------------------------------------
// Xdr class implementation
Xdr::Xdr (const std::string& name, const XdrMODE m) :
  mode(m),
#ifdef HAVE_XDR
  xdrs(NULL),
  fp(NULL),
#endif
  comm_len(xdr_MAX_STRING_LENGTH)
{
  this->open(name);
}



Xdr::~Xdr()
{
  close();
}



void Xdr::open (const std::string& name)
{
  if (name == "")
    return;

  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	fp = fopen(name.c_str(), (mode == ENCODE) ? "w" : "r");
	assert (fp);
	xdrs = new XDR;
	xdrstdio_create (xdrs, fp, (mode == ENCODE) ? XDR_ENCODE : XDR_DECODE);
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();
	
#endif
	return;

      }

    case READ:
      {
	in.open(name.c_str(), std::ios::in);
	assert (in.good());
	return;
      }

    case WRITE:
      {
	out.open(name.c_str(), std::ios::out);
	assert (out.good());
	return;
      }
      
    default:
      error();
    }  
}



void Xdr::close ()
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();
	
#endif
	return;
      }
      
    case READ:
      {
	if (in.is_open()) 
	  in.close();      
	return;
      }

    case WRITE:
      {
	if (out.is_open()) 
	  out.close();      
	return;
      }

    default:
      error();
    }
}



bool Xdr::is_open() const
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	if (fp)
	  if (xdrs)
	    return true;

	return false;

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

	return false;	

#endif

      }
      
    case READ:
      {
	return in.good();
      }

    case WRITE:
      {
	return out.good();
      }

    default:
      error();
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
#ifdef HAVE_XDR

	assert (is_open());

	xdr_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());
	
	out << a << "\t " << comment << '\n';

	return;
      }

    default:
      error();
    }
}



void Xdr::data (unsigned int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

	xdr_u_int(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}



void Xdr::data (short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}



void Xdr::data (unsigned short int& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_u_short(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}



void Xdr::data (float& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_float(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}



void Xdr::data (double& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	xdr_double(xdrs, &a);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	in >> a; in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a << "\t " << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}



#ifdef USE_COMPLEX_NUMBERS

void Xdr::data (std::complex<double>& a, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());
	double
	  _r=a.real(),
	  _i=a.imag();
	xdr_double(xdrs, &_r);
	xdr_double(xdrs, &_i);
	a = std::complex<double>(_r,_i);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());
	
	double _r, _i;
	in >> _r;
	in >> _i;
	a = std::complex<double>(_r,_i);
        in.getline(comm, comm_len);

	return;
      }

    case WRITE:
      {
	assert (out.good());

	out << a.real() << "\t " 
	    << a.imag() << "\t "
	    << comment << '\n';
	
	return;
      }

    default:
      error();
    }
}

#endif // USE_COMPLEX_NUMBERS



void Xdr::data (std::vector<int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(int),
		   (xdrproc_t) xdr_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    out << v[i] << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



void Xdr::data (std::vector<unsigned int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned int),
		   (xdrproc_t) xdr_u_int);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif	
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    out << v[i] << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



void Xdr::data (std::vector<short int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(short int),
		   (xdrproc_t) xdr_short);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    out << v[i] << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



void Xdr::data (std::vector<unsigned short int>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(unsigned short int),
		   (xdrproc_t) xdr_u_short);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    out << v[i] << " ";
	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



void Xdr::data (std::vector<float>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(float),
		   (xdrproc_t) xdr_float);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    OFSRealscientific(out,17,v[i]) << " ";
 	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



void Xdr::data (std::vector<double>& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

	unsigned int length = v.size();

	this->data(length, "# vector length");

	xdr_vector(xdrs, 
		   (char*) &v[0],
		   length,
		   sizeof(double),
		   (xdrproc_t) xdr_double);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (this->is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	this->data(length, "# vector length");
	
	// If you were expecting to read in a vector at this
	// point, it's not going to happen if length == 0!
	// assert (length != 0);
	
	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	    in >> v[i];
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length");

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (out.good());
	    OFSRealscientific(out,17,v[i]) << " ";
 	  }



	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}




#ifdef USE_COMPLEX_NUMBERS

void Xdr::data (std::vector< std::complex<double> >& v, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length = v.size();

	data(length, "# vector length x 2 (complex)");

	std::vector< std::complex<double> >::iterator iter = v.begin();
	
	for (; iter != v.end(); ++iter)
	    data(*iter, "");

// Alternative code
// 	/*
// 	 * save complex values as two std::vectors<double>.
// 	 * Using just one buffer increases time for copying,
// 	 * but reduces memory consumption
// 	 */
// 	std::vector<double> buf;
// 	buf.resize(length);

// 	// real part
// 	std::vector< std::complex<double> >::iterator c_iter   = v.begin();
// 	std::vector<double>::iterator                 buf_iter = buf.begin();
// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *buf_iter = c_iter->real();
// 	  ++buf_iter;
// 	}
// 	data(buf, "");

// 	// imaginary part
// 	c_iter   = v.begin();
// 	buf_iter = buf.begin();
// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *buf_iter = c_iter->real();
// 	  ++buf_iter;
// 	}
// 	data(buf, "");

// 	buf.clear();


// did not work...
// 	// with null pointer, let XDR dynamically allocate?
// 	xdr_vector(xdrs, 
// 		   (char*) &v[0],
// 		   length,
// 		   sizeof(std::complex<double>),
// 		   (xdrproc_t) 0);

#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif	
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

	unsigned int length=0;

	data(length, "# vector length x 2 (complex)");

	v.resize(length);

	std::vector< std::complex<double> >::iterator iter = v.begin();
	
	for (; iter != v.end(); ++iter)
	    data(*iter, "");

// alternative code
// 	/*
// 	 * load complex values as two std::vector<double>
// 	 * one after the other, store them in two buffers,
// 	 * since we have @e no chance to get the real and complex
// 	 * part one after the other into the std::complex
// 	 * (apart from messing with += or so, which i don't want to)
// 	 */
// 	std::vector<double> real_buf, imag_buf;
// 	real_buf.resize(length);
// 	imag_buf.resize(length);

// 	// get real & imaginary part
// 	data(real_buf, "");
// 	data(imag_buf, "");

// 	// copy into vector
// 	std::vector< std::complex<double> >::iterator c_iter   = v.begin();
// 	std::vector<double>::iterator                 real_buf_iter = real_buf.begin();
// 	std::vector<double>::iterator                 imag_buf_iter = imag_buf.begin();

// 	for (; c_iter != v.end(); ++c_iter)
// 	{
// 	  *c_iter = std::complex<double>(*real_buf_iter, *imag_buf_iter);
// 	  ++real_buf_iter;
// 	  ++imag_buf_iter;
// 	}

// 	// clear up
// 	real_buf.clear();
// 	imag_buf.clear();


// did not work...
// 	xdr_vector(xdrs, 
// 		   (char*) &v[0],
// 		   length,
// 		   sizeof(std::complex<double>),
// 		   (xdrproc_t) 0);
	
#else
	
	std::cerr << "ERROR: Functionality is not available." << std::endl
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	unsigned int length=0;

	data(length, "# vector length x 2 (complex)");

	v.resize(length);

	for (unsigned int i=0; i<v.size(); i++)
	  {
	    assert (in.good());
	
	    double _r, _i;
	    in >> _r;
	    in  >> _i;
	    v[i] = std::complex<double>(_r,_i);
	  }

	in.getline(comm, comm_len);

	return;	
      }

    case WRITE:
      {
	assert (out.good());

	unsigned int length=v.size();

	data(length, "# vector length x 2 (complex)");

 	for (unsigned int i=0; i<v.size(); i++)
 	  {
 	    assert (out.good());
	    OFSNumberscientific(out,17,v[i]) << " ";
 	  }

	out << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}

#endif // ifdef USE_COMPLEX_NUMBERS




void Xdr::data (std::string& s, const char* comment)
{
  switch (mode)
    {
    case ENCODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case DECODE:
      {
#ifdef HAVE_XDR

	assert (is_open());

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
		  << "Make sure HAVE_XDR is defined at build time" 
		  << std::endl
		  << "The XDR interface is not available in this installation"
		  << std::endl;

	error();

#endif
	return;
      }

    case READ:
      {
	assert (in.good());

	in.getline(comm, comm_len);

//#ifndef BROKEN_IOSTREAM
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
	assert (out.good());

	out << s << "\t " << comment << '\n';

	return;	
      }

    default:
      error();
    }
}



#undef xdr_REAL
