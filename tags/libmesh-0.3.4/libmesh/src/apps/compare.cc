#include "mesh_config.h"

// C++ includes
#include <iostream>
#include <vector>
#include <string>
#ifdef HAVE_GETOPT_H
// GCC 2.95.3 (and maybe others) do not include
// getopt.h in unistd.h...  Hower IBM xlC has no
// getopt.h!  This works around that.
#include <getopt.h>
#endif
#include <stdio.h>
#include <fstream>

// Local Includes
#include "libmesh.h"
#include "equation_systems.h"
#include "general_system.h"
#include "frequency_system.h"
#include "thin_system.h"
#include "mesh.h"
#include "perfmon.h"


/** 
 * convenient enums
 */
enum EqnSystemsFormat {EQNSYS_ASCII=0, 
		       EQNSYS_BINARY, 
		       EQNSYS_INVALID_IO};

enum SystemTypeHandled {GENERAL_SYSTEM=0,
#if defined(USE_COMPLEX_NUMBERS)
			FREQUENCY_SYSTEM, 
#endif
			THIN_SYSTEM,
			INVALID_SYSTEM};



/** 
 * how to use this, and command line processor
 */
void usage(char *progName)
{
  std::string baseName;
  static std::string helpList =
    "usage:\n"
    "	%s [options] ...\n"
    "\n"
    "options:\n"
    "    -m <string>                   Mesh file name\n"
    "    -d <dim>                      <dim>-dimensional mesh\n"
    "    -l <string>                   Left Equationn Systems file name\n"
    "    -r <string>                   Right Equation Systems file name\n"
    "    -t <float>                    threshold\n"
    "    -a                            ASCII format (default)\n"
    "    -b                            binary format\n"
    "    -G                            GeneralSystem system type (default)\n"
#if defined(USE_COMPLEX_NUMBERS)
    "    -F                            FrequencySystem system type\n"
#endif
    "    -T                            ThinSystem system type\n"
    "    -v                            Verbose\n"
    "    -q                            really quiet\n"
    "    -h                            Print help menu\n"
    "\n"
    "\n"
    " This program is used to compare equation systems to a user-specified\n"
    " threshold.  Equation systems are imported in the libMesh format\n"
    " provided through the read and write methods in class EquationSystems<T>.\n"
    " \n"
    "  ./compare -d 3 -m grid.xda -l leftfile.dat -r rightfile.dat -b -t 1.e-8\n"
    "\n"
    " will read in the mesh grid.xda, the equation systems leftfile.dat and\n"
    " rightfile.dat in binary format and compare systems, and especially the\n"
    " floats stored in vectors.  The comparison is said to be passed when the\n"
    " floating point values agree up to the given threshold.  When no threshold\n"
    " is set the default libMesh tolerance is used.  If neither -a or -b are set,\n"
    " ASCII format is assumed.  The same holds for -G,"
#if defined(USE_COMPLEX_NUMBERS)
" -F,"
#endif
    " and -T, where the\n"
    " default is -G.  For mesh import formats, please consult the meshtool help.\n"
    "\n"
    " Direct questions to:\n"
    " benkirk@cfdlab.ae.utexas.edu\n";

    
  if (progName == NULL)
    baseName = "UNKNOWN";
  else 
    baseName = progName;

  
  fprintf(stderr, helpList.c_str(), baseName.c_str());
  fflush(stderr);

  abort();
};



void process_cmd_line(int argc, char **argv,
		      std::vector<std::string>& names,
		      unsigned int& dim,
		      double& threshold,
		      EqnSystemsFormat& format,
		      SystemTypeHandled& type,
		      bool& verbose,
		      bool& quiet)
{
#if defined(USE_COMPLEX_NUMBERS)
  char optionStr[] =
    "d:m:l:r:t:abGFTvq?h";
#else
  char optionStr[] =
    "d:m:l:r:t:abGTvq?h";
#endif

  int opt;

  bool format_set    = false;
  bool type_set      = false;
  bool left_name_set = false;

  if (argc < 3)
    usage(argv[0]);


  while ((opt = getopt(argc, argv, optionStr)) != -1)
    {
      switch (opt)
	{
	  
	  /**
	   * Get mesh file name
	   */
	case 'm':
	  {
	    if (names.empty())
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Mesh file name must preceed left file name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }

	  /**
	   * Get the mesh dimension
	   */
	case 'd':
	  {
	    dim = atoi(optarg);
	    break;
	  }
	  
	  /**
	   * Get left file name
	   */
	case 'l':
	  {
	    if (!left_name_set)
	      {
	        names.push_back(optarg);
		left_name_set = true;
	      }
	    else
	      {
		std::cout << "ERROR: Mesh file name must preceed right file name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }
	  
	  /**
	   * Get output file name
	   */
	case 'r':
	  {
	    if ((!names.empty()) && (left_name_set))
	      names.push_back(optarg);
	    else
	      {
		std::cout << "ERROR: Mesh file name and left file name must preceed "
			  << "right file name!"
			  << std::endl;
		exit(1);
	      }
	    break;
	  }

	  /**
	   * Get the comparison threshold
	   */
	case 't':
	  {
	    threshold = atof(optarg);
	    break;
	  }

	  /**
	   * Use ascii format
	   */
	case 'a':
	  {
	    if (format_set)
	      {
		std::cout << "ERROR: Equation system file format already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
		format = EQNSYS_ASCII;
		format_set = true;
	      }
	    break;
	  }
	  
	  /**
	   * Use binary format
	   */
	case 'b':
	  {
	    if (format_set)
	      {
		std::cout << "ERROR: Equation system file format already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        format = EQNSYS_BINARY;
		format_set = true;
	      }
	    break;
	  }

	  
	  /**
	   * Load GeneralSystem
	   */
	case 'G':
	  {
	    if (type_set)
	      {
		std::cout << "ERROR: Equation system type already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        type = GENERAL_SYSTEM;
		type_set = true;
	      }
	    break;
	  }
	  

#if defined(USE_COMPLEX_NUMBERS)
	  /**
	   * Load FrequencySystem
	   */
	case 'F':
	  {
	    if (type_set)
	      {
		std::cout << "ERROR: Equation system type already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        type = FREQUENCY_SYSTEM;
		type_set = true;
	      }
	    break;
	  }
#endif	  


	  /**
	   * Load ThinSystem
	   */
	case 'T':
	  {
	    if (type_set)
	      {
		std::cout << "ERROR: Equation system type already set!"
			  << std::endl;
		exit(1);
	      }
	    else
	      {
	        type = THIN_SYSTEM;
		type_set = true;
	      }
	    break;
	  }
	  
	  /**
	   * Be verbose
	   */
	case 'v':
	  {
	    verbose = true;
	    break;
	  }
	  
	  /**
	   * Be totally quiet, no matter what -v says
	   */
	case 'q':
	  {
	    quiet = true;
	    break;
	  }
	  
	case 'h':
	case '?':
	  usage(argv[0]);
	  
	default:
	  return;
	}
    }

}





/**
 * template function for the different systems;
 * should get automatically instantiated due to
 * the calls below
 *
 * everything that is identical for the systems, and
 * should _not_ go into EquationSystems<T_sys>::compare(),
 * can go in this do_compare().
 */
template <typename T_sys>
bool do_compare (EquationSystems<T_sys>& les,
		 EquationSystems<T_sys>& res,
		 double threshold,
		 bool verbose)
{

  if (verbose)
    {
      std::cout << std::endl
		<< "*********   LEFT SYSTEM    *********" << std::endl;
      les.print_info  ();
      std::cout << "*********   RIGHT SYSTEM   *********" << std::endl;
      res.print_info ();
      std::cout << "********* COMPARISON PHASE *********" << std::endl;
    }
 
  /**
   * start comparing
   */
  bool result = les.compare(res, threshold, verbose);
  if (verbose)
    {
      std::cout << std::endl
		<< "*********     FINISHED     *********" << std::endl;
    }
  return result;
}










int main (int argc, char** argv)
{
  libMesh::init (argc, argv);
  
  // these should better be not contained in the following braces
  bool quiet                      = false;
  bool are_equal;

  {
    PerfMon perfmon(argv[0]);
    
    // default values
    std::vector<std::string> names;
    unsigned int dim                = static_cast<unsigned int>(-1);
    double threshold                = TOLERANCE;
    EqnSystemsFormat format         = EQNSYS_ASCII;
    SystemTypeHandled type          = GENERAL_SYSTEM;
    bool verbose                    = false;
 
    // get commands
    process_cmd_line(argc, argv, 
		     names,
		     dim,
		     threshold,
		     format,
		     type,
		     verbose,
		     quiet);


    if (dim == static_cast<unsigned int>(-1))
      {
	std::cout << "ERROR:  you must specify the dimension on "
		  << "the command line!\n\n"
		  << argv[0] << " -d 3 ... for example\n\n";
	error();
      }

    if (quiet)
      verbose = false;

    if (verbose)
      {
	  std::cout << "Settings:" << std::endl
		    << " dimensionality = " << dim << std::endl
		    << " mesh           = " << names[0] << std::endl
		    << " left system    = " << names[1] << std::endl
		    << " right system   = " << names[2] << std::endl
		    << " threshold      = " << threshold << std::endl
		    << " read format    = " << format << std::endl 
		    << std::endl;
      }	  


    /**
     * build the left and right mesh for left, inut them
     */
    Mesh left_mesh  (dim);
    Mesh right_mesh (dim);


    if (!names.empty())
      {
	left_mesh.read  (names[0]);
	right_mesh.read (names[0]);

	if (verbose)
	  left_mesh.print_info();
      }
    
    else
      {
	std::cout << "No input specified." << std::endl;
	return 1;
      }


    /**
     * build an appropriate EquationSystems object, read them
     */
    switch (type)
      {
	  
      case GENERAL_SYSTEM:
        {
	  EquationSystems<GeneralSystem> left_system  (left_mesh);
	  EquationSystems<GeneralSystem> right_system (right_mesh);

	  if (names.size() == 3)
	    {
	      switch (format)
	        {	  
		  case EQNSYS_ASCII:
		  {
		      left_system.read  (names[1], Xdr::READ);
		      right_system.read (names[2], Xdr::READ);
		      break;
		  }
		  
		  case EQNSYS_BINARY:
		  {
		      left_system.read  (names[1], Xdr::DECODE);
		      right_system.read (names[2], Xdr::DECODE);
		      break;
		  }

		  default:
		  {
		      std::cout << "Bad system type." << std::endl;
		      error();
		  }
		}
	  
	    }    
	  else
	    {
	      std::cout << "Bad input specified." << std::endl;
	      error();
	    }

	  are_equal = do_compare (left_system, right_system, threshold, verbose);

	  break;
	}


#if defined(USE_COMPLEX_NUMBERS)

      case FREQUENCY_SYSTEM:
        {
	  EquationSystems<FrequencySystem> left_system  (left_mesh);
	  EquationSystems<FrequencySystem> right_system (right_mesh);

	  if (names.size() == 3)
	    {
	      switch (format)
	        {	  
		  case EQNSYS_ASCII:
		  {
		      left_system.read  (names[1], Xdr::READ);
		      right_system.read (names[2], Xdr::READ);
		      break;
		  }
		  
		  case EQNSYS_BINARY:
		  {
		      left_system.read  (names[1], Xdr::DECODE);
		      right_system.read (names[2], Xdr::DECODE);
		      break;
		  }

		  default:
		  {
		      std::cout << "Bad system type." << std::endl;
		      error();
		  }
		}
	  
	    }    
	  else
	    {
	      std::cout << "Bad input specified." << std::endl;
	      error();
	    }

	  are_equal = do_compare (left_system, right_system, threshold, verbose);

	  break;
	}
	  
#endif


      case THIN_SYSTEM:
        {
	  EquationSystems<ThinSystem> left_system  (left_mesh);
	  EquationSystems<ThinSystem> right_system (right_mesh);

	  if (names.size() == 3)
	    {
	      switch (format)
	        {	  
		  case EQNSYS_ASCII:
		  {
		      left_system.read  (names[1], Xdr::READ);
		      right_system.read (names[2], Xdr::READ);
		      break;
		  }
		  
		  case EQNSYS_BINARY:
		  {
		      left_system.read  (names[1], Xdr::DECODE);
		      right_system.read (names[2], Xdr::DECODE);
		      break;
		  }

		  default:
		  {
		      std::cout << "Bad system type." << std::endl;
		      error();
		  }
		}
	  
	    }    
	  else
	    {
	      std::cout << "Bad input specified." << std::endl;
	      error();
	    }

	  are_equal = do_compare (left_system, right_system, threshold, verbose);

	  break;
	}


      default:
        {
	  std::cout << "Bad system type." << std::endl;
	  error();
	}
      }    
  }

  /**
   * let's see what do_compare found out
   */
  unsigned int our_result;

  if (are_equal)
    {
      if (!quiet)
	  std::cout << " Congrat's, up to the defined threshold, the two"  << std::endl
		    << " are identical." << std::endl;
      our_result=0;
    }
  else
    {
      if (!quiet)
        std::cout << " Oops, differences occured!"  << std::endl
		  << " Use -v to obtain more information where differences occured."
		  << std::endl;
      our_result=1;
    }

//  return libMesh::close();
  return our_result;
}
