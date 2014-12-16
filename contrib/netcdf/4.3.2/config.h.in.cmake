#ifndef CONFIG_H
#define CONFIG_H

/* Eliminate a number of warnings which come up based on deprecated
   POSIX naming conventions. */
#ifdef _MSC_VER

/* Define O_BINARY so that the appropriate flags
are set when opening a binary file on Windows. */

/* Disable a few warnings under Visual Studio, for the
   time being. */
  #include <io.h>
  #pragma warning( disable: 4018 4996 4244 4305 )
  #define unlink _unlink
  #define open _open
  #define close _close
  #define read _read
  #define lseek _lseeki64
  
  #define fstat _fstat64

  #define off_t __int64
  #define _off_t __int64
  
  #ifndef _OFF_T_DEFINED
  #define _OFF_T_DEFINED
  #endif
  
#ifdef _WIN32
	#ifndef strcasecmp
  		#define strcasecmp _stricmp
		#define snprintf _snprintf
  	#endif
#endif


  #define strdup _strdup
  #define fdopen _fdopen
  #define write _write
  #define snprintf _snprintf
  #define strtoll _strtoi64
#endif


/* Cache Size, other variables for HDF5 */
#define DEFAULT_CHUNK_SIZE	${DEFAULT_CHUNK_SIZE}
#define DEFAULT_CHUNKS_IN_CACHE	${DEFAULT_CHUNKS_IN_CACHE}
#define CHUNK_CACHE_SIZE 	${CHUNK_CACHE_SIZE}
#define CHUNK_CACHE_NELEMS 	${CHUNK_CACHE_NELEMS}
#define CHUNK_CACHE_PREEMPTION 	${CHUNK_CACHE_PREEMPTION}
#define MAX_DEFAULT_CACHE_SIZE	${MAX_DEFAULT_CACHE_SIZE}
#define NCIO_MINBLOCKSIZE	${NCIO_MINBLOCKSIZE}

#ifndef _FILE_OFFSET_BITS
#cmakedefine _FILE_OFFSET_BITS ${_FILE_OFFSET_BITS}
#cmakedefine _LARGEFILE64_SOURCE
#cmakedefine _LARGEFILE_SOURCE
#endif

#define PACKAGE_VERSION "${VERSION}"
#cmakedefine VERSION "${VERSION}"
#cmakedefine NC_VERSION "${VERSION}"
/* For HDF5 use. */
#cmakedefine H5_USE_16_API 1

/* Enable Logging, only valid for netcdf 4. */
#cmakedefine LOGGING 1

/* Various other options. */
#cmakedefine BUILD_V2 1
#cmakedefine NO_NETCDF_2 1
#cmakedefine USE_FSYNC 1
#cmakedefine JNA 1
#cmakedefine ENABLE_DOXYGEN 1
#cmakedefine ENABLE_INTERNAL_DOCS 1
#cmakedefine VALGRIND_TESTS 1
#cmakedefine ENABLE_CDMREMOTE 1
#cmakedefine USE_DAP 1
#cmakedefine ENABLE_DAP 1
#cmakedefine ENABLE_DAP_GROUPS 1
#cmakedefine ENABLE_DAP_REMOTE_TESTS 1
#cmakedefine EXTRA_TESTS
#cmakedefine USE_NETCDF4 1
#cmakedefine USE_LIBDL 1
#cmakedefine USE_HDF4 1
#cmakedefine USE_HDF5 1
#cmakedefine USE_FFIO 1
#cmakedefine USE_PARALLEL_POSIX 1
#cmakedefine USE_PARALLEL_MPIO 1
#cmakedefine USE_PARALLEL 1
#cmakedefine USE_PNETCDF 1
#cmakedefine USE_MMAP 1
#cmakedefine TEST_PARALLEL ${TEST_PARALLEL}
#cmakedefine BUILD_RPC 1
#cmakedefine USE_DISKLESS 1
#cmakedefine USE_SZIP 1
#cmakedefine USE_ZLIB 1
#cmakedefine USE_X_GETOPT 1
#cmakedefine LARGE_FILE_TESTS 1
#cmakedefine HAVE_DECL_ISFINITE 1
#cmakedefine HAVE_DECL_ISNAN 1
#cmakedefine HAVE_CURLOPT_USERNAME 1
#cmakedefine HAVE_CURLOPT_PASSWORD 1
#cmakedefine HAVE_CURLOPT_KEYPASSWD 1
#cmakedefine HAVE_CURLINFO_RESPONSE_CODE 1
#cmakedefine HAVE_DECL_SIGNBIT 1
#cmakedefine HAVE_DOPRNT
#cmakedefine HAVE_ALLOCA
#cmakedefine HAVE_SSIZE_T 1
#cmakedefine HAVE_LIBPNETCDF 1
#cmakedefine HAVE_LIBDL 1

/* Define to 1 if you have the <alloca.h> header file. */
#cmakedefine HAVE_ALLOCA_H @HAVE_ALLOCA_H@

/* Define to 1 if you have the <ctype.h> header file. */
#cmakedefine HAVE_CTYPE_H @HAVE_CTYPE_H@

/* Define to 1 if you have the <dirent> header file. */
#cmakedefine HAVE_DIRENT_H @HAVE_DIRENT_H@

/* Define to 1 if you have the <unistd.h> header file. */
#cmakedefine HAVE_UNISTD_H @HAVE_UNISTD_H@
#cmakedefine YY_NO_UNISTD_H @YY_NO_UNISTD_H@

/* Define to 1 if you have the <dlfcn.h> header file. */
#cmakedefine HAVE_DLFCN_H @HAVE_DLFCN_H@

/* Define to 1 if you have the <errno.h> header file. */
#cmakedefine HAVE_ERRNO_H @HAVE_ERRNO_H@

/* Define to 1 if you have the <fcntl.h> header file. */
#cmakedefine HAVE_FCNTL_H @HAVE_FCNTL_H@

/* Define to 1 if you have the <getopt.h> header file. */
#cmakedefine HAVE_GETOPT_H @HAVE_GETOPT_H@

/* Define to 1 if you have the <stdarg.h> header file. */
#cmakedefine HAVE_STDARG_H @HAVE_STDARG_H@

/* Define to 1 if you have the <hdf5.h> header file. */
#cmakedefine HAVE_HDF5_H @HAVE_HDF5_H@

/* Define to 1 if you have the <hdf5_hl.h> header file. */
#cmakedefine HAVE_HDF5_HL_H @HAVE_HDF5_HL_H@

/* Define to 1 if you have the <stdbool.h> header file. */
#cmakedefine HAVE_STDBOOL_H @HAVE_STDBOOL_H@

/* Define to 1 if you have the <locale.h> header file. */
#cmakedefine HAVE_LOCAL_H @HAVE_LOCAL_H@

/* Define to 1 if you have the <stdint.h> header file. */
#cmakedefine HAVE_STDINT_H @HAVE_STDINT_H@

/* Define to 1 if you have the <stdio.h> header file. */
#cmakedefine HAVE_STDIO_H @HAVE_STDIO_H@

/* Define to 1 if you have the <stdlib.h> header file. */
#cmakedefine HAVE_STDLIB_H @HAVE_STDLIB_H@

/* Define to 1 if you have the <strings.h> header file. */
#cmakedefine HAVE_STRINGS_H @HAVE_STRINGS_H@

/* Define to 1 if you have the <signal.h> header file. */
#cmakedefine HAVE_SIGNAL_H @HAVE_SIGNAL_H@

/* Define to 1 if you have the <sys/dir.h> header file. */
#cmakedefine HAVE_SYS_DIR_H @HAVE_SYS_DIR_H@

/* Define to 1 if you have the <sys/ndir.h> header file. */
#cmakedefine HAVE_SYS_NDIR_H @HAVE_SYS_NDIR_H@

/* Define to 1 if you have the <sys/param.h> header file. */
#cmakedefine HAVE_SYS_PARAM_H @HAVE_SYS_PARAM_H@

/* Define to 1 if you have the <sys/stat.h> header file. */
#cmakedefine HAVE_SYS_STAT_H @HAVE_SYS_STAT_H@

/* Define to 1 if you have the <sys/time.h> header file. */
#cmakedefine HAVE_SYS_TIME_H @HAVE_SYS_TIME_H@

/* Define to 1 if you have the <sys/resource.h> header file. */
#cmakedefine HAVE_SYS_RESOURCE_H @HAVE_SYS_RESOURCE_H@

/* Define to 1 if you have the <sys/types.h> header file. */
#cmakedefine HAVE_SYS_TYPES_H @HAVE_SYS_TYPES_H@

/* Define to 1 if you have the <sys/wait.h> header file. */
#cmakedefine HAVE_SYS_WAIT_H @HAVE_SYS_WAIT_H@

/* Define to 1 if you have the <inttypes.h> header file. */
#cmakedefine HAVE_INTTYPES_H @HAVE_INTTYPES_H@

/* Define to 1 if you have the <fcntl.h> header file. */
#cmakedefine HAVE_FCNTL_H @HAVE_FCNTL_H@

/* Define to 1 if you have the <malloc.h> header file. */
#cmakedefine HAVE_MALLOC_H @HAVE_MALLOC_H@

/* Define to 1 if you have the BaseTsd.h header file. */
#cmakedefine HAVE_BASETSD_H @HAVE_BASETSD_H@

/* Define if we have filelengthi64. */
#cmakedefine HAVE_FILE_LENGTH_I64 @HAVE_FILE_LENGTH_I64@

/* The size of `char` as computed by sizeof. */
#cmakedefine SIZEOF_CHAR @SIZEOF_CHAR@
/* The size of `double` as computed by sizeof. */
#cmakedefine SIZEOF_DOUBLE @SIZEOF_DOUBLE@ 
/* The size of `float` as computed by sizeof. */
#cmakedefine SIZEOF_FLOAT @SIZEOF_FLOAT@
/* The size of `int` as computed by sizeof. */
#cmakedefine SIZEOF_INT @SIZEOF_INT@
/* The size of `long` as computed by sizeof. */
#cmakedefine SIZEOF_LONG @SIZEOF_LONG@
/* The size of `long long` as computed by sizeof. */
#cmakedefine SIZEOF_LONG_LONG @SIZEOF_LONG_LONG@
/* The size of `off_t` as computed by sizeof. */
#cmakedefine SIZEOF_OFF_T @SIZEOF_OFF_T@
/* The size of `short` as computed by sizeof. */
#cmakedefine SIZEOF_OFF64_T @SIZEOF_OFF64_T@
#cmakedefine SIZEOF_SHORT @SIZEOF_SHORT@
/* The size of `size_t` as computed by sizeof. */
#cmakedefine SIZEOF_SIZE_T @SIZEOF_SIZE_T@
/* The size of `ssize_t` as computed by sizeof. */
#cmakedefine SIZEOF_SSIZE_T @SIZEOF_SSIZE_T@ 
/* The size of `uchar` as computed by sizeof. */
#cmakedefine SIZEOF_UCHAR @SIZEOF_UCHAR@
/* The size of `__int64` found on Windows systems. */
#cmakedefine SIZEOF___INT64 @SIZEOF___INT64@



#cmakedefine TEMP_LARGE "@TEMP_LARGE@"

/* Set if we have strdup */
#cmakedefine HAVE_STRDUP
#cmakedefine HAVE_STRLCAT
#cmakedefine HAVE_STRERROR
#cmakedefine HAVE_SNPRINTF
#cmakedefine HAVE_STRCHR
#cmakedefine HAVE_STRRCHR
#cmakedefine HAVE_STRCAT
#cmakedefine HAVE_STRCPY
#cmakedefine HAVE_STRDUP
#cmakedefine HAVE_STRCASECMP
#cmakedefine HAVE_STRTOD
#cmakedefine HAVE_STRTOLL
#cmakedefine HAVE_STROULL
#cmakedefine HAVE_STRSTR
#cmakedefine HAVE_MKSTEMP
#cmakedefine HAVE_RAND
#cmakedefine HAVE_GETTIMEOFDAY
#cmakedefine HAVE_MPI_COMM_F2C
#cmakedefine HAVE_MEMMOVE
#cmakedefine HAVE_MMAP
#cmakedefine HAVE_GETPAGESIZE
#cmakedefine HAVE_SYSCONF
#cmakedefine HAVE_MREMAP
#cmakedefine HAVE_DECL_ISINF

#cmakedefine HAVE_GETRLIMIT
#cmakedefine HAVE_FSYNC

#cmakedefine HAVE_H5PGET_FAPL_MPIPOSIX 1
#cmakedefine HAVE_H5PSET_DEFLATE
#cmakedefine HAVE_H5Z_SZIP


/* Specifies if various libraries are present. */
#cmakedefine HAVE_LIBM 1

/* Define to 1 if the system has the type `uchar'.*/ 
#cmakedefine HAVE_UCHAR 

/* Misc defines copied from autotools config.h.in */
#cmakedefine CRAY_STACKSEG_END
#cmakedefine DLL_EXPORT
#cmakedefine DLL_NETCDF

#include <ncconfigure.h>
#endif
