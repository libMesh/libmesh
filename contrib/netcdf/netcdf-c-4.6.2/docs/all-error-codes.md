NetCDF Error Code Listing {#nc-error-codes}
==================

\tableofcontents

# NetCDF-3 Error Codes {#nc3-error-codes}

~~~~
#define NC_NOERR        0       // No Error
#define NC_EBADID       (-33)   // Not a netcdf id
#define NC_ENFILE       (-34)   // Too many netcdfs open
#define NC_EEXIST       (-35)   // netcdf file exists && NC_NOCLOBBER
#define NC_EINVAL       (-36)   // Invalid Argument
#define NC_EPERM        (-37)   // Write to read only
#define NC_ENOTINDEFINE (-38)   // Operation not allowed in data mode
#define NC_EINDEFINE    (-39)   // Operation not allowed in define mode
#define NC_EINVALCOORDS (-40)   // Index exceeds dimension bound
#define NC_EMAXDIMS     (-41)   // NC_MAX_DIMS exceeded  [not enforced after 4.5.0]
#define NC_ENAMEINUSE   (-42)   // String match to name in use
#define NC_ENOTATT      (-43)   // Attribute not found
#define NC_EMAXATTS     (-44)   // NC_MAX_ATTRS exceeded  [not enforced after 4.5.0]
#define NC_EBADTYPE     (-45)   // Not a netcdf data type
#define NC_EBADDIM      (-46)   // Invalid dimension id or name
#define NC_EUNLIMPOS    (-47)   // NC_UNLIMITED in the wrong index
#define NC_EMAXVARS     (-48)   // NC_MAX_VARS exceeded  [not enforced after 4.5.0]
#define NC_ENOTVAR      (-49)   // Variable not found
#define NC_EGLOBAL      (-50)   // Action prohibited on NC_GLOBAL varid
#define NC_ENOTNC       (-51)   // Not a netcdf file
#define NC_ESTS         (-52)   // In Fortran, string too short
#define NC_EMAXNAME     (-53)   // NC_MAX_NAME exceeded
#define NC_EUNLIMIT     (-54)   // NC_UNLIMITED size already in use
#define NC_ENORECVARS   (-55)   // nc_rec op when there are no record vars
#define NC_ECHAR        (-56)   // Attempt to convert between text & numbers
#define NC_EEDGE        (-57)   // Edge+start exceeds dimension bound
#define NC_ESTRIDE      (-58)   // Illegal stride
#define NC_EBADNAME     (-59)   // Attribute or variable name contains illegal characters

// N.B. following must match value in ncx.h

#define NC_ERANGE       (-60)   // Math result not representable
#define NC_ENOMEM       (-61)   // Memory allocation (malloc) failure
#define NC_EVARSIZE     (-62)   // One or more variable sizes violate format constraints
#define NC_EDIMSIZE     (-63)   // Invalid dimension size
#define NC_ETRUNC       (-64)   // File likely truncated or possibly corrupted
#define NC_EAXISTYPE    (-65)   // Unknown axis type
~~~~

# DAP Error Codes {#dap-error-codes}

If the DAP client is enabled, then the following additional error codes
may occur.

~~~~
#define NC_EDAP         (-66)   // Generic DAP error
#define NC_ECURL        (-67)   // Generic libcurl error
#define NC_EIO          (-68)   // Generic IO error
#define NC_ENODATA      (-69)   // Attempt to access variable with no data
#define NC_EDAPSVC      (-70)   // DAP Server side error
#define NC_EDAS         (-71)   // Malformed or inaccessible DAS
#define NC_EDDS         (-72)   // Malformed or inaccessible DDS
#define NC_EDATADDS     (-73)   // Malformed or inaccessible DATADDS
#define NC_EDAPURL      (-74)   // Malformed DAP URL
#define NC_EDAPCONSTRAINT (-75) // Malformed DAP Constraint
#define NC_ETRANSLATION (-76)   // Untranslatable construct
#define NC_EACCESS      (-77)   // Access Failure
#define NC_EAUTH        (-78)   // Authorization Failure
~~~~

# Misc. additional errors
~~~~
#define NC_ENOTFOUND     (-90)  // No such file
#define NC_ECANTREMOVE   (-91)  // Cannot remove file
#define NC_EINTERNAL     (-92)  // NetCDF Library Internal Error
#define NC_EPNETCDF      (-93)  // Error at PnetCDF layer
~~~~

# NetCDF-4 Error Codes {#nc4-error-codes}

NetCDF-4 uses all error codes from NetCDF-3 (see section [NetCDF-3 Error
Codes](#NetCDF_002d3-Error-Codes)). The following additional error codes
were added for new errors unique to netCDF-4.

~~~~
#define NC_EHDFERR       (-101)    // Error at HDF5 layer.
#define NC_ECANTREAD     (-102)    // Cannot read.
#define NC_ECANTWRITE    (-103)    // Cannot write.
#define NC_ECANTCREATE   (-104)    // Cannot create.
#define NC_EFILEMETA     (-105)    // Problem with file metadata.
#define NC_EDIMMETA      (-106)    // Problem with dimension metadata.
#define NC_EATTMETA      (-107)    // Problem with attribute metadata.
#define NC_EVARMETA      (-108)    // Problem with variable metadata.
#define NC_ENOCOMPOUND   (-109)    // Not a compound type.
#define NC_EATTEXISTS    (-110)    // Attribute already exists.
#define NC_ENOTNC4       (-111)    // Attempting netcdf-4 operation on netcdf-3 file.
#define NC_ESTRICTNC3    (-112)    // Attempting netcdf-4 operation on strict nc3 netcdf-4 file.
#define NC_ENOTNC3       (-113)    // Attempting netcdf-3 operation on netcdf-4 file.
#define NC_ENOPAR        (-114)    // Parallel operation on file opened for non-parallel access.
#define NC_EPARINIT      (-115)    // Error initializing for parallel access.
#define NC_EBADGRPID     (-116)    // Bad group ID.
#define NC_EBADTYPID     (-117)    // Bad type ID.
#define NC_ETYPDEFINED   (-118)    // Type has already been defined and may not be edited.
#define NC_EBADFIELD     (-119)    // Bad field ID.
#define NC_EBADCLASS     (-120)    // Bad class.
#define NC_EMAPTYPE      (-121)    // Mapped access for atomic types only.
#define NC_ELATEFILL     (-122)    // Attempt to define fill value when data already exists.
#define NC_ELATEDEF      (-123)    // Attempt to define var properties, like deflate, after enddef.
#define NC_EDIMSCALE     (-124)    // Problem with HDF5 dimscales.
#define NC_ENOGRP        (-125)    // No group found.
#define NC_ESTORAGE      (-126)    // Cannot specify both contiguous and chunking.
#define NC_EBADCHUNK     (-127)    // Bad chunksize.
#define NC_ENOTBUILT     (-128)    // Attempt to use feature that was not turned on when netCDF was built.
#define NC_EDISKLESS     (-129)    // Error in using diskless  access.
#define NC_ECANTEXTEND   (-130)    // Attempt to extend dataset during ind. I/O operation.
#define NC_EMPI          (-131)    // MPI operation failed.
#define NC_EFILTER       (-132)    // Filter operation failed.
#define NC_ERCFILE       (-133)    // RC file failure
#define NC_ENULLPAD      (-134)    // Header Bytes not Null-Byte padded
#define NC_EINMEMORY     (-135)    // In-memory file error
~~~~

# PnetCDF Error Codes {#pnetcdf-error-codes}

~~~~
#define NC_ESMALL                       (-201) // size of MPI_Offset too small for format
#define NC_ENOTINDEP                    (-202) // Operation not allowed in collective data mode
#define NC_EINDEP                       (-203) // Operation not allowed in independent data mode
#define NC_EFILE                        (-204) // Unknown error in file operation
#define NC_EREAD                        (-205) // Unknown error in reading file
#define NC_EWRITE                       (-206) // Unknown error in writing to file
#define NC_EOFILE                       (-207) // file open/creation failed
#define NC_EMULTITYPES                  (-208) // Multiple etypes used in MPI datatype
#define NC_EIOMISMATCH                  (-209) // Input/Output data amount mismatch
#define NC_ENEGATIVECNT                 (-210) // Negative count is specified
#define NC_EUNSPTETYPE                  (-211) // Unsupported etype in memory MPI datatype
#define NC_EINVAL_REQUEST               (-212) // invalid nonblocking request ID
#define NC_EAINT_TOO_SMALL              (-213) // MPI_Aint not large enough to hold requested value
#define NC_ENOTSUPPORT                  (-214) // feature is not yet supported
#define NC_ENULLBUF                     (-215) // trying to attach a NULL buffer
#define NC_EPREVATTACHBUF               (-216) // previous attached buffer is found
#define NC_ENULLABUF                    (-217) // no attached buffer is found
#define NC_EPENDINGBPUT                 (-218) // pending bput is found, cannot detach buffer
#define NC_EINSUFFBUF                   (-219) // attached buffer is too small
#define NC_ENOENT                       (-220) // File does not exist
#define NC_EINTOVERFLOW                 (-221) // Overflow when type cast to 4-byte integer
#define NC_ENOTENABLED                  (-222) // feature is not enabled
#define NC_EBAD_FILE                    (-223) // Invalid file name (e.g., path name too long)
#define NC_ENO_SPACE                    (-224) // Not enough space
#define NC_EQUOTA                       (-225) // Quota exceeded
#define NC_ENULLSTART                   (-226) // argument start is a NULL pointer
#define NC_ENULLCOUNT                   (-227) // argument count is a NULL pointer
#define NC_EINVAL_CMODE                 (-228) // Invalid file create mode
#define NC_ETYPESIZE                    (-229) // MPI derived data type size error (bigger than the variable size)
#define NC_ETYPE_MISMATCH               (-230) // element type of the MPI derived data type mismatches the variable type
#define NC_ETYPESIZE_MISMATCH           (-231) // file type size mismatches buffer type size
#define NC_ESTRICTCDF2                  (-232) // Attempting CDF-5 operation on CDF-2 file
#define NC_ENOTRECVAR                   (-233) // Attempting operation only for record variables
#define NC_ENOTFILL                     (-234) // Attempting to fill a variable when its fill mode is off
#define NC_EINVAL_OMODE                 (-235) // Invalid file open mode
#define NC_EPENDING                     (-236) // Pending nonblocking request is found at file close
#define NC_EMAX_REQ                     (-237) // Size of I/O request exceeds INT_MAX
#define NC_EBADLOG                      (-238) // Unrecognized log file format
#define NC_EFLUSHED                     (-239) // Nonblocking request has already been flushed to the PFS. It is too late to cancel
#define NC_EMULTIDEFINE                 (-250) // NC definitions inconsistent among processes
#define NC_EMULTIDEFINE_OMODE           (-251) // inconsistent file open modes among processes
#define NC_EMULTIDEFINE_DIM_NUM         (-252) // inconsistent number of dimensions
#define NC_EMULTIDEFINE_DIM_SIZE        (-253) // inconsistent size of dimension
#define NC_EMULTIDEFINE_DIM_NAME        (-254) // inconsistent dimension names
#define NC_EMULTIDEFINE_VAR_NUM         (-255) // inconsistent number of variables
#define NC_EMULTIDEFINE_VAR_NAME        (-256) // inconsistent variable name
#define NC_EMULTIDEFINE_VAR_NDIMS       (-257) // inconsistent variable number of dimensions
#define NC_EMULTIDEFINE_VAR_DIMIDS      (-258) // inconsistent variable dimension IDs
#define NC_EMULTIDEFINE_VAR_TYPE        (-259) // inconsistent variable data type
#define NC_EMULTIDEFINE_VAR_LEN         (-260) // inconsistent variable size
#define NC_EMULTIDEFINE_NUMRECS         (-261) // inconsistent number of records
#define NC_EMULTIDEFINE_VAR_BEGIN       (-262) // inconsistent variable file begin offset (internal use)
#define NC_EMULTIDEFINE_ATTR_NUM        (-263) // inconsistent number of attributes
#define NC_EMULTIDEFINE_ATTR_SIZE       (-264) // inconsistent memory space used by attribute (internal use)
#define NC_EMULTIDEFINE_ATTR_NAME       (-265) // inconsistent attribute name
#define NC_EMULTIDEFINE_ATTR_TYPE       (-266) // inconsistent attribute type
#define NC_EMULTIDEFINE_ATTR_LEN        (-267) // inconsistent attribute length
#define NC_EMULTIDEFINE_ATTR_VAL        (-268) // inconsistent attribute value
#define NC_EMULTIDEFINE_FNC_ARGS        (-269) // inconsistent function arguments used in collective API
#define NC_EMULTIDEFINE_FILL_MODE       (-270) // inconsistent dataset fill mode
#define NC_EMULTIDEFINE_VAR_FILL_MODE   (-271) // inconsistent variable fill mode
#define NC_EMULTIDEFINE_VAR_FILL_VALUE  (-272) // inconsistent variable fill value
#define NC_EMULTIDEFINE_CMODE           (-273) // inconsistent file create modes among processes
~~~~

