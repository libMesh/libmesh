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
#define NC_EMAXDIMS     (-41)   // NC_MAX_DIMS exceeded
#define NC_ENAMEINUSE   (-42)   // String match to name in use
#define NC_ENOTATT      (-43)   // Attribute not found
#define NC_EMAXATTS     (-44)   // NC_MAX_ATTRS exceeded
#define NC_EBADTYPE     (-45)   // Not a netcdf data type
#define NC_EBADDIM      (-46)   // Invalid dimension id or name
#define NC_EUNLIMPOS    (-47)   // NC_UNLIMITED in the wrong index
#define NC_EMAXVARS     (-48)   // NC_MAX_VARS exceeded
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
~~~~

# NetCDF-4 Error Codes {#nc4-error-codes}

NetCDF-4 uses all error codes from NetCDF-3 (see section [NetCDF-3 Error
Codes](#NetCDF_002d3-Error-Codes)). The following additional error codes
were added for new errors unique to netCDF-4.

~~~~
#define NC_EHDFERR       (-101)
#define NC_ECANTREAD     (-102)
#define NC_ECANTWRITE    (-103)
#define NC_ECANTCREATE   (-104)
#define NC_EFILEMETA     (-105)
#define NC_EDIMMETA      (-106)
#define NC_EATTMETA      (-107)
#define NC_EVARMETA      (-108)
#define NC_ENOCOMPOUND   (-109)
#define NC_EATTEXISTS    (-110)
#define NC_ENOTNC4       (-111) // Attempting netcdf-4 operation on netcdf-3 file.
#define NC_ESTRICTNC3    (-112) // Attempting netcdf-4 operation on strict nc3 netcdf-4 file.
#define NC_EBADGRPID     (-113) // Bad group id. Bad!
#define NC_EBADTYPEID    (-114) // Bad type id.
#define NC_EBADFIELDID   (-115) // Bad field id.
#define NC_EUNKNAME      (-116)
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
~~~~
