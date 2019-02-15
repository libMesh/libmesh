Appendix B. File Format Specifications {#file_format_specifications}
=======================================

[TOC]

\section classic_format_spec The NetCDF Classic Format Specification

To present the format more formally, we use a BNF grammar notation. In
this notation:
- Non-terminals (entities defined by grammar rules) are in lower case.
- Terminals (atomic entities in terms of which the format
specification is written) are in upper case, and are specified
literally as US-ASCII characters within single-quote characters or are
described with text between angle brackets (‘\<’ and ‘\>’).
- Optional entities are enclosed between braces (‘[’ and ‘]’).
- A sequence of zero or more occurrences of an entity is denoted by
  ‘[entity ...]’.
- A vertical line character (‘|’) separates alternatives. Alternation
  has lower precedence than concatenation.
- Comments follow ‘//’ characters.
- A single byte that is not a printable character is denoted using a
hexadecimal number with the notation ‘\\xDD’, where each D is a
hexadecimal digit.
- A literal single-quote character is denoted by ‘\'’, and a literal
back-slash character is denoted by ‘\\’.

Following the grammar, a few additional notes are included to specify
format characteristics that are impractical to capture in a BNF
grammar, and to note some special cases for implementers. Comments in
the grammar point to the notes and special cases, and help to clarify
the intent of elements of the format.

<h1>The Format in Detail</h1>

<pre>
     netcdf_file  = header  data
     header       = magic  numrecs  dim_list  gatt_list  var_list
     magic        = 'C'  'D'  'F'  VERSION
     VERSION      = \\x01 |                      // classic format
                    \\x02                        // 64-bit offset format
     numrecs      = NON_NEG | STREAMING         // length of record dimension
     dim_list     = ABSENT | NC_DIMENSION  nelems  [dim ...]
     gatt_list    = att_list                    // global attributes
     att_list     = ABSENT | NC_ATTRIBUTE  nelems  [attr ...]
     var_list     = ABSENT | NC_VARIABLE   nelems  [var ...]
     ABSENT       = ZERO  ZERO                  // Means list is not present
     ZERO         = \\x00 \\x00 \\x00 \\x00         // 32-bit zero
     NC_DIMENSION = \\x00 \\x00 \\x00 \\x0A         // tag for list of dimensions
     NC_VARIABLE  = \\x00 \\x00 \\x00 \\x0B         // tag for list of variables
     NC_ATTRIBUTE = \\x00 \\x00 \\x00 \\x0C         // tag for list of attributes
     nelems       = NON_NEG       // number of elements in following sequence
     dim          = name  dim_length
     name         = nelems  namestring
                         // Names a dimension, variable, or attribute.
                         // Names should match the regular expression
                         // ([a-zA-Z0-9_]|{MUTF8})([^\\x00-\\x1F/\\x7F-\\xFF]|{MUTF8})*
                         // For other constraints, see "Note on names", below.
     namestring   = ID1 [IDN ...] padding
     ID1          = alphanumeric | '_'
     IDN          = alphanumeric | special1 | special2
     alphanumeric = lowercase | uppercase | numeric | MUTF8
     lowercase    = 'a'|'b'|'c'|'d'|'e'|'f'|'g'|'h'|'i'|'j'|'k'|'l'|'m'|
                    'n'|'o'|'p'|'q'|'r'|'s'|'t'|'u'|'v'|'w'|'x'|'y'|'z'
     uppercase    = 'A'|'B'|'C'|'D'|'E'|'F'|'G'|'H'|'I'|'J'|'K'|'L'|'M'|
                    'N'|'O'|'P'|'Q'|'R'|'S'|'T'|'U'|'V'|'W'|'X'|'Y'|'Z'
     numeric      = '0'|'1'|'2'|'3'|'4'|'5'|'6'|'7'|'8'|'9'
                                  // special1 chars have traditionally been
                                  // permitted in netCDF names.
     special1     = '_'|'.'|'@'|'+'|'-'
                                  // special2 chars are recently permitted in
                                  // names (and require escaping in CDL).
                                  // Note: '/' is not permitted.
     special2     = ' ' | '!' | '"' | '#'  | '$' | '%' | '&' | '\'' |
                    '(' | ')' | '*' | ','  | ':' | ';' | '<' | '='  |
                    '>' | '?' | '[' | '\\' | ']' | '^' | '`' | '{'  |
                    '|' | '}' | '~'
     MUTF8        = <multibyte UTF-8 encoded, NFC-normalized Unicode character>
     dim_length   = NON_NEG       // If zero, this is the record dimension.
                                  // There can be at most one record dimension.
     attr         = name  nc_type  nelems  [values ...]
     nc_type      = NC_BYTE   |
                    NC_CHAR   |
                    NC_SHORT  |
                    NC_INT    |
                    NC_FLOAT  |
                    NC_DOUBLE
     var          = name  nelems  [dimid ...]  vatt_list  nc_type  vsize  begin
                                  // nelems is the dimensionality (rank) of the
                                  // variable: 0 for scalar, 1 for vector, 2
                                  // for matrix, ...
     dimid        = NON_NEG       // Dimension ID (index into dim_list) for
                                  // variable shape.  We say this is a "record
                                  // variable" if and only if the first
                                  // dimension is the record dimension.
     vatt_list    = att_list      // Variable-specific attributes
     vsize        = NON_NEG       // Variable size.  If not a record variable,
                                  // the amount of space in bytes allocated to
                                  // the variable's data.  If a record variable,
                                  // the amount of space per record.  See "Note
                                  // on vsize", below.
     begin        = OFFSET        // Variable start location.  The offset in
                                  // bytes (seek index) in the file of the
                                  // beginning of data for this variable.
     data         = non_recs  recs
     non_recs     = [vardata ...] // The data for all non-record variables,
                                  // stored contiguously for each variable, in
                                  // the same order the variables occur in the
                                  // header.
     vardata      = [values ...]  // All data for a non-record variable, as a
                                  // block of values of the same type as the
                                  // variable, in row-major order (last
                                  // dimension varying fastest).
     recs         = [record ...]  // The data for all record variables are
                                  // stored interleaved at the end of the
                                  // file.
     record       = [varslab ...] // Each record consists of the n-th slab
                                  // from each record variable, for example
                                  // x[n,...], y[n,...], z[n,...] where the
                                  // first index is the record number, which
                                  // is the unlimited dimension index.
     varslab      = [values ...]  // One record of data for a variable, a
                                  // block of values all of the same type as
                                  // the variable in row-major order (last
                                  // index varying fastest).
     values       = bytes | chars | shorts | ints | floats | doubles
     string       = nelems  [chars]
     bytes        = [BYTE ...]  padding
     chars        = [CHAR ...]  padding
     shorts       = [SHORT ...]  padding
     ints         = [INT ...]
     floats       = [FLOAT ...]
     doubles      = [DOUBLE ...]
     padding      = <0, 1, 2, or 3 bytes to next 4-byte boundary>
                                  // Header padding uses null (\\x00) bytes.  In
                                  // data, padding uses variable's fill value.
                                  // See "Note on padding", below, for a special
                                  // case.
     NON_NEG      = <non-negative INT>
     STREAMING    = \\xFF \\xFF \\xFF \\xFF   // Indicates indeterminate record
                                          // count, allows streaming data
     OFFSET       = <non-negative INT> |  // For classic format or
                    <non-negative INT64>  // for 64-bit offset format
     BYTE         = <8-bit byte>          // See "Note on byte data", below.
     CHAR         = <8-bit byte>          // See "Note on char data", below.
     SHORT        = <16-bit signed integer, Bigendian, two's complement>
     INT          = <32-bit signed integer, Bigendian, two's complement>
     INT64        = <64-bit signed integer, Bigendian, two's complement>
     FLOAT        = <32-bit IEEE single-precision float, Bigendian>
     DOUBLE       = <64-bit IEEE double-precision float, Bigendian>
                                  // following type tags are 32-bit integers
     NC_BYTE      = \\x00 \\x00 \\x00 \\x01       // 8-bit signed integers
     NC_CHAR      = \\x00 \\x00 \\x00 \\x02       // text characters
     NC_SHORT     = \\x00 \\x00 \\x00 \\x03       // 16-bit signed integers
     NC_INT       = \\x00 \\x00 \\x00 \\x04       // 32-bit signed integers
     NC_FLOAT     = \\x00 \\x00 \\x00 \\x05       // IEEE single precision floats
     NC_DOUBLE    = \\x00 \\x00 \\x00 \\x06       // IEEE double precision floats
                                  // Default fill values for each type, may be
                                  // overridden by variable attribute named
                                  // '_FillValue'. See "Note on fill values",
                                  // below.
     FILL_CHAR    = \\x00                      // null byte
     FILL_BYTE    = \\x81                      // (signed char) -127
     FILL_SHORT   = \\x80 \\x01                 // (short) -32767
     FILL_INT     = \\x80 \\x00 \\x00 \\x01       // (int) -2147483647
     FILL_FLOAT   = \\x7C \\xF0 \\x00 \\x00       // (float) 9.9692099683868690e+36
     FILL_DOUBLE  = \\x47 \\x9E \\x00 \\x00 \\x00 \\x00 \\x00 \\x00 //(double)9.9692099683868690e+36
</pre>

Note on vsize: This number is the product of the dimension lengths
(omitting the record dimension) and the number of bytes per value
(determined from the type), increased to the next multiple of 4, for
each variable. If a record variable, this is the amount of space per
record (except that, for backward compatibility, it always includes
padding to the next multiple of 4 bytes, even in the exceptional case
noted below under “Note on padding”). The netCDF “record size” is
calculated as the sum of the vsize's of all the record variables.

The vsize field is actually redundant, because its value may be
computed from other information in the header. The 32-bit vsize field
is not large enough to contain the size of variables that require more
than 2^32 - 4 bytes, so 2^32 - 1 is used in the vsize field for such
variables.

Note on names: Earlier versions of the netCDF C-library reference
implementation enforced a more restricted set of characters in
creating new names, but permitted reading names containing arbitrary
bytes. This specification extends the permitted characters in names to
include multi-byte UTF-8 encoded Unicode and additional printing
characters from the US-ASCII alphabet. The first character of a name
must be alphanumeric, a multi-byte UTF-8 character, or '_' (reserved
for special names with meaning to implementations, such as the
“_FillValue” attribute). Subsequent characters may also include
printing special characters, except for '/' which is not allowed in
names. Names that have trailing space characters are also not
permitted.

Implementations of the netCDF classic and 64-bit offset format must
ensure that names are normalized according to Unicode NFC
normalization rules during encoding as UTF-8 for storing in the file
header. This is necessary to ensure that gratuitous differences in the
representation of Unicode names do not cause anomalies in comparing
files and querying data objects by name.

Note on streaming data: The largest possible record count, 2^32 - 1,
is reserved to indicate an indeterminate number of records. This means
that the number of records in the file must be determined by other
means, such as reading them or computing the current number of records
from the file length and other information in the header. It also
means that the numrecs field in the header will not be updated as
records are added to the file. [This feature is not yet implemented].

Note on padding: In the special case when there is only one record
variable and it is of type character, byte, or short, no padding is
used between record slabs, so records after the first record do not
necessarily start on four-byte boundaries. However, as noted above
under “Note on vsize”, the vsize field is computed to include padding
to the next multiple of 4 bytes. In this case, readers should ignore
vsize and assume no padding. Writers should store vsize as if padding
were included.

Note on byte data: It is possible to interpret byte data as either
signed (-128 to 127) or unsigned (0 to 255). When reading byte data
through an interface that converts it into another numeric type, the
default interpretation is signed. There are various attribute
conventions for specifying whether bytes represent signed or unsigned
data, but no standard convention has been established. The variable
attribute “_Unsigned” is reserved for this purpose in future
implementations.

Note on char data: Although the characters used in netCDF names must
be encoded as UTF-8, character data may use other encodings. The
variable attribute “_Encoding” is reserved for this purpose in future
implementations.

Note on fill values: Because data variables may be created before
their values are written, and because values need not be written
sequentially in a netCDF file, default “fill values” are defined for
each type, for initializing data values before they are explicitly
written. This makes it possible to detect reading values that were
never written. The variable attribute “_FillValue”, if present,
overrides the default fill value for a variable. If _FillValue is
defined then it should be scalar and of the same type as the variable.

Fill values are not required, however, because netCDF libraries have
traditionally supported a “no fill” mode when writing, omitting the
initialization of variable values with fill values. This makes the
creation of large files faster, but also eliminates the possibility of
detecting the inadvertent reading of values that haven't been written.

\section computing_offsets Notes on Computing File Offsets

The offset (position within the file) of a specified data value in a
classic format or 64-bit offset data file is completely determined by
the variable start location (the offset in the begin field), the
external type of the variable (the nc_type field), and the dimension
indices (one for each of the variable's dimensions) of the value
desired.

The external size in bytes of one data value for each possible netCDF
type, denoted extsize below, is:
- NC_BYTE 	1
- NC_CHAR 	1
- NC_SHORT 	2
- NC_INT 	4
- NC_FLOAT 	4
- NC_DOUBLE 	8

The record size, denoted by recsize below, is the sum of the vsize
fields of record variables (variables that use the unlimited
dimension), using the actual value determined by dimension sizes and
variable type in case the vsize field is too small for the variable
size.

To compute the offset of a value relative to the beginning of a
variable, it is helpful to precompute a “product vector” from the
dimension lengths. Form the products of the dimension lengths for the
variable from right to left, skipping the leftmost (record) dimension
for record variables, and storing the results as the product vector
for each variable.

For example:

\code
Non-record variable:

dimension lengths: [ 5 3 2 7] product vector: [210 42 14 7]

Record variable:

dimension lengths: [0 2 9 4] product vector: [0 72 36 4]
\endcode

At this point, the leftmost product, when rounded up to the next
multiple of 4, is the variable size, vsize, in the grammar above. For
example, in the non-record variable above, the value of the vsize
field is 212 (210 rounded up to a multiple of 4). For the record
variable, the value of vsize is just 72, since this is already a
multiple of 4.

Let coord be the array of coordinates (dimension indices, zero-based)
of the desired data value. Then the offset of the value from the
beginning of the file is just the file offset of the first data value
of the desired variable (its begin field) added to the inner product
of the coord and product vectors times the size, in bytes, of each
datum for the variable. Finally, if the variable is a record variable,
the product of the record number, 'coord[0]', and the record size,
recsize, is added to yield the final offset value.

A special case: Where there is exactly one record variable, we drop
the requirement that each record be four-byte aligned, so in this case
there is no record padding.

\subsection offset_examples Examples

By using the grammar above, we can derive the smallest valid netCDF
file, having no dimensions, no variables, no attributes, and hence, no
data. A CDL representation of the empty netCDF file is

\code
netcdf empty { }
\endcode

This empty netCDF file has 32 bytes. It begins with the four-byte
“magic number” that identifies it as a netCDF version 1 file: ‘C’,
‘D’, ‘F’, ‘\\x01’. Following are seven 32-bit integer zeros
representing the number of records, an empty list of dimensions, an
empty list of global attributes, and an empty list of variables.

Below is an (edited) dump of the file produced using the Unix command

\code
od -xcs empty.nc
\endcode

Each 16-byte portion of the file is displayed with 4 lines. The first
line displays the bytes in hexadecimal. The second line displays the
bytes as characters. The third line displays each group of two bytes
interpreted as a signed 16-bit integer. The fourth line (added by
human) presents the interpretation of the bytes in terms of netCDF
components and values.

\code
        4344    4601    0000    0000    0000    0000    0000    0000
       C   D   F 001  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
       17220   17921   00000   00000   00000   00000   00000   00000
     [magic number ] [  0 records  ] [  0 dimensions   (ABSENT)    ]

        0000    0000    0000    0000    0000    0000    0000    0000
      \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
       00000   00000   00000   00000   00000   00000   00000   00000
     [  0 global atts  (ABSENT)    ] [  0 variables    (ABSENT)    ]
\endcode

As a less trivial example, consider the CDL

\code
     netcdf tiny {
     dimensions:
             dim = 5;
     variables:
             short vx(dim);
     data:
             vx = 3, 1, 4, 1, 5 ;
     }
\endcode

which corresponds to a 92-byte netCDF file. The following is an edited dump of this file:

\code
        4344    4601    0000    0000    0000    000a    0000    0001
       C   D   F 001  \0  \0  \0  \0  \0  \0  \0  \n  \0  \0  \0 001
       17220   17921   00000   00000   00000   00010   00000   00001
     [magic number ] [  0 records  ] [NC_DIMENSION ] [ 1 dimension ]

        0000    0003    6469    6d00    0000    0005    0000    0000
      \0  \0  \0 003   d   i   m  \0  \0  \0  \0 005  \0  \0  \0  \0
       00000   00003   25705   27904   00000   00005   00000   00000
     [  3 char name = "dim"        ] [ size = 5    ] [ 0 global atts

        0000    0000    0000    000b    0000    0001    0000    0002
      \0  \0  \0  \0  \0  \0  \0 013  \0  \0  \0 001  \0  \0  \0 002
       00000   00000   00000   00011   00000   00001   00000   00002
      (ABSENT)     ] [NC_VARIABLE  ] [ 1 variable  ] [ 2 char name =

        7678    0000    0000    0001    0000    0000    0000    0000
       v   x  \0  \0  \0  \0  \0 001  \0  \0  \0  \0  \0  \0  \0  \0
       30328   00000   00000   00001   00000   00000   00000   00000
      "vx"         ] [1 dimension  ] [ with ID 0   ] [ 0 attributes

        0000    0000    0000    0003    0000    000c    0000    0050
      \0  \0  \0  \0  \0  \0  \0 003  \0  \0  \0  \f  \0  \0  \0   P
       00000   00000   00000   00003   00000   00012   00000   00080
      (ABSENT)     ] [type NC_SHORT] [size 12 bytes] [offset:    80]

        0003    0001    0004    0001    0005    8001
      \0 003  \0 001  \0 004  \0 001  \0 005 200 001
       00003   00001   00004   00001   00005  -32767
     [    3] [    1] [    4] [    1] [    5] [fill ]
\endcode

\section offset_format_spec The 64-bit Offset Format

The netCDF 64-bit offset format differs from the classic format only
in the VERSION byte, ‘\\x02’ instead of ‘\\x01’, and the OFFSET entity,
a 64-bit instead of a 32-bit offset from the beginning of the
file. This small format change permits much larger files, but there
are still some practical size restrictions. Each fixed-size variable
and the data for one record's worth of each record variable are still
limited in size to a little less that 4 GiB. The rationale for this
limitation is to permit aggregate access to all the data in a netCDF
variable (or a record's worth of data) on 32-bit platforms.

\section netcdf_4_spec The NetCDF-4 Format

The netCDF-4 format implements and expands the netCDF-3 data model by
using an enhanced version of HDF5 as the storage layer. Use is made of
features that are only available in HDF5 version 1.8 and later.

Using HDF5 as the underlying storage layer, netCDF-4 files remove many
of the restrictions for classic and 64-bit offset files. The richer
enhanced model supports user-defined types and data structures,
hierarchical scoping of names using groups, additional primitive types
including strings, larger variable sizes, and multiple unlimited
dimensions. The underlying HDF5 storage layer also supports
per-variable compression, multidimensional tiling, and efficient
dynamic schema changes, so that data need not be copied when adding
new variables to the file schema.

Creating a netCDF-4/HDF5 file with netCDF-4 results in an HDF5
file. The features of netCDF-4 are a subset of the features of HDF5,
so the resulting file can be used by any existing HDF5 application.

Although every file in netCDF-4 format is an HDF5 file, there are HDF5
files that are not netCDF-4 format files, because the netCDF-4 format
intentionally uses a limited subset of the HDF5 data model and file
format features. Some HDF5 features not supported in the netCDF
enhanced model and netCDF-4 format include non-hierarchical group
structures, HDF5 reference types, multiple links to a data object,
user-defined atomic data types, stored property lists, more permissive
rules for data object names, the HDF5 date/time type, and attributes
associated with user-defined types.

A complete specification of HDF5 files is beyond the scope of this
document. For more information about HDF5, see the HDF5 web site:
http://hdf.ncsa.uiuc.edu/HDF5/.

The specification that follows is sufficient to allow HDF5 users to
create files that will be accessible from netCDF-4.

\subsection creation_order Creation Order

The netCDF API maintains the creation order of objects that are
created in the file. The same is not true in HDF5, which maintains the
objects in alphabetical order. Starting in version 1.8 of HDF5, the
ability to maintain creation order was added. This must be explicitly
turned on in the HDF5 data file in several ways.

Each group must have link and attribute creation order set. The
following code (from libsrc4/nc4hdf.c) shows how the netCDF-4 library
sets these when creating a group.

\code
           /* Create group, with link_creation_order set in the group
            * creation property list. */
           if ((gcpl_id = H5Pcreate(H5P_GROUP_CREATE)) < 0)
              return NC_EHDFERR;
           if (H5Pset_link_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED|H5P_CRT_ORDER_INDEXED) < 0)
              BAIL(NC_EHDFERR);
           if (H5Pset_attr_creation_order(gcpl_id, H5P_CRT_ORDER_TRACKED|H5P_CRT_ORDER_INDEXED) < 0)
              BAIL(NC_EHDFERR);
           if ((grp->hdf_grpid = H5Gcreate2(grp->parent->hdf_grpid, grp->name,
                                            H5P_DEFAULT, gcpl_id, H5P_DEFAULT)) < 0)
              BAIL(NC_EHDFERR);
           if (H5Pclose(gcpl_id) < 0)
              BAIL(NC_EHDFERR);
\endcode

Each dataset in the HDF5 file must be created with a property list for
which the attribute creation order has been set to creation
ordering. The H5Pset_attr_creation_order function is used to set the
creation ordering of attributes of a variable.

The following example code (from libsrc4/nc4hdf.c) shows how the
creation ordering is turned on by the netCDF library.

\code
        /* Turn on creation order tracking. */
        if (H5Pset_attr_creation_order(plistid, H5P_CRT_ORDER_TRACKED|
                                       H5P_CRT_ORDER_INDEXED) < 0)
           BAIL(NC_EHDFERR);
\endcode

\subsection groups_spec Groups

NetCDF-4 groups are the same as HDF5 groups, but groups in a netCDF-4
file must be strictly hierarchical. In general, HDF5 permits
non-hierarchical structuring of groups (for example, a group that is
its own grandparent). These non-hierarchical relationships are not
allowed in netCDF-4 files.

In the netCDF API, the global attribute becomes a group-level
attribute. That is, each group may have its own global attributes.

The root group of a file is named “/” in the netCDF API, where names
of groups are used. It should be noted that the netCDF API (like the
HDF5 API) makes little use of names, and refers to entities by number.

\subsection dims_spec Dimensions with HDF5 Dimension Scales

Until version 1.8, HDF5 did not have any capability to represent
shared dimensions. With the 1.8 release, HDF5 introduced the dimension
scale feature to allow shared dimensions in HDF5 files.

The dimension scale is unfortunately not exactly equivalent to the
netCDF shared dimension, and this leads to a number of compromises in
the design of netCDF-4.

A netCDF shared dimension consists solely of a length and a name. An
HDF5 dimension scale also includes values for each point along the
dimension, information that is (optionally) included in a netCDF
coordinate variable.

To handle the case of a netCDF dimension without a coordinate
variable, netCDF-4 creates dimension scales of type char, and leaves
the contents of the dimension scale empty. Only the name and length of
the scale are significant. To distinguish this case, netCDF-4 takes
advantage of the NAME attribute of the dimension scale. (Not to be
confused with the name of the scale itself.) In the case of dimensions
without coordinate data, the HDF5 dimension scale NAME attribute is
set to the string: "This is a netCDF dimension but not a netCDF
variable."

In the case where a coordinate variable is defined for a dimension,
the HDF5 dimscale matches the type of the netCDF coordinate variable,
and contains the coordinate data.

A further difficulty arrises when an n-dimensional coordinate variable
is defined, where n is greater than one. NetCDF allows such coordinate
variables, but the HDF5 model does not allow dimension scales to be
attached to other dimension scales, making it impossible to completely
represent the multi-dimensional coordinate variables of the netCDF
model.

To capture this information, multidimensional coordinate variables
have an attribute named _Netcdf4Coordinates. The attribute is an array
of H5T_NATIVE_INT, with the netCDF dimension IDs of each of its
dimensions.

The _Netcdf4Coordinates attribute is otherwise hidden by the netCDF
API. It does not appear as one of the attributes for the netCDF
variable involved, except through the HDF5 API.

\subsection dim_spec2 Dimensions without HDF5 Dimension Scales

Starting with the netCDF-4.1 release, netCDF can read HDF5 files which
do not use dimension scales. In this case the netCDF library assigns
dimensions to the HDF5 dataset as needed, based on the length of the
dimension.

When an HDF5 file is opened, each dataset is examined in turn. The
lengths of all the dimensions involved in the shape of the dataset are
determined. Each new (i.e. previously unencountered) length results in
the creation of a phony dimension in the netCDF API.

This will not accurately detect a shared, unlimited dimension in the
HDF5 file, if different datasets have different lengths along this
dimension (possible in HDF5, but not in netCDF).

Note that this is a read-only capability for the netCDF library. When
the netCDF library writes HDF5 files, they always use a dimension
scale for every dimension.

Datasets must have either dimension scales for every dimension, or no
dimension scales at all. Partial dimension scales are not, at this
time, understood by the netCDF library.

\subsection dim_spec3 Dimension and Coordinate Variable Ordering

In order to preserve creation order, the netCDF-4 library writes
variables in their creation order. Since some variables are also
dimension scales, their order reflects both the order of the
dimensions and the order of the coordinate variables.

However, these may be different. Consider the following code:

\code
           /* Create a test file. */
           if (nc_create(FILE_NAME, NC_CLASSIC_MODEL|NC_NETCDF4, &ncid)) ERR;

           /* Define dimensions in order. */
           if (nc_def_dim(ncid, DIM0, NC_UNLIMITED, &dimids[0])) ERR;
           if (nc_def_dim(ncid, DIM1, 4, &dimids[1])) ERR;

           /* Define coordinate variables in a different order. */
           if (nc_def_var(ncid, DIM1, NC_DOUBLE, 1, &dimids[1], &varid[1])) ERR;
           if (nc_def_var(ncid, DIM0, NC_DOUBLE, 1, &dimids[0], &varid[0])) ERR;
\endcode

In this case the order of the coordinate variables will be different
from the order of the dimensions.

In practice, this should make little difference in user code, but if
the user is writing code that depends on the ordering of dimensions,
the netCDF library was updated in version 4.1 to detect this
condition, and add the attribute _Netcdf4Dimid to the dimension scales
in the HDF5 file. This attribute holds a scalar H5T_NATIVE_INT which
is the (zero-based) dimension ID for this dimension.

If this attribute is present on any dimension scale, it must be
present on all dimension scales in the file.

\subsection vars_spec Variables

Variables in netCDF-4/HDF5 files exactly correspond to HDF5
datasets. The data types match naturally between netCDF and HDF5.

In netCDF classic format, the problem of endianness is solved by
writing all data in big-endian order. The HDF5 library allows data to
be written as either big or little endian, and automatically reorders
the data when it is read, if necessary.

By default, netCDF uses the native types on the machine which writes
the data. Users may change the endianness of a variable (before any
data are written). In that case the specified endian type will be used
in HDF5 (for example, a H5T_STD_I16LE will be used for NC_SHORT, if
little-endian has been specified for that variable.)
- NC_BYTE = H5T_NATIVE_SCHAR
- NC_UBYTE = H5T_NATIVE_UCHAR
- NC_CHAR = H5T_C_S1
- NC_STRING = variable length array of H5T_C_S1
- NC_SHORT = H5T_NATIVE_SHORT
- NC_USHORT = H5T_NATIVE_USHORT
- NC_INT = H5T_NATIVE_INT
- NC_UINT = H5T_NATIVE_UINT
- NC_INT64 = H5T_NATIVE_LLONG
- NC_UINT64 = H5T_NATIVE_ULLONG
- NC_FLOAT = H5T_NATIVE_FLOAT
- NC_DOUBLE = H5T_NATIVE_DOUBLE

The NC_CHAR type represents a single character, and the NC_STRING an
array of characters. This can be confusing because a one-dimensional
array of NC_CHAR is used to represent a string (i.e. a scalar
NC_STRING).

An odd case may arise in which the user defines a variable with the
same name as a dimension, but which is not intended to be the
coordinate variable for that dimension. In this case the string
"_nc4_non_coord_" is pre-pended to the name of the HDF5 dataset, and
stripped from the name for the netCDF API.

\subsection atts_spec Attributes

Attributes in HDF5 and netCDF-4 correspond very closely. Each
attribute in an HDF5 file is represented as an attribute in the
netCDF-4 file, with the exception of the attributes below, which are
hidden by the netCDF-4 API.
- _Netcdf4Coordinates An integer array containing the dimension IDs of
  a variable which is a multi-dimensional coordinate variable.
- _nc3_strict When this (scalar, H5T_NATIVE_INT) attribute exists in
  the root group of the HDF5 file, the netCDF API will enforce the
  netCDF classic model on the data file.
- REFERENCE_LIST This attribute is created and maintained by the HDF5
  dimension scale API.
- CLASS This attribute is created and maintained by the HDF5 dimension
  scale API.
- DIMENSION_LIST This attribute is created and maintained by the HDF5
  dimension scale API.
- NAME This attribute is created and maintained by the HDF5 dimension
  scale API.
- _Netcdf4Dimid Holds a scalar H5T_NATIVE_INT that is the (zero-based)
  dimension ID for this dimension, needed when dimensions and
  coordinate variables are defined in different orders.
- _NCProperties Holds provenance information about a file at the time
  it was created. It specifies the versions of the netCDF and HDF5
  libraries used to create the file.

\subsection user_defined_spec User-Defined Data Types

Each user-defined data type in an HDF5 file exactly corresponds to a
user-defined data type in the netCDF-4 file. Only base data types
which correspond to netCDF-4 data types may be used. (For example, no
HDF5 reference data types may be used.)

\subsection compression_spec Compression

The HDF5 library provides data compression using the zlib library and
the szlib library. NetCDF-4 only allows users to create data with the
zlib library (due to licensing restrictions on the szlib
library). Since HDF5 supports the transparent reading of the data with
either compression filter, the netCDF-4 library can read data
compressed with szlib (if the underlying HDF5 library is built to
support szlib), but has no way to write data with szlib compression.

With zlib compression (a.k.a. deflation) the user may set a deflation
factor from 0 to 9. In our measurements the zero deflation level does
not compress the data, but does incur the performance penalty of
compressing the data. The netCDF API does not allow the user to write
a variable with zlib deflation of 0 - when asked to do so, it turns
off deflation for the variable instead. NetCDF can read an HDF5 file
with deflation of zero, and correctly report that to the user.

\section netcdf_4_classic_spec The NetCDF-4 Classic Model Format

Every classic and 64-bit offset file can be represented as a netCDF-4
file, with no loss of information. There are some significant benefits
to using the simpler netCDF classic model with the netCDF-4 file
format. For example, software that writes or reads classic model data
can write or read netCDF-4 classic model format data by
recompiling/relinking to a netCDF-4 API library, with no or only
trivial changes needed to the program source code. The netCDF-4
classic model format supports this usage by enforcing rules on what
functions may be called to store data in the file, to make sure its
data can be read by older netCDF applications (when relinked to a
netCDF-4 library).

Writing data in this format prevents use of enhanced model features
such as groups, added primitive types not available in the classic
model, and user-defined types. However performance features of the
netCDF-4 formats that do not require additional features of the
enhanced model, such as per-variable compression and chunking,
efficient dynamic schema changes, and larger variable size limits,
offer potentially significant performance improvements to readers of
data stored in this format, without requiring program changes.

When a file is created via the netCDF API with a CLASSIC_MODEL mode
flag, the library creates an attribute (_nc3_strict) in the root
group. This attribute is hidden by the netCDF API, but is read when
the file is later opened, and used to ensure that no enhanced model
features are written to the file.

\section hdf4_sd_format HDF4 SD Format

Starting with version 4.1, the netCDF libraries can read HDF4 SD
(Scientific Dataset) files. Access is limited to those HDF4 files
created with the Scientific Dataset API. Access is read-only.

Dataset types are translated between HDF4 and netCDF in a
straightforward manner.
- DFNT_CHAR = NC_CHAR
- DFNT_UCHAR, DFNT_UINT8 = NC_UBYTE
- DFNT_INT8 = NC_BYTE
- DFNT_INT16 = NC_SHORT
- DFNT_UINT16 = NC_USHORT
- DFNT_INT32 = NC_INT
- DFNT_UINT32 = NC_UINT
- DFNT_FLOAT32 = NC_FLOAT
- DFNT_FLOAT64 = NC_DOUBLE
