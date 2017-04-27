# NetCDF Programming Notes {#programming_notes}

[TOC]

<H2>See Also:</H2>

* \subpage nc-error-codes

# Ignored if NULL {#ignored_if_null}

Many of the argurments of netCDF functions are pointers. For example,
the nc_inq() functions takes four pointers:

~~~.C
int nc_inq(int ncid, int *ndimsp, int *nvarsp, int *nattsp, int *unlimdimidp);
~~~

A NULL may be passed for any of these pointers, and it will be ignored. For example, interested in the number of dimensions only, the following code will work:

~~~.C
int ndims;
...
if (nc_inq(ncid, &ndims, NULL, NULL, NULL))
   return SOME_ERROR;
~~~

# Allocating Storage for the Result {#allocating_storage_for_the_result}

User must allocate space for the result of an inq function before the function is called.

# Specify a Hyperslab {#specify_hyperslab}

The NetCDF allows specification of hyperslabs to be read or written
with vectors which specify the start, count, stride, and mapping.

## A Vector Specifying Start Index for Each Dimension {#start_vector}

A vector of size_t integers specifying the index in the
variable where the first of the data values will be read.

The indices are relative to 0, so for example, the first data value of
a variable would have index (0, 0, ... , 0).

The length of start vector must be the same as the number of
dimensions of the specified variable. The elements of start
correspond, in order, to the variable's dimensions.

## A Vector Specifying Count for Each Dimension {#count_vector}

A vector of size_t integers specifying the edge lengths
along each dimension of the block of data values to be read.

To read a single value, for example, specify count as (1, 1, ... , 1).

The length of count is the number of dimensions of the specified
variable. The elements of count correspond, in order, to the
variable's dimensions.

Setting any element of the count array to zero causes the function to
exit without error, and without doing anything.

## A Vector Specifying Stride for Each Dimension {#stride_vector}

A vector of size_t integers specifying the interval between selected
indices.

A value of 1 accesses adjacent values of the netCDF variable in the
corresponding dimension; a value of 2 accesses every other value of
the netCDF variable in the corresponding dimension; and so on.

The elements of the stride vector correspond, in order, to the
variable's dimensions.

A NULL stride argument is treated as (1, 1, ... , 1).

## A Vector Specifying Mapping for Each Dimension {#map_vector}

A vector of integers that specifies the mapping between the dimensions
of a netCDF variable and the in-memory structure of the internal data
array.

imap[0] gives the distance between elements of the internal array
corresponding to the most slowly varying dimension of the netCDF
variable. imap[n-1] (where n is the rank of the netCDF variable) gives
the distance between elements of the internal array corresponding to
the most rapidly varying dimension of the netCDF variable. Intervening
imap elements correspond to other dimensions of the netCDF variable in
the obvious way. Distances between elements are specified in
type-independent units of elements.

\note The distance between internal elements that occupy adjacent
memory locations is 1 and not the element's byte-length as in netCDF
2.

# NetCDF ID {#ncid}

Most netCDF function require the netCDF ID as a first parameter.

In the netCDF classic model, the netCDF ID is associated with an open
file. Each file, when opened by nc_open(), or created by nc_create(),
is assigned an ncid, which it retains until nc_close() is called.

In the netCDF enhanced model, the ncid refers to a group with a
file. (Each file contains at least the root group, which is the ncid
that is returned by nc_open() and nc_create().)

For netCDF-4/HDF5 files, netCDF IDs can come not just from nc_open()
and nc_create(), but also from nc_def_grp(), nc_inq_grps(),
nc_inq_ncid(), nc_inq_grp_parent(), and nc_inq_grp_full_ncid().
