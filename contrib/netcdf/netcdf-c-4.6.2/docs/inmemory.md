NetCDF In-Memory Support
====================================

<!-- double header is needed to workaround doxygen bug -->

NetCDF In-Memory Support {#inmemory}
====================================

[TOC]

Introduction {#inmemory_intro}
--------------

It can be convenient to operate on a netcdf file whose
content is held in memory instead of in a disk file.
The netcdf API has been modified in a number of ways
to support this capability.

Actually, three distinct but related capabilities are provided.

1. DISKLESS -- Read a file into memory, operate on it, and optionally
write it back out to disk when nc_close() is called.
2. INMEMORY -- Tell the netcdf-c library to treat a provided block
of memory as if it were a netcdf file. At close, it is possible to ask
for the final contents of the memory chunk. Be warned that there is
some complexity to this as described below.
4. MMAP -- (deprecated) Tell the netcdf-c library to use the *mmap()* operating
system functionality to access a file.

The first two capabilities are intertwined in the sense that the
*diskless* capability makes use internally of the *inmemory*
capability (for netcdf classic only). But, the *inmemory*
capability can be used independently of the *diskless*
capability.

The *mmap()* capability provides a capability similar to *diskless* but
using special capabilities of the underlying operating system. It turns out
that the mmap capability has seen no significant use, so its use is deprecated
and will be removed at some point in the future.

Note also that *diskless* and *inmemory* can be used for both
*netcdf-3* (classic) and *netcdf-4* (enhanced) data. The *mmap*
capability can only be used with *netcdf-3*.

Enabling Diskless File Access {#Enable_Diskless}
--------------
The *diskless* capability can be used relatively transparently
using the *NC_DISKLESS* mode flag.

Note that since the file is stored in memory, size limitations apply.
If you are on using a 32-bit pointer then the file size must be less than 2^32
bytes in length. On a 64-bit machine, the size must be less than 2^64 bytes.

Also note that for a diskless file, there are two notions of
*write* with respect to the file. The first notion is that the
file is read-only through the netCDF API. For example, if the file
is read-only, then a call to, for example, _nc_def_dim()_ will fail.
The second notion of *write* refers to the file on disk to which 
the contents of memory might be persisted.

WARNING: control of the two kinds of *write* has changed since
release 4.6.1.

The mode flag NC_WRITE determines the first kind of *write*.
If set, then NC_WRITE means that the file can be modified through
the netCDF API, otherwise it is read-only. This is a change since
release 4.6.1.

The new mode flag NC_PERSIST now determines the second kind of
*write*.  If set, then NC_PERSIST means that the memory contents
will be persisted to disk, possibly overwriting the previous
file contents.  Otherwise, the default is to throw away the
in-memory contents.

### Diskless File Open
Calling *nc_open()* using the mode flag *NC_DISKLESS* will cause
the file being opened to be read into memory. When calling *nc_close()*,
the file will optionally be re-written (aka "persisted") to disk. This
persist capability will be invoked if and only if *NC_PERSIST* is specified
in the mode flags at the call to *nc_open()*.

### Diskless File Create
Calling *nc_create()* using the mode flag *NC_DISKLESS* will cause
the file to initially be created and kept in memory.
When calling *nc_close()*, the file will be written
to disk if and only if *NC_PERSIST* is specified
in the mode flags at the call to *nc_create()*.

Enabling Inmemory File Access {#Enable_Inmemory}
--------------

The netcdf API has been extended to support the inmemory capability.
The relevant API is defined in the file `netcdf_mem.h`.

The important data structure to use is `NC_memio`.
````
typedef struct NC_memio {
    size_t size;
    void* memory;
    int flags;
} NC_memio;

````
An instance of this data structure is used when providing or
retrieving a block of data. It specifies the memory and its size
and also some relevant flags that define how to manage the memory.

Current only one flag is defined -- *NC_MEMIO_LOCKED*.
This tells the netcdf library that it should never try to
*realloc()* the memory nor to *free()* the memory. Note
that this does not mean that the memory cannot be modified, but
only that the modifications will be within the confines of the provided
memory. If doing such modifications is impossible without
reallocating the memory, then the modification will fail.

### In-Memory API

The new API consists of the following functions.
````
int nc_open_mem(const char* path, int mode, size_t size, void* memory, int* ncidp);

int nc_create_mem(const char* path, int mode, size_t initialsize, int* ncidp);

int nc_open_memio(const char* path, int mode, NC_memio* info, int* ncidp);

int nc_close_memio(int ncid, NC_memio* info);

````
### The **nc_open_mem** Function

The *nc_open_mem()* function is actually a convenience
function that internally invokes *nc_open_memio()*.
It essentially provides simple read-only access to a chunk of memory
of some specified size.

### The **nc_open_memio** Function

This function provides a more general read/write capability with respect
to a chunk of memory. It has a number of constraints and its
semantics are somewhat complex. This is primarily due to limitations
imposed by the underlying HDF5 library.

The constraints are as follows.

1. If the *NC_MEMIO_LOCKED* flag is set, then the netcdf library will
make no attempt to reallocate or free the provided memory.
If the caller invokes the *nc_close_memio()* function to retrieve the
final memory block, it should be the same
memory block as was provided when *nc_open_memio* was called.
Note that it is still possible to modify the in-memory file if the NC_WRITE
mode flag was set. However, failures can occur if an operation
cannot complete because the memory needs to be expanded.
2. If the *NC_MEMIO_LOCKED* flag is <b>not</b> set, then
the netcdf library will take control of the incoming memory.
This means that the user should not make any attempt to free
or even read the incoming memory block in this case.
The newcdf library is free to reallocate the incomming
memory block to obtain a larger block when an attempt to modify
the in-memory file requires more space. Note that implicit in this
is that the old block -- the one originally provided -- may be
free'd as a side effect of re-allocating the memory using the
*realloc()* function.
The caller may invoke the *nc_close_memio()* function to retrieve the
final memory block, which may not be the same as the originally block
provided by the caller. In any case, the returned block must always be freed
by the caller and the original block should not be freed.

### The **nc_create_mem** Function

This function allows a user to create an in-memory file, write to it,
and then retrieve the final memory using *nc_close_memio()*.
The *initialsize* argument to *nc_create_mem()* tells the library
how much initial memory to allocate. Technically, this is advisory only
because it may be ignored by the underlying HDF5 library.
It is used, however, for netcdf-3 files. 

### The **nc_close_memio** Function

The ordinary *nc_close()* function can be called to close an in-memory file.
However, it is often desirable to obtain the final size and memory block
for the in-memory file when that file has been modified.
The *nc_close_memio()* function provides a means to do this.
Its second argument is a pointer to an *NC_memio* object
into which the final memory and size are stored. WARNING,
the returned memory is owned by the caller and so the caller
is responsible for calling *free()* on that returned memory.

### Support for Writing with *NC_MEMIO_LOCKED*

When the NC_MEMIO_LOCKED flag is set in the *NC_memio* object
passed to *nc_open_memio()*, it is still possible to modify
the opened in-memory file (using the NC_WRITE mode flag).

The big problem is that any changes must fit into the memory provided
by the caller via the *NC_memio* object. This problem can be
mitigated, however, by using the "trick" of overallocating
the caller supplied memory. That is, if the original file is, say, 300 bytes,
then it is possible to allocate, say, 65000 bytes and copy the original file
into the first 300 bytes of the larger memory block. This will allow
the netcdf-c library to add to the file up to that 65000 byte limit.
In this way, it is possible to avoid memory reallocation while still
allowing modifications to the file. You will still need to call
*nc_close_memio()* to obtain the size of the final, modified, file.

Enabling MMAP File Access (Deprecated) {#Enable_MMAP}
--------------

The MMAP functionality is deprecated.

Some operating systems provide a capability called MMAP.
This allows disk files to automatically be mapped to chunks of memory.
It operates in a fashion somewhat similar to operating system virtual
memory, except with respect to a file.

By setting mode flag NC_MMAP, it is possible to do the equivalent
of NC_DISKLESS but using the operating system's mmap capabilities.

Currently, MMAP support is only available when using netcdf-3 or cdf5
files.

Known Bugs {#Inmemory_Bugs}
--------------

1. If you are modifying a locked memory chunk (using
   NC_MEMIO_LOCKED) and are accessing it as a netcdf-4 file, and
   you overrun the available space, then the HDF5 library will
   fail with a segmentation fault.

2. You will get an HDF5 error under the following conditions.

   1. You call nc_open on a file with the flags NC_DISKLESS|NC_WRITE
      but without NC_PERSIST.
   2. The file to be read is read-only (i.e. mode 0444).

   Note that this should be ok because the modifications to the file
   are not intended to pushed back into the disk file. However, the
   HDF5 core driver does not allow this.

References {#Inmemory_References}
--------------

1. https://support.hdfgroup.org/HDF5/doc1.8/Advanced/FileImageOperations/HDF5FileImageOperations.pdf

Point of Contact
--------------

__Author__: Dennis Heimbigner<br>
__Email__: dmh at ucar dot edu
__Initial Version__: 2/3/2018<br>
__Last Revised__: 2/5/2018


