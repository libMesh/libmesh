## IMPORTANT -- READ BEFORE BUILDING!
Do not download a GitHub-generated "ZIP" archive. These do not contain the required submodules, and therefore cannot be used to build libmesh. Use only git clones or "release" tarballs when following these instructions.

## Build Instructions
The default is to build libmesh "out of tree," i.e. within a separate `build` directory, rather than in the source tree itself. This simplifies the process of having multiple, independently-configured builds.
1. `cd` to location of libmesh clone or extracted release tarball.
1. (Only if using a git clone) `git submodule update --init --recursive`
1. `mkdir build`
1. `cd build`
1. `../configure --prefix=/path/to/libmesh/install`
1. `make`
1. `make check` (optional, runs the example programs and unit tests when possible)
1. `make install`


## METHODS
libMesh supports the notion of multiple methods, that is, configuration
settings used to build the library.  The major methods supported by
the library are:

* opt: Fully optimized mode, with little to no error checking. No debugging
  symbols are included in the resulting library.  Aggressive optimization
  flags are used.

* dbg: Full debugging mode. All useful compiler warnings are enabled,
  as well as robust internal state checking. The asymptotic complexity
  of some algorithms is allowed to be worse than the design spec states.

* devel: Use high levels of compiler optimization, but also enable internal
  state checking.  Debugging symbols are included, but the resulting
  code is not always easy to navigate in a debugger because of
  compiler optimizations.

* pro: Optimized code path with compiler flags suitable for use with gprof.

* oprof: Optimized code path with compiler flags suitable for use with oprofile.

To select a set of methods, you can pass them to configure in one of two ways:

    ../configure --with-methods="opt dbg devel"

or

    ../configure METHODS="devel oprof"

If unspecified, `METHODS="opt dbg devel"` is the default.


## Multiple Builds with Different Compilers
libMesh fully supports out-of-tree builds, and users are encouraged to use this
feature when needed to support multiple compilers. For example, on a system
where multiple compilers are available and accessible via modules, you can share
the same source tree by creating a subdirectory for each compiler build:

    export LIBMESH_SRC=/local/libmesh
    cd $LIBMESH_SRC
    module load gcc/9.3
    cd $LIBMESH_SRC && mkdir gcc-9.3 && cd gcc-9.3 && ../configure && make && make install
    module swap gcc/9.3 intel/19.0
    cd $LIBMESH_SRC && mkdir intel-19.0 && cd intel-19.0 && ../configure && make && make install


## Dependencies
libMesh has no required dependencies other than a C++ compiler which
fully supports the C++11 standard. To run on distributed memory
platforms in parallel, you will also need MPI.

## Optional Packages
We support a [large number](http://libmesh.github.io/externalsoftware.html) of
optional packages. Some of these are distributed inside the `contrib` directory
and compiled directly with libMesh, others can be used via third-party installations.

## License
LibMesh is open source software distributed under the LGPL license, version 2.1.
