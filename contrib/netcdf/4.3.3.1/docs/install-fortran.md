Building the NetCDF-4.2 and later Fortran libraries {#building_netcdf_fortran}
===================================================

[TOC]

In versions before 4.2, the Fortran netCDF library source was bundled
with the C library source in one distribution, and it was possible to
combine the libraries in a single library file. With version 4.2, the
Fortran netCDF library for Fortran77 and Fortran90 APIs has been
separated into its own source distribution, and should now be built as a
separate library, after the C library is built and installed. This
separation simplifies the building and use of the C and Fortran netCDF
libraries and allows them to evolve independently.

Please note that in the example commands below, we assume use of a
POSIX-standard shell, such as sh, bash, ksh, or zsh. If you are using
csh instead, you will have to use the

       setenv ENV_VARIABLE  value

syntax to set environment variables instead of the

       ENV_VARIABLE=value

syntax used in the examples that use a POSIX-standard shell. In either
case, <I>${DIR1}</I> is the value of the environment variable <I>DIR1</I>.

It will be easier to build the netCDF Fortran library if the C (and if
needed, HDF5) libraries are built as shared libraries (the default), but
you can also use static libraries, as described in a later section.

Building with shared libraries {#building_fortran_shared_libraries}
==============================

1.  First make sure the netCDF C library has been built, tested, and
    installed under directory <I>${DIR1}</I>, as specified by
    --prefix=<I>${DIR1}</I> to the C library configure script, or under
    directory /usr/local by default.
2.  For the Fortran netCDF library, use the same C compiler as used to
    create the netCDF C library, specified with the CC environment
    variable, if necessary.
3.  If the netCDF C library was installed as a shared library in a
    location that is not searched by default, you will need to set the
    LD\_LIBRARY\_PATH environment variable (or DYLD\_LIBRARY\_PATH on
    OSX) to specify that directory before running the configure script,
    for example:

        export LD_LIBRARY_PATH=${DIR1}/lib:${LD_LIBRARY_PATH}

4.  If you set the LD\_LIBRARY\_PATH (or DYLD\_LIBRARY\_PATH)
    environment variable in the previous step, don't use the "sudo"
    command before the following "configure" or "make check" commands.
    Using "sudo" causes the LD\_\* environment variables to be ignored,
    as a security precaution. You can use "sudo make install" as the
    last step, but you shouldn't need to use "sudo" before that.
5.  For the configure script, set CPPFLAGS and LDFLAGS variables to
    specify the include and lib directories for the netCDF C library.
    For example, to install the Fortran libraries in the same directory
    *\${DIR1}* where the C netCDF library is installed:

          CPPFLAGS=-I${DIR1}/include LDFLAGS=-L${DIR1}/lib ./configure --prefix=${DIR1}

    If you are cross-compiling, you should also include the configure
    option "--disable-fortran-type-check", as in:

          CPPFLAGS=-I${DIR1}/include LDFLAGS=-L${DIR1}/lib \
            ./configure --disable-fortran-type-check --prefix=${DIR1}

6.  If that succeeds, run "make check".
7.  If that succeeds, run "make install" or "sudo make install".

Building with static libraries {#building_fortran_with_static_libraries}
==============================

If you can't build the C netCDF library as a shared library or if it has
already been installed by someone else only as a static library (which
means there are no libnetcdf.so files in the library directory where the
netCDF C library was installed), then building and installing the
Fortran netCDF libraries will be somewhat more complicated.

If you need to set the LD\_LIBRARY\_PATH (or DYLD\_LIBRARY\_PATH)
environment variable, don't use the "sudo" command before the following
"configure" or "make check" commands. Using "sudo" causes the LD\_\*
environment variables to be ignored. You can use "sudo make install" as
the last step, but you shouldn't need to use "sudo" before that.

1.  Assume the static netCDF C library is installed under *\${DIR1}*,
    and the other needed shared libraries for HDF5, zlib, and curl are
    installed under *\${DIR2}* (which might be the same as *\${DIR1}*).
2.  Use the same C compiler as used to create the netCDF C library,
    specified with the CC environment variable, if necessary.
3.  Set the CPPFLAGS, LDFLAGS, and LD\_LIBRARY\_PATH environment
    variables to specify where the netCDF C library is installed and
    where the other shared libraries may be found, before running the
    configure script. For example:

          CPPFLAGS="-I${DIR1}/include -I${DIR2}/include" \
          LD_LIBRARY_PATH=${DIR1}/lib:${DIR2}/lib:${LD_LIBRARY_PATH} \
          LDFLAGS="-L${DIR1}/lib -L${DIR2}/lib" \
          LIBS="-lnetcdf -lhdf5_hl -lhdf5 -lz -lcurl" \
          ./configure --disable-shared --prefix=${DIR1}

    If you are cross-compiling, you should also include the configure
    option "--disable-fortran-type-check".

4.  For parallel I/O: The configure script sets CFLAGS appropriately for
    standard compilers, but if you are building with parallel I/O using
    wrappers such as mpicc and mpif90, you sometimes have to set CFLAGS
    to indicate which Fortran compiler is wrapped by mpif90. For
    example, if "mpicc --show" and "mpif90 --show" indicate gcc and
    gfortran are being used, then set CFLAGS=-DgFortran, and similarly
    set CFLAGS=-DpgiFortran for Portland Group compilers.
5.  If that succeeds, run "make check".
6.  If that succeeds, run "make install" or "sudo make install".

Linking your programs with netCDF Fortran libraries {#linking_against_netcdf_fortran}
==============================

If you built the shared libraries, you can link with something like

       fortran_compiler my_prog.f -o my_prog -I${DIR1}/include -L${DIR1}/lib -lnetcdff

to link your Fortran software with the installed netCDF Fortran and C
libraries. If you didn't install the shared libraries in a standard
place, you may need to set LD\_LIBRARY\_PATH (or DYLD\_LIBRARY\_PATH for
OSX) before running the resulting program.

If you built static libraries, you will need to use something like

       fortran_compiler my_prog.f -o my_prog -I${DIR1}/include \
        -L${DIR1}/lib -lnetcdff -lnetcdf -L${DIR2}/lib -lhdf5_hl -lhdf5 -lz -lcurl -lm

to link Fortran software with the installed Fortran library and the
libraries on which it depends.

A simpler alternative that should work for either shared or static
libraries is to use the "nf-config" utility installed in *${DIR1}*/bin:

       `nf-config --fc` my_prog.f -o my_prog `nf-config --fflags --flibs`

or the more general "pkg-config" utility, if you have it:

       fortran_compiler my_prog.f -o my_prog `pkg-config --cflags --libs netcdf-fortran`
       

Specifying The Environment for Building {#specify_build_env_fortran}
========================================


The netCDF configure script searches your path to find the compilers and tools it needed. To use compilers that can't be found in your path, set their environment variables.

The configure script will use gcc and associated GNU tools if they are found. Many users, especially those with performance concerns, will wish to use a vendor supplied compiler.

For example, on an AIX system, users may wish to use xlc (the AIX compiler) in one of its many flavors. Set environment variables before the build to achieve this.

For example, to change the C compiler, set CC to xlc (in sh: export CC=xlc). (But don't forget to also set CXX to xlC, or else configure will try to use g++, the GNU C++ compiler to build the netCDF C++ API. Similarly set FC to xlf90 so that the Fortran APIs are built properly.)

By default, the netCDF library is built with assertions turned on. If you wish to turn off assertions, set CPPFLAGS to -DNDEBUG (csh ex: setenv CPPFLAGS -DNDEBUG).

If GNU compilers are used, the configure script sets CPPFLAGS to “-g -O2”. If this is not desired, set CPPFLAGS to nothing, or to whatever other value you wish to use, before running configure.

For cross-compiles, the following environment variables can be used to override the default fortran/C type settings like this (in sh):

     export NCBYTE_T=''integer(selected_int_kind(2))''
     export NCSHORT_T=''integer*2''
     export NF_INT1_T=''integer(selected_int_kind(2))''
     export NF_INT2_T=''integer*2''
     export NF_INT1_IS_C_SHORT=1
     export NF_INT2_IS_C_SHORT=1
     export NF_INT_IS_C_INT=1
     export NF_REAL_IS_C_FLOAT=1
     export NF_DOUBLEPRECISION_IS_C_DOUBLE=1
     
In this case you will need to run configure with –disable-fortran-compiler-check and –disable-fortran-type-check.

Variable Description Notes
--------------------------

Variable | Usage | Description
---|---|---
CC	| C compiler	| If you don't specify this, the configure script will try to find a suitable C compiler. The default choice is gcc. If you wish to use a vendor compiler you must set CC to that compiler, and set other environment variables (as described below) to appropriate settings.
FC	| Fortran compiler (if any)| 	If you don't specify this, the configure script will try to find a suitable Fortran and Fortran 77 compiler. Set FC to "" explicitly, or provide the –disable-f77 option to configure, if no Fortran interface (neither F90 nor F77) is desired. Use –disable-f90 to disable the netCDF Fortran 90 API, but build the netCDF Fortran 77 API.
F77	| Fortran 77 compiler (if any)	| Only specify this if your platform explicitly needs a different Fortran 77 compiler. Otherwise use FC to specify the Fortran compiler. If you don't specify this, the configure script will try to find a suitable Fortran compiler. For vendor compilers, make sure you're using the same vendor's Fortran 90 compiler. Using Fortran compilers from different vendors, or mixing vendor compilers with g77, the GNU F77 compiler, is not supported and may not work.
CXX	| C++ compiler	| If you don't specify this, the configure script will try to find a suitable C++ compiler. Set CXX to "" explicitly, or use the –disable-cxx configure option, if no C++ interface is desired. If using a vendor C++ compiler, use that vendor's C compiler to compile the C interface. Using different vendor compilers for C and C++ may not work.
CFLAGS	| C compiler flags	| "-O" or "-g", for example.
CPPFLAGS	| C preprocessor options	| "-DNDEBUG" to omit assertion checks, for example.
FCFLAGS| 	Fortran 90 compiler flags	| "-O" or "-g", for example. These flags will be used for FORTRAN 90. If setting these you may also need to set FFLAGS for the FORTRAN 77 test programs.
FFLAGS	| Fortran 77 compiler flags	| "-O" or "-g", for example. If you need to pass the same arguments to the FORTRAN 90 build, also set FCFLAGS.
CXXFLAGS	| C++ compiler flags	| "-O" or "-g", for example.
ARFLAGS, NMFLAGS, FPP, M4FLAGS, LIBS, FLIBS, FLDFLAGS	| Miscellaneous	| One or more of these were needed for some platforms, as specified below. Unless specified, you should not set these environment variables, because that may interfere with the configure script. 