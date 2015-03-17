#! /bin/sh
#
# This forms the basis for the nc-config utility, which tells you
# various things about the netCDF installation. This code was
# contributed by netCDF user Arlindo DaSilva. Thanks Arlindo!

prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=@CMAKE_INSTALL_PREFIX@
libdir=@CMAKE_INSTALL_PREFIX@/lib
includedir=@CMAKE_INSTALL_PREFIX@/include

cc="@CMAKE_C_COMPILER@"
cflags="-I@CMAKE_INSTALL_PREFIX@/include"
libs="-L@CMAKE_INSTALL_PREFIX@/lib @NC_LIBS@"

has_dap="@USE_DAP@"
if [ -z $has_dap ]; then
    has_dap="no"
else
    has_dap="yes"
fi

has_nc2="@BUILD_V2@"

if [ -z $has_nc2 ]; then
    has_nc2="no"
else
    has_nc2="yes"
fi

has_nc4="@USE_NETCDF4@"
if [ -z $has_nc4 ]; then
    has_nc4="no"
else
    has_nc4="yes"
fi

has_hdf4="@USE_HDF4@"
if [ -z $has_hdf4 ]; then
    has_hdf4="no"
else
    has_hdf4="yes"
fi

has_pnetcdf="@USE_PNETCDF@"
if [ -z $has_pnetcdf ]; then
    has_pnetcdf="no"
else
    has_pnetcdf="yes"
fi

has_hdf5="@USE_HDF5@"
if [ -z $has_hdf5 ]; then
    has_hdf5="no"
else
    has_hdf5="yes"
fi

has_szlib="@USE_SZLIB@"
if [ -z $has_szlib ]; then
    has_szlib="no"
else
    has_szlib="yes"
fi


version="@PACKAGE@ @VERSION@"

has_f90="no"
if type -p nf-config > /dev/null 2>&1; then
  fc=`nf-config --fc`
  fflags=`nf-config --fflags`
  flibs=`nf-config --flibs`
  has_f90=`nf-config --has-f90`
fi

has_cxx="no"
has_cxx4="no"
if type -p ncxx4-config > /dev/null 2>&1; then
  cxx4=`ncxx4-config --cxx`
  has_cxx4="yes"
elif type -p ncxx-config > /dev/null 2>&1; then
  cxx=`ncxx-config --cxx`
  has_cxx="yes"
fi

usage()
{
    cat <<EOF
Usage: nc-config [OPTION]

Available values for OPTION include:

  --help        display this help message and exit
  --all         display all options
  --cc          C compiler
  --cflags      pre-processor and compiler flags
  --has-dap     whether OPeNDAP is enabled in this build
  --has-nc2     whether NetCDF-2 API is enabled
  --has-nc4     whether NetCDF-4/HDF-5 is enabled in this build
  --has-hdf5    whether HDF5 is used in build (always the same as --has-nc4)
  --has-hdf4    whether HDF4 was used in build
  --has-pnetcdf whether parallel-netcdf (a.k.a. pnetcdf) was used in build
  --has-szlib   whether szlib is included in build
  --libs        library linking information for netcdf
  --prefix      Install prefix
  --includedir  Include directory
  --version     Library version

EOF
# When supported by ncxx4-config and ncxx-config, add
#  --cxxflags    flags needed to compile a netCDF-4 C++ program
#  --cxxlibs     libraries needed to link a netCDF-4 C++ program
if type -p ncxx4-config > /dev/null 2>&1; then
    cat <<EOF
  --cxx4         C++ compiler for netCDF-4 C++ library
  --has-c++4     whether netCDF-4 C++ API is installed
EOF
elif type -p ncxx-config > /dev/null 2>&1; then
    cat <<EOF
  --cxx         C++ compiler
  --has-c++     whether C++ API is installed

EOF
fi
if type -p nf-config > /dev/null 2>&1; then
    cat <<EOF
  --fc          Fortran compiler
  --fflags      flags needed to compile a Fortran program
  --flibs       libraries needed to link a Fortran program
  --has-f90     whether Fortran 90 API is installed

EOF
fi
    exit $1
}

all()
{
        echo
        echo "This $version has been built with the following features: "
        echo
        echo "  --cc        -> $cc"
        echo "  --cflags    -> $cflags"
        echo "  --libs      -> $libs"
        echo
        echo "  --has-c++   -> $has_cxx"
        echo "  --cxx       -> $cxx"
        echo "  --has-c++4  -> $has_cxx4"
        echo "  --cxx4      -> $cxx4"
        echo
        echo "  --fc        -> $fc"
        echo "  --fflags    -> $fflags"
        echo "  --flibs     -> $flibs"
        echo "  --has-f90   -> $has_f90"
        echo
        echo "  --has-dap   -> $has_dap"
        echo "  --has-nc2   -> $has_nc2"
        echo "  --has-nc4   -> $has_nc4"
        echo "  --has-hdf5  -> $has_hdf5"
        echo "  --has-hdf4  -> $has_hdf4"
        echo "  --has-pnetcdf-> $has_pnetcdf"
        echo "  --has-szlib -> $has_szlib"
	echo
        echo "  --prefix    -> $prefix"
        echo "  --includedir-> $includedir"
        echo "  --version   -> $version"
        echo
}

if test $# -eq 0; then
    usage 1
fi

while test $# -gt 0; do
    case "$1" in
    # this deals with options in the style
    # --option=value and extracts the value part
    # [not currently used]
    -*=*) value=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) value= ;;
    esac

    case "$1" in

    --help)
	usage 0
	;;

    --all)
	all
	;;

    --cc)
	echo $cc
	;;

    --cflags)
	echo $cflags
	;;

    --has-dap)
       	echo $has_dap
       	;;

    --has-nc2)
       	echo $has_nc2
       	;;

    --has-nc4)
       	echo $has_nc4
       	;;

    --has-hdf5)
       	echo $has_hdf5
       	;;

    --has-hdf4)
       	echo $has_hdf4
       	;;

    --has-pnetcdf)
       	echo $has_pnetcdf
       	;;

    --has-szlib)
       	echo $has_szlib
       	;;

     --libs)
       	echo $libs
       	;;

    --prefix)
       	echo "${prefix}"
       	;;

    --includedir)
       	echo "${includedir}"
       	;;

    --version)
	echo $version
	;;

    --has-c++)
       	echo $has_cxx
       	;;

    --cxx)
	echo $cxx
	;;

    --has-c++4)
       	echo $has_cxx4
       	;;

    --cxx4)
	echo $cxx4
	;;

#    --cxxflags)
#	echo $cxxflags
#	;;
#
#    --cxxlibs)
#	echo $cxxlibs
#	;;

    --fc)
	echo $fc
	;;

    --fflags)
	echo $fflags
	;;

    --flibs)
       	echo $flibs
       	;;

    --has-f90)
       	echo $has_f90
       	;;

    *)
        echo "unknown option: $1"
	usage
	exit 1
	;;
    esac
    shift
done

exit 0
