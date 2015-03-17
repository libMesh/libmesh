# Cmake Package & Configuration file support.
# Based on code from 'Mastering Cmake'

# Compute installation prefix relative to this file.
get_filename_component (_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component (_prefix "${_dir}/../.." ABSOLUTE)

# Import the targets
include ("${_prefix}/share/cmake/netcdf-targets.cmake")

# Report other information
set (netcdf_INCLUDE_DIRS "${_prefix}/include/")



