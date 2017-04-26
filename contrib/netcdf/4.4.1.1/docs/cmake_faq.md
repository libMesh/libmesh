CMake-Related Frequently-Asked-Questions (FAQ) {#cmake_faq}
==============================================

Below are a list of commonly-asked questions regarding NetCDF and CMake.

## How can I see the options available to CMake? {#listoptions}

        $ cmake [path to source tree] -L	- This will show the basic options.
        $ cmake [path to source tree] -LA	- This will show the basic and advanced options.


## How do I specify how to build a shared or static library? {#sharedstatic}

    This is controlled with the internal `cmake` option, `BUILD_SHARED_LIBS`.

        $ cmake [Source Directory] -DBUILD_SHARED_LIBS=[ON/OFF]


## Can I build both shared and static libraries at the same time with cmake? {#sharedstaticboth}

    Not at this time; it is required to instead build first one version, and then the other, if you need both.

    ## What if I want to link against multiple libraries in a non-standard location? {#nonstdloc}

        You can specify the path to search when looking for dependencies and header files using the `CMAKE_PREFIX_PATH` variable:

    * Windows:

            $ cmake [Source Directory] -DCMAKE_PREFIX_PATH=c:\shared\libs\


    * Linux/Unix/OSX:

            $ cmake [Source Directory] -DCMAKE_PREFIX_PATH=/usr/custom_library_locations/		


## How can I specify linking against a particular library? {#partlib}

    It depends on the library.  To specify a custom `ZLib`, for example, you would do the following:

        $ cmake [Source Directory] -DZLIB_LIBRARY=/path/to/my/zlib.lib


    `HDF5` is more complex, since it requires both the `hdf5` and `hdf5_hl` libraries. You would specify custom `HDF5` libraries as follows:

        $ cmake [Source Directory] -DHDF5_LIB=/path/to/hdf5.lib \
            -DHDF5_HL_LIB=/path/to/hdf5_hl.lib \
            -DHDF5_INCLUDE_DIR=/path/to/hdf5/include


    Alternatively, you may specify:

        $ cmake [Source Directory] \
            -DHDF5_LIBRARIES="/path/to/hdf5.lib;/path/to/hdf5_hl.lib" \
            -DHDF5_INCLUDE_DIR=/path/to/hdf5/include/


## How can I specify a Parallel Build using HDF5 {#parallelhdf}

    If cmake is having problems finding the parallel `HDF5` install, you can specify the location manually:


        $ cmake [Source Directory] -DENABLE_PARALLEL=ON \
            -DHDF5_LIB=/usr/lib64/openmpi/lib/libhdf5.so \
            -DHDF5_HL_LIB=/usr/lib64/openmpi/lib/libhdf5.hl.so \
            -DHDF5_INCLUDE_DIR=/usr/include/openmpi-x86_64 \


    You will, of course, need to use the location of the libraries specific to your development environment.
