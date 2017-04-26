# - adds support for the 'make dist' and 'make distcheck' commands	-*- cmake -*-
#
# Usage:
#   add_makedist()    ... called exactly once per project in the top-level
#                         CMakeLists.txt; it adds the 'dist' and 'distcheck'
#                         targets
#
#   enable_makedist(...list-of-sources...)
#                     ... called in every CMakeLists.txt which has source
#                         files. Beside the listed sources, the CMakeLists.txt
#                         file itself and all files listed in ${EXTRA_DIST}
#                         will be added to the tarball
#
# This module implements the 'make dist' and 'make distcheck'
# commands. When sources are given to enable_makedist() which are
# having the 'GENERATED' property, they will be ignored unless
# they are in ${EXTRA_DIST} or are having the EXTRADIST property.
#
# It supports the following variables:
#
#   MAKEDIST_TMPDIR   ... directory for temporary files
#   MAKEDIST_PKGBASE  ... basename of the created tarball; defaults to
#                         ${PACKAGE}-${VERSION} (see below)
#   MAKEDIST_TARFLAGS ... flags which are used to create the tarball
#   MAKEDIST_CMAKEFLAGS
#                     ... flags which are given to 'cmake' by 'make distcheck'
#   MAKEDIST_BUILDTARGETS
#                     ... the build-targets tried by 'make distcheck';
#                         defaults to nothing (--> all)
#   MAKEDIST_CHECKTARGETS
#                     ... the check-targets tried by 'make distcheck';
#                         defaults to 'test'
#   MAKEDIST_INSTALLTARGETS
#                     ... the install-targets tried by 'make distcheck';
#                         defaults to 'install'
#
#   EXTRA_DIST        ... set per CMakeLists.txt file to add non-source or
#                         generated source files to the tarball
#
# Unless MAKEDIST_PKGBASE is modified, the following variables are
# required to be set:
#
#   PACKAGE           ... name of the package (e.g. 'foo')
#   VERSION           ... version of the package (e.g. '0.1.2')
#
#
# Example:
#   --- top-level CMakeLists.txt ---
#   set(PACKAGE foo)
#   set(VERSION 1.2.3)
#
#   find_package(MakeDist)
#   add_subdirectory(foo)
#
#   add_makedist()
#
#   set(EXTRA_DIST COPYING)
#   enable_makedist(config.h.in)
#
#       --> adds ./COPYING and ./config.h.in to the tarball
#
#   --- foo/CMakeLists.txt ---
#   set_source_files_properties(
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-A.h
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-B.h
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-C.h
#     PROPERTIES GENERATED TRUE)
#
#   set_source_files_properties(
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-A.h
#     PROPERTIES EXTRADIST TRUE)
#   set(EXTRA_DIST ${CMAKE_CURRENT_BINARY_DIR}/generated-file-B.h)
#
#   enable_makedist(
#     foo.c
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-A.h
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-B.h
#     ${CMAKE_CURRENT_BINARY_DIR}/generated-file-C.h)
#
#       --> adds ./foo/generated-file-{A,B}.h and ./foo/foo.c to the
#           tarball, but not ${CMAKE_CURRENT_BINARY_DIR}/generated-file-C.h
#
# TODO/limitations:
#   * works probably with the 'make' generator only
#   * requires at least a SUSv3/POSIX compliant system; some parts
#     (e.g. bzip2-tarballs) depend on some GNUisms; this should be
#     made configurable perhaps
#   * there could be added some more tests (e.g. gnit-style (check for
#     existence of required files and a proper VERSION format))
#   * it is somehow ugly, to see hundreds of 'Built target .distdir'
#     messages
#   * is is ugly that all sources must be specified for enable_makedist();
#     there should be added support to 'cmake' to enumerate targets and
#     their sources


# Copyright (C) 2006 Enrico Scholz <enrico.scholz@informatik.tu-chemnitz.de>
#
# Redistribution and use, with or without modification, are permitted
# provided that the following conditions are met:
# 
#    1. Redistributions must retain the above copyright notice, this
#       list of conditions and the following disclaimer.
#    2. The name of the author may not be used to endorse or promote
#       products derived from this software without specific prior
#       written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


set(MakeDist_FOUND 1)

set(MAKEDIST_TMPDIR         "${CMAKE_BINARY_DIR}/.make-dist"                  CACHE PATH   "directory for temporary files created by'make dist*'")
set(MAKEDIST_PKGBASE        "${PACKAGE}-${VERSION}"                           CACHE STRING "basename of the tarball created by 'make dist*'")
set(MAKEDIST_TARFLAGS       --bzip2 			                      CACHE STRING "flags used by 'make dist' to create the tarball")
set(MAKEDIST_TARBALL        "${CMAKE_BINARY_DIR}/${MAKEDIST_PKGBASE}.tar.bz2" CACHE PATH   "tarball created by 'make dist'")
set(MAKEDIST_CMAKEFLAGS     -DBUILD_SHARED_LIBS:BOOL=ON                       CACHE STRING "flags which are given to 'cmake' by 'make distcheck'")
set(MAKEDIST_BUILDTARGETS   ""                                                CACHE STRING "build-target(s) tried by 'make distcheck'")
set(MAKEDIST_CHECKTARGETS   test                                              CACHE STRING "check-target(s) tried by 'make distcheck'")
set(MAKEDIST_INSTALLTARGETS install                                           CACHE STRING "install-target(s) tried by 'make distcheck'")

mark_as_advanced(
  MAKEDIST_TMPDIR MAKEDIST_PKGBASE MAKEDIST_TARBALL
  MAKEDIST_BUILDTARGETS MAKEDIST_CHECKTARGETS MAKEDIST_INSTALLTARGETS)


set(MAKEDIST_TARDIR   "${MAKEDIST_TMPDIR}/tar/${MAKEDIST_PKGBASE}")
set(MAKEDIST_CHECKDIR "${MAKEDIST_TMPDIR}/check")

macro(enable_makedist)
  add_custom_target(.distdir)
  set_source_files_properties(.distdir PROPERTIES
    SYMBOLIC TRUE
    GENERATED TRUE)

  foreach(__makedist_i ${EXTRA_DIST})
    set_source_files_properties("${__makedist_i}" PROPERTIES EXTRADIST TRUE)
  endforeach(__makedist_i)

  foreach(__makedist_source CMakeLists.txt ${EXTRA_DIST} ${ARGN})
    get_filename_component(__makedist_full_cur_path "${__makedist_source}" ABSOLUTE)
    get_source_file_property(__makedist_generated   "${__makedist_source}" GENERATED)
    get_source_file_property(__makedist_extradist   "${__makedist_source}" EXTRADIST)

    if(NOT __makedist_generated)
      file(RELATIVE_PATH __makedist_cur_path "${CMAKE_SOURCE_DIR}" "${__makedist_full_cur_path}")
    elseif(__makedist_extradist)
      file(RELATIVE_PATH __makedist_cur_path "${CMAKE_BINARY_DIR}" "${__makedist_full_cur_path}")
    else(__makedist_extradist)
      set(__makedist_cur_path)
    endif(NOT __makedist_generated)

    if(__makedist_cur_path)
      get_filename_component(__makedist_tmp "${MAKEDIST_TARDIR}/${__makedist_cur_path}" PATH)

      add_custom_command(TARGET .distdir
	COMMAND mkdir -p "${__makedist_tmp}"
	COMMAND cp -a "${__makedist_full_cur_path}" "${__makedist_tmp}/"
	DEPENDS ${__makedist_source})
    else(__makedist_cur_path)
      message(STATUS "skipping ${__makedist_source}")
    endif(__makedist_cur_path)
  endforeach(__makedist_source)
endmacro(enable_makedist)


macro(add_makedist)
  set_source_files_properties(dist distcheck .distcheck PROPERTIES
    SYMBOLIC TRUE
    GENERATED TRUE)

  add_custom_target(dist
    # cleanup tmpdir
    COMMAND chmod -R u+Xw "${MAKEDIST_TMPDIR}" 2>/dev/null || :
    COMMAND rm -rf "${MAKEDIST_TMPDIR}"

    # execute 'make .distdir'
    COMMAND ${CMAKE_MAKE_PROGRAM} .distdir
    # set proper permissions
    COMMAND chmod -R go-w,a+rX,u+w "${MAKEDIST_TARDIR}"

    # create the tarball
    COMMAND cd "${MAKEDIST_TARDIR}/.." && tar cf "${MAKEDIST_TARBALL}\${MAKEDIST_TARBALL_TMP}" ${MAKEDIST_TARFLAGS} ${MAKEDIST_PKGBASE}

    # cleanup tmpdir
    COMMAND rm -rf "${MAKEDIST_TMPDIR}"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )


  add_custom_target(distcheck
    COMMAND rm -f "${MAKEDIST_TARBALL}.tmp"
    COMMAND ${CMAKE_MAKE_PROGRAM} dist       MAKEDIST_TARBALL_TMP=.tmp
    COMMAND ${CMAKE_MAKE_PROGRAM} .distcheck MAKEDIST_TARBALL_TMP=.tmp
    COMMAND rm -f "${MAKEDIST_TARBALL}"
    COMMAND mv -f "${MAKEDIST_TARBALL}.tmp" "${MAKEDIST_TARBALL}")

    
  add_custom_target(.distcheck
    # cleanup tmpdir and create directory layout
    COMMAND chmod -R u+Xw "${MAKEDIST_CHECKDIR}" 2>/dev/null || :
    COMMAND rm -rf "${MAKEDIST_CHECKDIR}"
    COMMAND mkdir -p "${MAKEDIST_CHECKDIR}/source" "${MAKEDIST_CHECKDIR}/build"

    # extract tarball
    COMMAND cd "${MAKEDIST_CHECKDIR}/source" && tar xf "${MAKEDIST_TARBALL}\${MAKEDIST_TARBALL_TMP}" ${MAKEDIST_TARFLAGS}
    # write-protect sources to detect modifies-sourcetree bugs
    COMMAND chmod -R a-w "${MAKEDIST_CHECKDIR}/source"

    # invoke 'cmake'
    COMMAND ${CMAKE_COMMAND} -E echo "executing initial cmake"
    COMMAND cd "${MAKEDIST_CHECKDIR}/build"  &&
		    ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX:PATH="${MAKEDIST_CHECKDIR}/install"
		    ${MAKEDIST_CMAKEFLAGS} "${MAKEDIST_CHECKDIR}/source/${MAKEDIST_PKGBASE}"
    COMMAND ${CMAKE_COMMAND} -E echo "initial cmake succeeded"

    # execute 'make build' and 'make check'
    COMMAND ${CMAKE_COMMAND} -E echo "building project"
    COMMAND cd "${MAKEDIST_CHECKDIR}/build"  && ${CMAKE_MAKE_PROGRAM} ${MAKEDIST_BUILDTARGETS}
    COMMAND cd "${MAKEDIST_CHECKDIR}/build"  && ${CMAKE_MAKE_PROGRAM} ${MAKEDIST_CHECKTARGETS}

    # execute 'make install' without DESTDIR
    COMMAND cd "${MAKEDIST_CHECKDIR}/build"  && ${CMAKE_MAKE_PROGRAM} ${MAKEDIST_INSTALLTARGETS} DESTDIR=
    # write protect installation path to detect writing outside of DESTDIR
    COMMAND chmod -R a-w "${MAKEDIST_CHECKDIR}/install"
    # execute 'make install' with DESTDIR and move the files to a better location
    COMMAND cd "${MAKEDIST_CHECKDIR}/build"  && ${CMAKE_MAKE_PROGRAM} ${MAKEDIST_INSTALLTARGETS} DESTDIR="${MAKEDIST_CHECKDIR}/install-tmp"
    COMMAND mv "${MAKEDIST_CHECKDIR}/install-tmp/${MAKEDIST_CHECKDIR}/install" "${MAKEDIST_CHECKDIR}/install-destdir"

    # generate list of files which were installed by the both 'make
    # install' commands above and compare them
    COMMAND cd "${MAKEDIST_CHECKDIR}/install"         && find -type f | sort > ../files.install
    COMMAND cd "${MAKEDIST_CHECKDIR}/install-destdir" && find -type f | sort > ../files.destdir
    COMMAND cd "${MAKEDIST_CHECKDIR}" && diff files.install files.destdir

    # cleanup tmpdir
    COMMAND chmod -R u+Xw "${MAKEDIST_TMPDIR}" 2>/dev/null || :
    COMMAND rm -rf "${MAKEDIST_TMPDIR}"
    )
endmacro(add_makedist)
