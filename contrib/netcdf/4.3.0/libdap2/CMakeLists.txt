SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(dap2_SOURCES constraints3.c common34.c dapcvt.c dapalign.c dapodom.c daputil.c ncdaperr.c cdf3.c cache.c dapdump.c dapdebug.c dapattr3.c ncd3dispatch.c ncdap3.c ncdap3a.c getvara3.c dceconstraints.c dcetab.c dceparse.c dcelex.c)

add_library(dap2 OBJECT ${dap2_SOURCES}) 

###
# Options related to the man page generation.
###

###
# The C API man page.
###

set(man_MANS netcdf.3)

set(ARGS_MANPAGE -DAPI=C)

IF (USE_NETCDF4)
	SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DNETCDF4=TRUE)
ENDIF(USE_NETCDF4)

IF (BUILD_DAP)
	SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DDAP=TRUE)
ENDIF(BUILD_DAP)

IF (BUILD_PARALLEL)
	SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DPARALLEL_IO=TRUE)
ENDIF(BUILD_PARALLEL)

