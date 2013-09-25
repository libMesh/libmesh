SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

SET(cdmr_SOURCES nccrcvt.c nccr.c nccrdispatch.c crmeta.c nccrproto.c nccrgetvarx.c crdebug.c crutil.c nccrdump.c ncStreamx.c curlwrap.c ast_byteio.c ast_curl.c ast_internal.c ast_runtime.c ast_util.c ast_debug.c cceparse.c cceconstraints.c)

add_library(cdmr OBJECT ${cdmr_SOURCES}) 






