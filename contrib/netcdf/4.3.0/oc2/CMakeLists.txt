SET(CMAKE_INCLUDE_CURRENT_DIR ON)

INCLUDE_DIRECTORIES(".")

#SET(oc_SOURCES occurlfunctions.c dapparse.c daplex.c ocdata.c ocuri.c ochttp.c ocbytes.c oc.c ocdebug.c ocdump.c ocinternal.c occompile.c oclist.c oclog.c ocnode.c ocutil.c occlientparams.c ocrc.c ocread.c xxdr.c daptab.c)
SET(oc_SOURCES oc.c daplex.c dapparse.c daptab.c occlientparams.c occompile.c occurlfunctions.c ocdata.c ocdebug.c ocdump.c ocinternal.c ocnode.c ochttp.c ocrc.c ocread.c ocutil.c ocbytes.c oclist.c ocuri.c oclog.c xxdr.c)

add_library(oc2 OBJECT ${oc_SOURCES}) 



