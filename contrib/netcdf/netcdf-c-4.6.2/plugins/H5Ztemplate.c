#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <hdf5.h>
/* Older versions of the hdf library may define H5PL_type_t here */
#include <H5PLextern.h>
#include "xxxx.h"

/*
Provide a textual template (not a C++ template)
from which one can construct a new filter.
The filter "name is marked with "XXXX" or "xxxx"

*/

const H5Z_class2_t H5Z_XXXX[1] = {{
    H5Z_CLASS_T_VERS,       /* H5Z_class_t version */
    (H5Z_filter_t)H5Z_FILTER_XXXX, /* Filter id number             */
    1,                             /* encoder_present flag (set to true) */
    1,                             /* decoder_present flag (set to true) */
    "xxxx",                        /* Filter name for debugging    */
    NULL,                          /* The "can apply" callback     */
    NULL,                          /* The "set local" callback     */
    (H5Z_func_t)H5Z_filter_xxxx,   /* The actual filter function   */
}};

/* External Discovery Functions */
H5PL_type_t
H5PLget_plugin_type(void)
{
    return H5PL_TYPE_FILTER;
}

const void*
H5PLget_plugin_info(void)
{
    return H5Z_XXXX;
}

size_t H5Z_filter_xxxx(unsigned int flags, size_t cd_nelmts,
		     const unsigned int cd_values[], size_t nbytes,
		     size_t *buf_size, void **buf);


size_t H5Z_filter_xxxx(unsigned int flags, size_t cd_nelmts,
                     const unsigned int cd_values[], size_t nbytes,
                     size_t *buf_size, void **buf)
{
  char *outbuf = NULL;
  size_t outbuflen, outdatalen;
  int ret;

  if (flags & H5Z_FLAG_REVERSE) {

    /** Decompress data.
     **
     ** This process is troublesome since the size of uncompressed data
     ** is unknown, so the low-level interface must be used.
     ** Data is decompressed to the output buffer (which is sized
     ** for the average case); if it gets full, its size is doubled
     ** and decompression continues.  This avoids repeatedly trying to
     ** decompress the whole block, which could be really inefficient.
     **/

    char *newbuf = NULL;
    size_t newbuflen;

  } else {

    /** Compress data.
     **
     ** This is quite simple, since the size of compressed data in the worst
     ** case is known and it is not much bigger than the size of uncompressed
     ** data.  This allows us to use the simplified one-shot interface to
     ** compression.
     **/

    unsigned int odatalen;  /* maybe not the same size as outdatalen */

    /* Prepare the output buffer. */
    outbuflen = M;  /* worst case */
    outbuf = H5allocate_memory(outbuflen,0);
    if (outbuf == NULL) {
      fprintf(stderr, "memory allocation failed for xxxx compression\n");
      goto cleanupAndFail;
    }

    /* Compress data. */

  }

  /* Always replace the input buffer with the output buffer. */
  H5free_memory(*buf);
  *buf = outbuf;
  *buf_size = outbuflen;
  return outdatalen;

 cleanupAndFail:
  if (outbuf)
    H5free_memory(outbuf);
  return 0;
}
