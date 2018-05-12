AC_DEFUN([CONFIGURE_XDR],
[
  AC_LINK_IFELSE([AC_LANG_PROGRAM([@%:@include <stdio.h>
                                   @%:@include <rpc/rpc.h>
                                   @%:@include <rpc/xdr.h>],
                                  [
                                    XDR * xdr;
                                    FILE * fp;
                                    xdrstdio_create(xdr, fp, XDR_ENCODE);
                                  ])],
                 [
                   AC_MSG_RESULT(yes)
                   AC_DEFINE(HAVE_XDR, 1, [Flag indicating headers and libraries for XDR IO are available])
                   enablexdr=yes
                 ],
                 [
                   AC_MSG_RESULT(no)
                   enablexdr=no
                 ])
])
