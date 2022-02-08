/*
 * This code is in public domain
 */

#ifndef _WIN_MKSTEMP_H_
#define _WIN_MKSTEMP_H_

#ifdef _MSC_VER
#include <sys/stat.h>
#include <fcntl.h>
#include <io.h>

inline int mkstemp(char *tmpl)
{
    char *fn = _mktemp(tmpl);
    if (fn == NULL)
        return -1;

    return _open(fn, _O_RDWR | _O_CREAT | _O_EXCL, _S_IREAD | _S_IWRITE);
}
#else
#error "mkstemp() is not implemented"
#endif

#endif // _WIN_MKSTEMP_H_
