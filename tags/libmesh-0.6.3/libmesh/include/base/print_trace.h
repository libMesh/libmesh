#ifndef __print_trace_h__
#define __print_trace_h__

#include <iostream>

#include "libmesh_config.h"

/*
 * Print a stack trace (for code compiled with gcc)
 */
void print_trace(std::ostream &out = std::cerr);

#endif
