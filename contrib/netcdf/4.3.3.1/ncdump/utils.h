/*********************************************************************
 *   Copyright 2011, University Corporation for Atmospheric Research
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   
 *   Stuff that's common to both ncdump and nccopy
 *
 *********************************************************************/

#include "config.h"

#define	STREQ(a, b)	(*(a) == *(b) && strcmp((a), (b)) == 0)

/* Delimiter for separating netCDF groups in absolute pathnames, same as for HDF5 */
#define NC_GRP_DELIM '/'

typedef int bool_t;
enum {false=0, true=1};

struct safebuf_t;
/* Buffer structure for implementing growable strings, used in
 * preventing buffer overflows when the size needed for a character
 * buffer cannot be easily predicted or limited to any specific
 * maximum, such as when used in recursive function calls for nested
 * vlens and nested compound types. */
typedef struct safebuf_t {
    size_t len;			/* current length of buffer */
    size_t cl;			/* current length of string in buffer, < len-1 */
    char *buf;
} safebuf_t;

/* structure for list of ids, such as varids or grpids specified with -v or -g option */
typedef struct idnode {
    struct idnode* next;
    int id;
} idnode_t;

/* node in stack of group ids */
typedef struct grpnode {
    int grpid;
    struct grpnode *next;
} grpnode_t;

/* 
 * The opaque structure to hold state of iteration over groups.
 * (Just implemented as a stack of group ids.)
 */
typedef struct {
    int ngrps;			/* number of groups left to visit */
    grpnode_t *top;		/* group ids left to visit */
} ncgiter_t;

extern char *progname;		/* for error messages */

#ifndef NO_NETCDF_2
#define NO_NETCDF_2		/* assert we aren't using any netcdf-2 stuff */
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define NC_CHECK(fncall) {int statnc=fncall;if(statnc!=NC_NOERR)check(statnc,__FILE__,__LINE__);}

/* Print error message to stderr and exit */
extern void	error ( const char *fmt, ... );

/* Check error on malloc and exit with message if out of memory */
extern void*    emalloc ( size_t size );

/* Check error return.  If bad, print error message and exit. */
extern void check(int err, const char* file, const int line);

/* Return malloced name with chars special to CDL escaped. */
char* escaped_name(const char* cp);

/* Print name of netCDF var, dim, att, group, type, member, or enum
 * symbol with escaped special chars */
void print_name(const char *name);

/* Get dimid from a full dimension path name that may include group
 * names */
extern int  nc_inq_dimid2(int ncid, const char *dimname, int *dimidp);

/* Test if variable is a record variable */
extern int  isrecvar ( int ncid, int varid );

/* Get a new, empty id list. */
extern idnode_t* newidlist(void);

/* Add id to id list  */
extern void idadd(idnode_t* idlist, int id);

/* Test if id is in id list */
extern bool_t	idmember ( const idnode_t* idlist, int id );

/* Test if a group id is in group list */
extern bool_t	group_wanted ( int grpid, int nlgrps, const idnode_t* grpids );

/* Check group list for missing groups */
extern int grp_matches(int ncid, int nlgrps, char** lgrps, idnode_t *grpids);

/* Returns 1 if string s1 ends with string s2, 0 otherwise. */
extern int strendswith(const char *s1, const char *s2);

/* Within group with id ncid, get varid of variable with name varname
 * using nested group syntax "gp1/gp2/var" */
extern int nc_inq_gvarid ( int ncid, const char *varname, int *varidp );

/* Get variable id varid within group grpid using absolute or relative pathname for variable */
    extern int nc_inq_gvarid(int grpid, const char *varname, int *varidp);

/* Return how many variables are named varname in any groups in ncid */
extern size_t nc_inq_varname_count(int ncid, char *varname);

/* Check if any variable names specified in the list lvars (of length nlvars) are missing. */
extern int missing_vars(int ncid, int nlvars, char **lvars);

/* Make list of variables from comma-delimited string */
extern void make_lvars(char *optarg, int *nlvarsp, char ***lvarsp);

/* Make list of groups from comma-delimited string */
extern void make_lgrps(char *optarg, int *nlgrpsp, char*** lgrpsp, idnode_t **grpidsp);

/* Release an id list */
extern void freeidlist(idnode_t *idlist);

/* 
 * Simplest interface for group iteration: get total number of groups
 * (including all descendant groups, recursively) and all group ids
 * for start group and its descendants, in preorder.  Note that this
 * loses information about subgroup relationships, just flattening all
 * groups into a serial list.
 */
extern int nc_inq_grps_full(int ncid, int *numgrps, int *ncids);

/* 
 * More complex iterator interface: get group iterator for start group
 * ncid and all its descendant groups.
 */
extern int nc_get_giter(int ncid, ncgiter_t **iterp);

/* 
 * Get group id of next group.  On first call returns start group,
 * subsequently returns other subgroup ids in preorder.  Returns grpid
 * of 0 (never an actual group number) when no more groups.
 */
extern int nc_next_giter(ncgiter_t *iterp, int *grpid);

/*
 * Release memory allocated for group iterator.
 */
extern void nc_free_giter(ncgiter_t *iterp);

#ifdef __cplusplus
}
#endif
