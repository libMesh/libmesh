/*********************************************************************
 *   Copyright 2011, University Corporation for Atmospheric Research
 *   See netcdf/README file for copying and redistribution conditions.
 *   $Id$
 *********************************************************************/

#include "config.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"

/*
 * Print error message to stderr and exit
 */
void
error(const char *fmt, ...)
{
    va_list args ;

    (void) fprintf(stderr,"%s: ", progname);
    va_start(args, fmt) ;
    (void) vfprintf(stderr,fmt,args) ;
    va_end(args) ;

    (void) fprintf(stderr, "\n") ;
    (void) fflush(stderr);	/* to ensure log files are current */
    exit(EXIT_FAILURE);
}

void *
emalloc (			/* check return from malloc */
	size_t size)
{
    void   *p;

    p = (void *) malloc (size==0 ? 1 : size); /* malloc(0) not portable */
    if (p == 0) {
	error ("out of memory\n");
    }
    return p;
}


void
check(int err, const char* file, const int line)
{
    fprintf(stderr,"%s\n",nc_strerror(err));
    fprintf(stderr,"Location: file %s; line %d\n", file,line);
    fflush(stderr); fflush(stdout);
    exit(1);
}


/* 
 * Returns malloced name with chars special to CDL escaped.
 * Caller should free result when done with it.
 */
char*
escaped_name(const char* cp) {
    char *ret;			/* string returned */
    char *sp;
    assert(cp != NULL);

    /* For some reason, and on some machines (e.g. tweety)
       utf8 characters such as \343 are considered control character. */
/*    if(*cp && (isspace(*cp) | iscntrl(*cp)))*/
    if((*cp >= 0x01 && *cp <= 0x20) || (*cp == 0x7f))
    {
	error("name begins with space or control-character: %c",*cp);
    }

    ret = emalloc(4*strlen(cp) + 1); /* max if every char escaped */
    sp = ret;
    *sp = 0;			    /* empty name OK */
    /* Special case: leading number allowed, but we must escape it for CDL */
    if((*cp >= '0' && *cp <= '9'))
    {
	*sp++ = '\\';
    }
    for (; *cp; cp++) {
	if (isascii((int)*cp)) {
	    if(iscntrl((int)*cp)) {	/* render control chars as two hex digits, \%xx */
		snprintf(sp, 4,"\\%%%.2x", *cp);
		sp += 4;
	    } else {
		switch (*cp) {
		case ' ':
		case '!':
		case '"':
		case '#':
		case '$':
		case '&':
		case '\'':
		case '(':
		case ')':
		case '*':
		case ',':
		case ':':
		case ';':
		case '<':
		case '=':
		case '>':
		case '?':
		case '[':
		case ']':
		case '\\':
		case '^':
		case '`':
		case '{':
		case '|':
		case '}':
		case '~':
		    *sp++ = '\\';
		    *sp++ = *cp;
		    break;
		default:		/* includes '/' */
		    *sp++ = *cp;
		    break;
		}
	    }
	} else { 		/* not ascii, assume just UTF-8 byte */
	    *sp++ = *cp;
	}
    }
    *sp = 0;
    return ret;
}


/* 
 * Print name with escapes for special characters
 */
void
print_name(const char* name) {
    char *ename = escaped_name(name);
    fputs(ename, stdout);
    free(ename);
}

/* Missing functionality that should be in nc_inq_dimid(), to get
 * dimid from a full dimension path name that may include group
 * names */
int 
nc_inq_dimid2(int ncid, const char *dimname, int *dimidp) {
    int ret = NC_NOERR;
    /* If '/' doesn't occur in dimname, just return id found by
     * nc_inq_dimid() */
    char *sp = strrchr(dimname, '/');
    if(!sp) { /* No '/' in dimname, so return nc_inq_dimid() result */
	ret = nc_inq_dimid(ncid, dimname, dimidp);
    } 
#ifdef USE_NETCDF4
    else {  /* Parse group name out and get dimid using that */
      size_t grp_namelen = sp - dimname;
      char *grpname = emalloc(grp_namelen+1);
      
      int grpid;
      strncpy(grpname, dimname, grp_namelen+1);
      grpname[grp_namelen] = '\0';
      ret = nc_inq_grp_full_ncid(ncid, grpname, &grpid);
      if(ret == NC_NOERR) {
	ret = nc_inq_dimid(grpid, dimname, dimidp);
      }
      free(grpname);
    }	
#endif	/* USE_NETCDF4 */
    return ret;
}


/*
 * return 1 if varid identifies a record variable
 * else return 0
 */
int
isrecvar(int ncid, int varid)
{
    int ndims;
    int is_recvar = 0;
    int *dimids;

    NC_CHECK( nc_inq_varndims(ncid, varid, &ndims) );
#ifdef USE_NETCDF4
    if (ndims > 0) {
	int nunlimdims;
	int *recdimids;
	int dim, recdim;
	dimids = (int *) emalloc((ndims + 1) * sizeof(int));
	NC_CHECK( nc_inq_vardimid(ncid, varid, dimids) );
	NC_CHECK( nc_inq_unlimdims(ncid, &nunlimdims, NULL) );
	recdimids = (int *) emalloc((nunlimdims + 1) * sizeof(int));
	NC_CHECK( nc_inq_unlimdims(ncid, NULL, recdimids) );
	for (dim = 0; dim < ndims && is_recvar == 0; dim++) {
	    for(recdim = 0; recdim < nunlimdims; recdim++) {
		if(dimids[dim] == recdimids[recdim]) {
		    is_recvar = 1;
		    break;
		}		
	    }
	}
	free(dimids);
	free(recdimids);
    }
#else
    if (ndims > 0) {
	int recdimid;
	dimids = (int *) emalloc((ndims + 1) * sizeof(int));
	NC_CHECK( nc_inq_vardimid(ncid, varid, dimids) );
	NC_CHECK( nc_inq_unlimdim(ncid, &recdimid) );
	if(dimids[0] == recdimid)
	    is_recvar = 1;
	free(dimids);
    }
#endif /* USE_NETCDF4 */
    return is_recvar;
}

static idnode_t*
newidnode(void) {
    idnode_t *newvp = (idnode_t*) emalloc(sizeof(idnode_t));
    return newvp;
}

/*
 * Get a new, empty variable list.
 */
idnode_t*
newidlist(void) {
    idnode_t *vp = newidnode();

    vp -> next = 0;
    vp -> id = -1;		/* bad id */

    return vp;
}

void
idadd(idnode_t* vlist, int varid) {
    idnode_t *newvp = newidnode();
    
    newvp -> next = vlist -> next;
    newvp -> id = varid;
    vlist -> next = newvp;
}

/* 
 * return true if id is member of list that idlist points to.
 */
bool_t
idmember(const idnode_t* idlist, int id)
{
    idnode_t *vp = idlist -> next;

    for (; vp ; vp = vp->next)
      if (vp->id == id)
	return true;
    return false;    
}

/*
 * Release a variable list.
 */
void
freeidlist(idnode_t *idlist)
{
   while(idlist) {
      idnode_t *vp = idlist->next;
      free(idlist);
      idlist = vp;
   }
}

/* 
 * Return true if group identified by grpid is member of grpids, a list of groups.
 * nlgrps is number of groups in the list.
 */
bool_t
group_wanted(int grpid, int nlgrps, const idnode_t* grpids)
{
    /* If -g not specified, all groups are wanted */
    if(nlgrps == 0) return true;
    /* if -g specified, look for match in group id list */
    return idmember(grpids, grpid);
}

/* Determine whether a group named formatting_specs.lgrps[igrp] exists
 * in a netCDF file or group with id ncid.  If so, return the count of
 * how many matching groups were found, else return a count of 0.  If
 * the name begins with "/", it is interpreted as an absolute group
 * name, in which case only 0 or 1 is returned.  Otherwise, interpret
 * it as a relative name, and the total number of occurrences within
 * the file/group identified by ncid is returned.  
 *
 * Also has side effect of updating the ngrpids and the associate
 * grpids array that represent the group list specified by the -g
 * option.  TODO: put this in its own function instead.
 */
static size_t
nc_inq_grpname_count(int ncid, int igrp, char **lgrps, idnode_t *grpids) {
    size_t count = 0;
#ifdef USE_NETCDF4
    int numgrps;
    int *ncids;
    int g;
    int grpid;
    int status;
#endif
    char *grpname = lgrps[igrp];

    /* permit empty string to also designate root group */
    if(grpname[0] == '\0' || STREQ(grpname,"/")) { 
	count = 1;
	idadd(grpids, ncid);
	return count;
    }
#ifdef USE_NETCDF4
    /* Handle absolute group names */
    if(grpname[0] == '/') {
	int grpid;
	status = nc_inq_grp_full_ncid(ncid, grpname, &grpid);
	if(status == NC_NOERR) {
	    count = 1;
	    idadd(grpids, grpid);
	} else if(status == NC_ENOGRP) {
	    count = 0;
	} else {
	    error("when looking up group %s: %s ", grpname, nc_strerror(status));
	}
	return count;
    }
    
    /* look in this group */
    status = nc_inq_grp_ncid(ncid, grpname, &grpid);
    if (status == NC_NOERR) {
	count++;
	idadd(grpids, grpid);
    }
    /* if this group has subgroups, call recursively on each of them */
    NC_CHECK( nc_inq_grps(ncid, &numgrps, NULL) );
    if(numgrps > 0) {
	/* Allocate memory to hold the list of group ids. */
	ncids = emalloc(numgrps * sizeof(int));
	/* Get the list of group ids. */
	NC_CHECK( nc_inq_grps(ncid, NULL, ncids) );
	/* Call this function recursively for each group. */
	for (g = 0; g < numgrps; g++) {
	    count += nc_inq_grpname_count(ncids[g], igrp, lgrps, grpids);
	}
	free(ncids);
    }
#endif /* USE_NETCDF4 */
    return count;    
}

/* Check if any group names specified with "-g grp1,...,grpn" are
 * missing.  Returns total number of matching groups if no missing
 * groups detected, otherwise exits. */
int
grp_matches(int ncid, int nlgrps, char** lgrps, idnode_t *grpids) {
    int ig;
    size_t total = 0;

    for (ig=0; ig < nlgrps; ig++) {
	size_t count = nc_inq_grpname_count(ncid, ig, lgrps, grpids);
	if(count == 0) {
	    error("%s: No such group", lgrps[ig]);
	    return 0;
	}
	total += count;
    }
    return total;
}

/* Returns 1 if string s1 ends with string s2, 0 otherwise. */
int
strendswith(const char *s1, const char *s2) {
    size_t m1 = strlen(s1);
    size_t m2 = strlen(s2);
    if (m1 < m2)
	return 0;
    return (strcmp(s1 + (m1 - m2), s2) == 0);
}

/* Get varid of variable with name using nested group syntax
 * "gp1/gp2/var" or "/gp1/gp2/var".  In the former case, grpname of
 * grp corresponding to grpid must end in "gp1/gp2".  In the latter
 * case, grpname for grpid must be exactly "/gp1/gp2".  If variable
 * named "var" is not in group grpid, returns NC_ENOTVAR, else sets
 * varid and returns NC_NOERR.  */
int 
nc_inq_gvarid(int grpid, const char *varname, int *varidp) {
    /* if varname has no "/" chars, then
          return varidp from nc_inq_varid(grpid, varname, varidp)
       if varname begins with "/"
          
       else
          get groupname corresponding to grpid
          get vargroup = substring of varname up to last "/"
          get relname = substring of varname after last "/"
          if (varname starts with "/" and groupname == vargroup) ||
             (groupname ends with vargroup)
             return nc_inq_varid(grpid, relname, varidp)
          else
             return NC_ENOTVAR
    */
    
#ifdef USE_NETCDF4
    char *vargroup;
    char *relname;
    char *groupname;
    int status;
    if (varname[0] == '\0')
	return NC_ENOTVAR;
    vargroup = strdup(varname);
    if (vargroup == NULL) 
	return NC_ENOMEM;
    relname = strrchr(vargroup, NC_GRP_DELIM);
    if (relname != NULL) {	/* name has a "/" in it */
	size_t len;		/* length of full group name for grpid */
	*relname++ = '\0';	/* split vargroup string in two,
				 * vargroup and relname */
	if ( (status = nc_inq_grpname_full(grpid, &len, NULL)) != NC_NOERR ) {
	    free(vargroup);
	    return status;
	}
	groupname = (char *)emalloc(len + 1);
	if ( (status = nc_inq_grpname_full(grpid, &len, groupname)) == NC_NOERR ) {
	    if(varname[0] == NC_GRP_DELIM) {
		if( strcmp(groupname, vargroup) == 0)
		    status = nc_inq_varid(grpid, relname, varidp);
		else
		    status = NC_ENOTVAR;
	    } else {
		if(strendswith(groupname, vargroup))
		    status = nc_inq_varid(grpid, relname, varidp);
		else
		    status = NC_ENOTVAR;
	    }
	}
	free(vargroup);
	free(groupname);
	return status;
    }
    free(vargroup);
#endif	/* USE_NETCDF4 */
    return nc_inq_varid(grpid, varname, varidp);
}

/* Determine whether a variable named varname exists in any group in
   an open netCDF file with id ncid.  If so, return the count of how
   many matching variables were found, else return a count of 0.  The
   variable name can be absolute such as "/foo" or "/GRP1/GRP1A/foo",
   in which case there is only one group to look in, given by the path
   from the root group.  Alternatively, the variable name can be
   relative, such as "foo" or "GRPA/GRPB/foo", in which case every
   group is examined for a variable with that relative name.  */
size_t
nc_inq_varname_count(int ncid, char *varname) {
    /* 
       count = 0;
       status = nc_inq_gvarid(ncid, varname, varid);
       if (status == NC_NOERR)
          count++;
       for each subgroup gid {
          count += nc_inq_varname_count(gid, varname);
       }
       return count;
    */
    size_t count = 0;
    int varid;
    /* look in this group */
    int status = nc_inq_gvarid(ncid, varname, &varid);
#ifdef USE_NETCDF4
    int numgrps;
    int *ncids;
    int g;
#endif

    if (status == NC_NOERR)
	count++;

#ifdef USE_NETCDF4
    /* if this group has subgroups, call recursively on each of them */
    NC_CHECK( nc_inq_grps(ncid, &numgrps, NULL) );
	 
    /* Allocate memory to hold the list of group ids. */
    ncids = emalloc((numgrps + 1) * sizeof(int));
	
    /* Get the list of group ids. */
    NC_CHECK( nc_inq_grps(ncid, NULL, ncids) );
	
    /* Call this function for each group. */
    for (g = 0; g < numgrps; g++) {
	count += nc_inq_varname_count(ncids[g], varname);
    }
    free(ncids);
#endif /* USE_NETCDF4 */
    return count;    
   
}

/* Check if any variable names specified with "-v var1,...,varn" are
 * missing.  Returns 0 if no missing variables detected, otherwise
 * exits. */
int
missing_vars(int ncid, int nlvars, char **lvars) {
    int iv;
    for (iv=0; iv < nlvars; iv++) {
	if(nc_inq_varname_count(ncid, lvars[iv]) == 0) {
	    error("%s: No such variable", lvars[iv]);
	}
    }
    return 0;
}

void
make_lvars(char *optarg, int *nlvarsp, char ***lvarsp)
{
    char *cp = optarg;
    int nvars = 1;
    char ** cpp;

    /* compute number of variable names in comma-delimited list */
    *nlvarsp = 1;
    while (*cp++)
      if (*cp == ',')
 	nvars++;
    *nlvarsp = nvars;
    *lvarsp = (char **) emalloc(nvars * sizeof(char*));
    cpp = *lvarsp;
    /* copy variable names into list */
    for (cp = strtok(optarg, ","); cp != NULL; cp = strtok((char *) NULL, ",")) {
	*cpp = strdup(cp);
	cpp++;
    }
}

void
make_lgrps(char *optarg, int *nlgrps, char ***lgrpsp, idnode_t **grpidsp)
{
    char *cp = optarg;
    int ngrps = 1;
    char ** cpp;

    /* compute number of group names in comma-delimited list */
    while (*cp++)
      if (*cp == ',')
 	ngrps++;
    *nlgrps = ngrps;
    *lgrpsp = (char **) emalloc(ngrps * sizeof(char*));
    cpp = *lgrpsp;
    /* copy group names into list */
    for (cp = strtok(optarg, ","); cp != NULL; cp = strtok((char *) NULL, ",")) {
	*cpp = strdup(cp);
	cpp++;
    }
    /* make empty list of grpids, to be filled in after input file opened */
    *grpidsp = newidlist();
}

/* initialize and return a new empty stack of grpids */
static ncgiter_t *
gs_init() {
    ncgiter_t *s = emalloc(sizeof(ncgiter_t));
    s->ngrps = 0;
    s->top = NULL;
    return s;
}

/* free a stack and all its nodes */
static void
gs_free(ncgiter_t *s) {
    grpnode_t *n0, *n1;
    n0 = s->top;
    while (n0) {
	n1 = n0->next;
	free(n0);
	n0 = n1;
    }
    free(s);
}

/* test if a stack is empty */
static int
gs_empty(ncgiter_t *s)
{
    return s->ngrps == 0;
}

/* push a grpid on stack */
static void
gs_push(ncgiter_t *s, int grpid)
{
    grpnode_t *node = emalloc(sizeof(grpnode_t));
 
    node->grpid = grpid;
    node->next = gs_empty(s) ? NULL : s->top;
    s->top = node;
    s->ngrps++;
}

/* pop value off stack and return */
static int 
gs_pop(ncgiter_t *s)
{
    if (gs_empty(s)) {
	return -1;		/* underflow, stack is empty */
    } else {			/* pop a node */
	grpnode_t *top = s->top;
	int value = top->grpid;
	s->top = top->next;
	/* TODO: first call to free gets seg fault with libumem */
	free(top);
	s->ngrps--;
	return value;
    }
}

#ifdef UNUSED
/* Return top value on stack without popping stack.  Defined for
 * completeness but not used (here). */
static int 
gs_top(ncgiter_t *s)
{
    if (gs_empty(s)) {
	return -1;		/* underflow, stack is empty */
    } else {			/* get top value */
	grpnode_t *top = s->top;
	int value = top->grpid;
	return value;
    }
}
#endif

/* Like netCDF-4 function nc_inq_grps(), but can be called from
 * netCDF-3 only code as well.  Maybe this is what nc_inq_grps()
 * should do if built without netCDF-4 data model support. */
static int
nc_inq_grps2(int ncid, int *numgrps, int *grpids)
{
    int stat = NC_NOERR;

    /* just check if ncid is valid id of open netCDF file */
    NC_CHECK(nc_inq(ncid, NULL, NULL, NULL, NULL));

#ifdef USE_NETCDF4
    NC_CHECK(nc_inq_grps(ncid, numgrps, grpids));
#else
    *numgrps = 0;
#endif
    return stat;
}

/* Initialize group iterator for start group and all its descendant
 * groups. */
int
nc_get_giter(int grpid,	       /* start group id */
	    ncgiter_t **iterp  /* returned opaque iteration state */
    ) 
{
    int stat = NC_NOERR;

    stat = nc_inq(grpid, NULL, NULL, NULL, NULL); /* check if grpid is valid */
    if(stat != NC_EBADGRPID && stat != NC_EBADID) {
	*iterp = gs_init();
	gs_push(*iterp, grpid);
    }

    return stat;
}

/* 
 * Get group id of next group.  On first call gets start group id,
 * subsequently returns other subgroup ids in preorder.  Returns zero
 * when no more groups left.
 */
int
nc_next_giter(ncgiter_t *iterp, int *grpidp) {
    int stat = NC_NOERR;
    int numgrps;
    int *grpids;
    int i;

    if(gs_empty(iterp)) {
	*grpidp = 0;		/* not a group, signals iterator is done */
    } else {
	*grpidp = gs_pop(iterp);
	NC_CHECK(nc_inq_grps2(*grpidp, &numgrps, NULL));
	if(numgrps > 0) {
	    grpids = (int *)emalloc(sizeof(int) * numgrps);
	    NC_CHECK(nc_inq_grps2(*grpidp, &numgrps, grpids));
	    for(i = numgrps - 1; i >= 0; i--) { /* push ids on stack in reverse order */
		gs_push(iterp, grpids[i]);
	    }
	    free(grpids);
	}
    }
    return stat;
}

/*
 * Release group iter.
 */
void
nc_free_giter(ncgiter_t *iterp)
{
    gs_free(iterp);
}

/* 
 * Get total number of groups (including the top-level group and all
 * descendant groups, recursively) and all descendant subgroup ids
 * (including the input rootid of the start group) for a group and
 * all its descendants, in preorder.
 *
 * If grpids or numgrps is NULL, it will be ignored.  So typical use
 * is to call with grpids NULL to get numgrps, allocate enough space
 * for the group ids, then call again to get them.
 */
int
nc_inq_grps_full(int rootid, int *numgrps, int *grpids) 
{
    int stat = NC_NOERR;
    ncgiter_t *giter;		/* pointer to group iterator */
    int grpid;
    size_t count;

    NC_CHECK(nc_get_giter(rootid, &giter));
    
    count = 0;
    NC_CHECK(nc_next_giter(giter, &grpid));
    while(grpid != 0) {
	if(grpids)
	    grpids[count] = grpid;
	count++;
	NC_CHECK(nc_next_giter(giter, &grpid));
    }
    if(numgrps)
	*numgrps = count;
    nc_free_giter(giter);
    return stat;
}
