#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>


static size_t var_chunksizes[4] = {4, 4, 4, 4} ;
static unsigned int var_filterparams[1] = {9U} ;

void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
        (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

int
main() {/* create bzip2.nc */

    int  stat;  /* return status */
    int  ncid;  /* netCDF id */

    /* group ids */
    int bzip2_grp;

    /* dimension ids */
    int dim0_dim;
    int dim1_dim;
    int dim2_dim;
    int dim3_dim;

    /* dimension lengths */
    size_t dim0_len = 4;
    size_t dim1_len = 4;
    size_t dim2_len = 4;
    size_t dim3_len = 4;

    /* variable ids */
    int var_id;

    /* rank (number of dimensions) for each variable */
#   define RANK_var 4

    /* variable shapes */
    int var_dims[RANK_var];

    /* enter define mode */
    stat = nc_create("bzip2.nc", NC_CLOBBER|NC_NETCDF4, &ncid);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_put_att_text(ncid, NC_GLOBAL, "_Format", 1, "netCDF-4");
    check_err(stat,__LINE__,__FILE__);
    bzip2_grp = ncid;

    /* define dimensions */
    stat = nc_def_dim(bzip2_grp, "dim0", dim0_len, &dim0_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(bzip2_grp, "dim1", dim1_len, &dim1_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(bzip2_grp, "dim2", dim2_len, &dim2_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(bzip2_grp, "dim3", dim3_len, &dim3_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    var_dims[0] = dim0_dim;
    var_dims[1] = dim1_dim;
    var_dims[2] = dim2_dim;
    var_dims[3] = dim3_dim;
    stat = nc_def_var(bzip2_grp, "var", NC_FLOAT, RANK_var, var_dims, &var_id);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_var_chunking(bzip2_grp, var_id, NC_CHUNKED, var_chunksizes);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_var_fill(bzip2_grp, var_id, NC_NOFILL, NULL);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_var_filter(bzip2_grp, var_id, 307, 1, var_filterparams);
    check_err(stat,__LINE__,__FILE__);

    /* leave define mode */
    stat = nc_enddef (bzip2_grp);
    check_err(stat,__LINE__,__FILE__);

    /* assign variable data */

    {
    float var_data[256] = {((float)0), ((float)1), ((float)2), ((float)3), ((float)4), ((float)5), ((float)6), ((float)7), ((float)8), ((float)9), ((float)10), ((float)11), ((float)12), ((float)13), ((float)14), ((float)15), ((float)16), ((float)17), ((float)18), ((float)19), ((float)20), ((float)21), ((float)22), ((float)23), ((float)24), ((float)25), ((float)26), ((float)27), ((float)28), ((float)29), ((float)30), ((float)31), ((float)32), ((float)33), ((float)34), ((float)35), ((float)36), ((float)37), ((float)38), ((float)39), ((float)40), ((float)41), ((float)42), ((float)43), ((float)44), ((float)45), ((float)46), ((float)47), ((float)48), ((float)49), ((float)50), ((float)51), ((float)52), ((float)53), ((float)54), ((float)55), ((float)56), ((float)57), ((float)58), ((float)59), ((float)60), ((float)61), ((float)62), ((float)63), ((float)64), ((float)65), ((float)66), ((float)67), ((float)68), ((float)69), ((float)70), ((float)71), ((float)72), ((float)73), ((float)74), ((float)75), ((float)76), ((float)77), ((float)78), ((float)79), ((float)80), ((float)81), ((float)82), ((float)83), ((float)84), ((float)85), ((float)86), ((float)87), ((float)88), ((float)89), ((float)90), ((float)91), ((float)92), ((float)93), ((float)94), ((float)95), ((float)96), ((float)97), ((float)98), ((float)99), ((float)100), ((float)101), ((float)102), ((float)103), ((float)104), ((float)105), ((float)106), ((float)107), ((float)108), ((float)109), ((float)110), ((float)111), ((float)112), ((float)113), ((float)114), ((float)115), ((float)116), ((float)117), ((float)118), ((float)119), ((float)120), ((float)121), ((float)122), ((float)123), ((float)124), ((float)125), ((float)126), ((float)127), ((float)128), ((float)129), ((float)130), ((float)131), ((float)132), ((float)133), ((float)134), ((float)135), ((float)136), ((float)137), ((float)138), ((float)139), ((float)140), ((float)141), ((float)142), ((float)143), ((float)144), ((float)145), ((float)146), ((float)147), ((float)148), ((float)149), ((float)150), ((float)151), ((float)152), ((float)153), ((float)154), ((float)155), ((float)156), ((float)157), ((float)158), ((float)159), ((float)160), ((float)161), ((float)162), ((float)163), ((float)164), ((float)165), ((float)166), ((float)167), ((float)168), ((float)169), ((float)170), ((float)171), ((float)172), ((float)173), ((float)174), ((float)175), ((float)176), ((float)177), ((float)178), ((float)179), ((float)180), ((float)181), ((float)182), ((float)183), ((float)184), ((float)185), ((float)186), ((float)187), ((float)188), ((float)189), ((float)190), ((float)191), ((float)192), ((float)193), ((float)194), ((float)195), ((float)196), ((float)197), ((float)198), ((float)199), ((float)200), ((float)201), ((float)202), ((float)203), ((float)204), ((float)205), ((float)206), ((float)207), ((float)208), ((float)209), ((float)210), ((float)211), ((float)212), ((float)213), ((float)214), ((float)215), ((float)216), ((float)217), ((float)218), ((float)219), ((float)220), ((float)221), ((float)222), ((float)223), ((float)224), ((float)225), ((float)226), ((float)227), ((float)228), ((float)229), ((float)230), ((float)231), ((float)232), ((float)233), ((float)234), ((float)235), ((float)236), ((float)237), ((float)238), ((float)239), ((float)240), ((float)241), ((float)242), ((float)243), ((float)244), ((float)245), ((float)246), ((float)247), ((float)248), ((float)249), ((float)250), ((float)251), ((float)252), ((float)253), ((float)254), ((float)255)} ;
    size_t var_startset[4] = {0, 0, 0, 0} ;
    size_t var_countset[4] = {4, 4, 4, 4};
    stat = nc_put_vara(bzip2_grp, var_id, var_startset, var_countset, var_data);
    check_err(stat,__LINE__,__FILE__);
    }


    stat = nc_close(bzip2_grp);
    check_err(stat,__LINE__,__FILE__);
    return 0;
}
