#include <stdio.h>
#include <stdlib.h>
#include <netcdf.h>


void
check_err(const int stat, const int line, const char *file) {
    if (stat != NC_NOERR) {
        (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
        fflush(stderr);
        exit(1);
    }
}

int
main() {/* create ctest0_64.nc */

    int  stat;  /* return status */
    int  ncid;  /* netCDF id */

    /* dimension ids */
    int Dr_dim;
    int D1_dim;
    int D2_dim;
    int D3_dim;
    int dim_MINUS_name_MINUS_dashes_dim;
    int dim_PERIOD_name_PERIOD_dots_dim;
    int dim_PLUS_name_PLUS_plusses_dim;
    int dim_ATSIGN_name_ATSIGN_ats_dim;

    /* dimension lengths */
    size_t Dr_len = NC_UNLIMITED;
    size_t D1_len = 1;
    size_t D2_len = 2;
    size_t D3_len = 3;
    size_t dim_MINUS_name_MINUS_dashes_len = 4;
    size_t dim_PERIOD_name_PERIOD_dots_len = 5;
    size_t dim_PLUS_name_PLUS_plusses_len = 6;
    size_t dim_ATSIGN_name_ATSIGN_ats_len = 7;

    /* variable ids */
    int c_id;
    int b_id;
    int s_id;
    int i_id;
    int f_id;
    int d_id;
    int cr_id;
    int br_id;
    int sr_id;
    int ir_id;
    int fr_id;
    int dr_id;
    int c1_id;
    int b1_id;
    int s1_id;
    int i1_id;
    int f1_id;
    int d1_id;
    int c2_id;
    int b2_id;
    int s2_id;
    int i2_id;
    int f2_id;
    int d2_id;
    int c3_id;
    int b3_id;
    int s3_id;
    int i3_id;
    int f3_id;
    int d3_id;
    int cr1_id;
    int br2_id;
    int sr3_id;
    int f11_id;
    int d12_id;
    int c13_id;
    int s21_id;
    int i22_id;
    int f23_id;
    int c31_id;
    int b32_id;
    int s33_id;
    int sr11_id;
    int ir12_id;
    int fr13_id;
    int cr21_id;
    int br22_id;
    int sr23_id;
    int fr31_id;
    int dr32_id;
    int cr33_id;
    int c111_id;
    int b112_id;
    int s113_id;
    int f121_id;
    int d122_id;
    int c123_id;
    int s131_id;
    int i132_id;
    int f133_id;
    int f211_id;
    int d212_id;
    int c213_id;
    int s221_id;
    int i222_id;
    int f223_id;
    int c231_id;
    int b232_id;
    int s233_id;
    int s311_id;
    int i312_id;
    int f313_id;
    int var_MINUS_name_MINUS_dashes_id;
    int var_PERIOD_name_PERIOD_dots_id;
    int var_PLUS_name_PLUS_plusses_id;
    int var_ATSIGN_name_ATSIGN_ats_id;

    /* rank (number of dimensions) for each variable */
#   define RANK_c 0
#   define RANK_b 0
#   define RANK_s 0
#   define RANK_i 0
#   define RANK_f 0
#   define RANK_d 0
#   define RANK_cr 1
#   define RANK_br 1
#   define RANK_sr 1
#   define RANK_ir 1
#   define RANK_fr 1
#   define RANK_dr 1
#   define RANK_c1 1
#   define RANK_b1 1
#   define RANK_s1 1
#   define RANK_i1 1
#   define RANK_f1 1
#   define RANK_d1 1
#   define RANK_c2 1
#   define RANK_b2 1
#   define RANK_s2 1
#   define RANK_i2 1
#   define RANK_f2 1
#   define RANK_d2 1
#   define RANK_c3 1
#   define RANK_b3 1
#   define RANK_s3 1
#   define RANK_i3 1
#   define RANK_f3 1
#   define RANK_d3 1
#   define RANK_cr1 2
#   define RANK_br2 2
#   define RANK_sr3 2
#   define RANK_f11 2
#   define RANK_d12 2
#   define RANK_c13 2
#   define RANK_s21 2
#   define RANK_i22 2
#   define RANK_f23 2
#   define RANK_c31 2
#   define RANK_b32 2
#   define RANK_s33 2
#   define RANK_sr11 3
#   define RANK_ir12 3
#   define RANK_fr13 3
#   define RANK_cr21 3
#   define RANK_br22 3
#   define RANK_sr23 3
#   define RANK_fr31 3
#   define RANK_dr32 3
#   define RANK_cr33 3
#   define RANK_c111 3
#   define RANK_b112 3
#   define RANK_s113 3
#   define RANK_f121 3
#   define RANK_d122 3
#   define RANK_c123 3
#   define RANK_s131 3
#   define RANK_i132 3
#   define RANK_f133 3
#   define RANK_f211 3
#   define RANK_d212 3
#   define RANK_c213 3
#   define RANK_s221 3
#   define RANK_i222 3
#   define RANK_f223 3
#   define RANK_c231 3
#   define RANK_b232 3
#   define RANK_s233 3
#   define RANK_s311 3
#   define RANK_i312 3
#   define RANK_f313 3
#   define RANK_var_MINUS_name_MINUS_dashes 0
#   define RANK_var_PERIOD_name_PERIOD_dots 0
#   define RANK_var_PLUS_name_PLUS_plusses 0
#   define RANK_var_ATSIGN_name_ATSIGN_ats 0

    /* variable shapes */
    int cr_dims[RANK_cr];
    int br_dims[RANK_br];
    int sr_dims[RANK_sr];
    int ir_dims[RANK_ir];
    int fr_dims[RANK_fr];
    int dr_dims[RANK_dr];
    int c1_dims[RANK_c1];
    int b1_dims[RANK_b1];
    int s1_dims[RANK_s1];
    int i1_dims[RANK_i1];
    int f1_dims[RANK_f1];
    int d1_dims[RANK_d1];
    int c2_dims[RANK_c2];
    int b2_dims[RANK_b2];
    int s2_dims[RANK_s2];
    int i2_dims[RANK_i2];
    int f2_dims[RANK_f2];
    int d2_dims[RANK_d2];
    int c3_dims[RANK_c3];
    int b3_dims[RANK_b3];
    int s3_dims[RANK_s3];
    int i3_dims[RANK_i3];
    int f3_dims[RANK_f3];
    int d3_dims[RANK_d3];
    int cr1_dims[RANK_cr1];
    int br2_dims[RANK_br2];
    int sr3_dims[RANK_sr3];
    int f11_dims[RANK_f11];
    int d12_dims[RANK_d12];
    int c13_dims[RANK_c13];
    int s21_dims[RANK_s21];
    int i22_dims[RANK_i22];
    int f23_dims[RANK_f23];
    int c31_dims[RANK_c31];
    int b32_dims[RANK_b32];
    int s33_dims[RANK_s33];
    int sr11_dims[RANK_sr11];
    int ir12_dims[RANK_ir12];
    int fr13_dims[RANK_fr13];
    int cr21_dims[RANK_cr21];
    int br22_dims[RANK_br22];
    int sr23_dims[RANK_sr23];
    int fr31_dims[RANK_fr31];
    int dr32_dims[RANK_dr32];
    int cr33_dims[RANK_cr33];
    int c111_dims[RANK_c111];
    int b112_dims[RANK_b112];
    int s113_dims[RANK_s113];
    int f121_dims[RANK_f121];
    int d122_dims[RANK_d122];
    int c123_dims[RANK_c123];
    int s131_dims[RANK_s131];
    int i132_dims[RANK_i132];
    int f133_dims[RANK_f133];
    int f211_dims[RANK_f211];
    int d212_dims[RANK_d212];
    int c213_dims[RANK_c213];
    int s221_dims[RANK_s221];
    int i222_dims[RANK_i222];
    int f223_dims[RANK_f223];
    int c231_dims[RANK_c231];
    int b232_dims[RANK_b232];
    int s233_dims[RANK_s233];
    int s311_dims[RANK_s311];
    int i312_dims[RANK_i312];
    int f313_dims[RANK_f313];

    /* enter define mode */
    stat = nc_create("ctest0_64.nc", NC_CLOBBER|NC_64BIT_OFFSET, &ncid);
    check_err(stat,__LINE__,__FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "Dr", Dr_len, &Dr_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "D1", D1_len, &D1_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "D2", D2_len, &D2_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "D3", D3_len, &D3_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "dim-name-dashes", dim_MINUS_name_MINUS_dashes_len, &dim_MINUS_name_MINUS_dashes_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "dim.name.dots", dim_PERIOD_name_PERIOD_dots_len, &dim_PERIOD_name_PERIOD_dots_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "dim+name+plusses", dim_PLUS_name_PLUS_plusses_len, &dim_PLUS_name_PLUS_plusses_dim);
    check_err(stat,__LINE__,__FILE__);
    stat = nc_def_dim(ncid, "dim@name@ats", dim_ATSIGN_name_ATSIGN_ats_len, &dim_ATSIGN_name_ATSIGN_ats_dim);
    check_err(stat,__LINE__,__FILE__);

    /* define variables */

    stat = nc_def_var(ncid, "c", NC_CHAR, RANK_c, 0, &c_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "b", NC_BYTE, RANK_b, 0, &b_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "s", NC_SHORT, RANK_s, 0, &s_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "i", NC_INT, RANK_i, 0, &i_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "f", NC_FLOAT, RANK_f, 0, &f_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "d", NC_DOUBLE, RANK_d, 0, &d_id);
    check_err(stat,__LINE__,__FILE__);

    cr_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "cr", NC_CHAR, RANK_cr, cr_dims, &cr_id);
    check_err(stat,__LINE__,__FILE__);

    br_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "br", NC_BYTE, RANK_br, br_dims, &br_id);
    check_err(stat,__LINE__,__FILE__);

    sr_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "sr", NC_SHORT, RANK_sr, sr_dims, &sr_id);
    check_err(stat,__LINE__,__FILE__);

    ir_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "ir", NC_INT, RANK_ir, ir_dims, &ir_id);
    check_err(stat,__LINE__,__FILE__);

    fr_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "fr", NC_FLOAT, RANK_fr, fr_dims, &fr_id);
    check_err(stat,__LINE__,__FILE__);

    dr_dims[0] = Dr_dim;
    stat = nc_def_var(ncid, "dr", NC_DOUBLE, RANK_dr, dr_dims, &dr_id);
    check_err(stat,__LINE__,__FILE__);

    c1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "c1", NC_CHAR, RANK_c1, c1_dims, &c1_id);
    check_err(stat,__LINE__,__FILE__);

    b1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "b1", NC_BYTE, RANK_b1, b1_dims, &b1_id);
    check_err(stat,__LINE__,__FILE__);

    s1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "s1", NC_SHORT, RANK_s1, s1_dims, &s1_id);
    check_err(stat,__LINE__,__FILE__);

    i1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "i1", NC_INT, RANK_i1, i1_dims, &i1_id);
    check_err(stat,__LINE__,__FILE__);

    f1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "f1", NC_FLOAT, RANK_f1, f1_dims, &f1_id);
    check_err(stat,__LINE__,__FILE__);

    d1_dims[0] = D1_dim;
    stat = nc_def_var(ncid, "d1", NC_DOUBLE, RANK_d1, d1_dims, &d1_id);
    check_err(stat,__LINE__,__FILE__);

    c2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "c2", NC_CHAR, RANK_c2, c2_dims, &c2_id);
    check_err(stat,__LINE__,__FILE__);

    b2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "b2", NC_BYTE, RANK_b2, b2_dims, &b2_id);
    check_err(stat,__LINE__,__FILE__);

    s2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "s2", NC_SHORT, RANK_s2, s2_dims, &s2_id);
    check_err(stat,__LINE__,__FILE__);

    i2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "i2", NC_INT, RANK_i2, i2_dims, &i2_id);
    check_err(stat,__LINE__,__FILE__);

    f2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "f2", NC_FLOAT, RANK_f2, f2_dims, &f2_id);
    check_err(stat,__LINE__,__FILE__);

    d2_dims[0] = D2_dim;
    stat = nc_def_var(ncid, "d2", NC_DOUBLE, RANK_d2, d2_dims, &d2_id);
    check_err(stat,__LINE__,__FILE__);

    c3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "c3", NC_CHAR, RANK_c3, c3_dims, &c3_id);
    check_err(stat,__LINE__,__FILE__);

    b3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "b3", NC_BYTE, RANK_b3, b3_dims, &b3_id);
    check_err(stat,__LINE__,__FILE__);

    s3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "s3", NC_SHORT, RANK_s3, s3_dims, &s3_id);
    check_err(stat,__LINE__,__FILE__);

    i3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "i3", NC_INT, RANK_i3, i3_dims, &i3_id);
    check_err(stat,__LINE__,__FILE__);

    f3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "f3", NC_FLOAT, RANK_f3, f3_dims, &f3_id);
    check_err(stat,__LINE__,__FILE__);

    d3_dims[0] = D3_dim;
    stat = nc_def_var(ncid, "d3", NC_DOUBLE, RANK_d3, d3_dims, &d3_id);
    check_err(stat,__LINE__,__FILE__);

    cr1_dims[0] = Dr_dim;
    cr1_dims[1] = D1_dim;
    stat = nc_def_var(ncid, "cr1", NC_CHAR, RANK_cr1, cr1_dims, &cr1_id);
    check_err(stat,__LINE__,__FILE__);

    br2_dims[0] = Dr_dim;
    br2_dims[1] = D2_dim;
    stat = nc_def_var(ncid, "br2", NC_BYTE, RANK_br2, br2_dims, &br2_id);
    check_err(stat,__LINE__,__FILE__);

    sr3_dims[0] = Dr_dim;
    sr3_dims[1] = D3_dim;
    stat = nc_def_var(ncid, "sr3", NC_SHORT, RANK_sr3, sr3_dims, &sr3_id);
    check_err(stat,__LINE__,__FILE__);

    f11_dims[0] = D1_dim;
    f11_dims[1] = D1_dim;
    stat = nc_def_var(ncid, "f11", NC_FLOAT, RANK_f11, f11_dims, &f11_id);
    check_err(stat,__LINE__,__FILE__);

    d12_dims[0] = D1_dim;
    d12_dims[1] = D2_dim;
    stat = nc_def_var(ncid, "d12", NC_DOUBLE, RANK_d12, d12_dims, &d12_id);
    check_err(stat,__LINE__,__FILE__);

    c13_dims[0] = D1_dim;
    c13_dims[1] = D3_dim;
    stat = nc_def_var(ncid, "c13", NC_CHAR, RANK_c13, c13_dims, &c13_id);
    check_err(stat,__LINE__,__FILE__);

    s21_dims[0] = D2_dim;
    s21_dims[1] = D1_dim;
    stat = nc_def_var(ncid, "s21", NC_SHORT, RANK_s21, s21_dims, &s21_id);
    check_err(stat,__LINE__,__FILE__);

    i22_dims[0] = D2_dim;
    i22_dims[1] = D2_dim;
    stat = nc_def_var(ncid, "i22", NC_INT, RANK_i22, i22_dims, &i22_id);
    check_err(stat,__LINE__,__FILE__);

    f23_dims[0] = D2_dim;
    f23_dims[1] = D3_dim;
    stat = nc_def_var(ncid, "f23", NC_FLOAT, RANK_f23, f23_dims, &f23_id);
    check_err(stat,__LINE__,__FILE__);

    c31_dims[0] = D3_dim;
    c31_dims[1] = D1_dim;
    stat = nc_def_var(ncid, "c31", NC_CHAR, RANK_c31, c31_dims, &c31_id);
    check_err(stat,__LINE__,__FILE__);

    b32_dims[0] = D3_dim;
    b32_dims[1] = D2_dim;
    stat = nc_def_var(ncid, "b32", NC_BYTE, RANK_b32, b32_dims, &b32_id);
    check_err(stat,__LINE__,__FILE__);

    s33_dims[0] = D3_dim;
    s33_dims[1] = D3_dim;
    stat = nc_def_var(ncid, "s33", NC_SHORT, RANK_s33, s33_dims, &s33_id);
    check_err(stat,__LINE__,__FILE__);

    sr11_dims[0] = Dr_dim;
    sr11_dims[1] = D1_dim;
    sr11_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "sr11", NC_SHORT, RANK_sr11, sr11_dims, &sr11_id);
    check_err(stat,__LINE__,__FILE__);

    ir12_dims[0] = Dr_dim;
    ir12_dims[1] = D1_dim;
    ir12_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "ir12", NC_INT, RANK_ir12, ir12_dims, &ir12_id);
    check_err(stat,__LINE__,__FILE__);

    fr13_dims[0] = Dr_dim;
    fr13_dims[1] = D1_dim;
    fr13_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "fr13", NC_FLOAT, RANK_fr13, fr13_dims, &fr13_id);
    check_err(stat,__LINE__,__FILE__);

    cr21_dims[0] = Dr_dim;
    cr21_dims[1] = D2_dim;
    cr21_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "cr21", NC_CHAR, RANK_cr21, cr21_dims, &cr21_id);
    check_err(stat,__LINE__,__FILE__);

    br22_dims[0] = Dr_dim;
    br22_dims[1] = D2_dim;
    br22_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "br22", NC_BYTE, RANK_br22, br22_dims, &br22_id);
    check_err(stat,__LINE__,__FILE__);

    sr23_dims[0] = Dr_dim;
    sr23_dims[1] = D2_dim;
    sr23_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "sr23", NC_SHORT, RANK_sr23, sr23_dims, &sr23_id);
    check_err(stat,__LINE__,__FILE__);

    fr31_dims[0] = Dr_dim;
    fr31_dims[1] = D3_dim;
    fr31_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "fr31", NC_FLOAT, RANK_fr31, fr31_dims, &fr31_id);
    check_err(stat,__LINE__,__FILE__);

    dr32_dims[0] = Dr_dim;
    dr32_dims[1] = D3_dim;
    dr32_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "dr32", NC_DOUBLE, RANK_dr32, dr32_dims, &dr32_id);
    check_err(stat,__LINE__,__FILE__);

    cr33_dims[0] = Dr_dim;
    cr33_dims[1] = D3_dim;
    cr33_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "cr33", NC_CHAR, RANK_cr33, cr33_dims, &cr33_id);
    check_err(stat,__LINE__,__FILE__);

    c111_dims[0] = D1_dim;
    c111_dims[1] = D1_dim;
    c111_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "c111", NC_CHAR, RANK_c111, c111_dims, &c111_id);
    check_err(stat,__LINE__,__FILE__);

    b112_dims[0] = D1_dim;
    b112_dims[1] = D1_dim;
    b112_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "b112", NC_BYTE, RANK_b112, b112_dims, &b112_id);
    check_err(stat,__LINE__,__FILE__);

    s113_dims[0] = D1_dim;
    s113_dims[1] = D1_dim;
    s113_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "s113", NC_SHORT, RANK_s113, s113_dims, &s113_id);
    check_err(stat,__LINE__,__FILE__);

    f121_dims[0] = D1_dim;
    f121_dims[1] = D2_dim;
    f121_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "f121", NC_FLOAT, RANK_f121, f121_dims, &f121_id);
    check_err(stat,__LINE__,__FILE__);

    d122_dims[0] = D1_dim;
    d122_dims[1] = D2_dim;
    d122_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "d122", NC_DOUBLE, RANK_d122, d122_dims, &d122_id);
    check_err(stat,__LINE__,__FILE__);

    c123_dims[0] = D1_dim;
    c123_dims[1] = D2_dim;
    c123_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "c123", NC_CHAR, RANK_c123, c123_dims, &c123_id);
    check_err(stat,__LINE__,__FILE__);

    s131_dims[0] = D1_dim;
    s131_dims[1] = D3_dim;
    s131_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "s131", NC_SHORT, RANK_s131, s131_dims, &s131_id);
    check_err(stat,__LINE__,__FILE__);

    i132_dims[0] = D1_dim;
    i132_dims[1] = D3_dim;
    i132_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "i132", NC_INT, RANK_i132, i132_dims, &i132_id);
    check_err(stat,__LINE__,__FILE__);

    f133_dims[0] = D1_dim;
    f133_dims[1] = D3_dim;
    f133_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "f133", NC_FLOAT, RANK_f133, f133_dims, &f133_id);
    check_err(stat,__LINE__,__FILE__);

    f211_dims[0] = D2_dim;
    f211_dims[1] = D1_dim;
    f211_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "f211", NC_FLOAT, RANK_f211, f211_dims, &f211_id);
    check_err(stat,__LINE__,__FILE__);

    d212_dims[0] = D2_dim;
    d212_dims[1] = D1_dim;
    d212_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "d212", NC_DOUBLE, RANK_d212, d212_dims, &d212_id);
    check_err(stat,__LINE__,__FILE__);

    c213_dims[0] = D2_dim;
    c213_dims[1] = D1_dim;
    c213_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "c213", NC_CHAR, RANK_c213, c213_dims, &c213_id);
    check_err(stat,__LINE__,__FILE__);

    s221_dims[0] = D2_dim;
    s221_dims[1] = D2_dim;
    s221_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "s221", NC_SHORT, RANK_s221, s221_dims, &s221_id);
    check_err(stat,__LINE__,__FILE__);

    i222_dims[0] = D2_dim;
    i222_dims[1] = D2_dim;
    i222_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "i222", NC_INT, RANK_i222, i222_dims, &i222_id);
    check_err(stat,__LINE__,__FILE__);

    f223_dims[0] = D2_dim;
    f223_dims[1] = D2_dim;
    f223_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "f223", NC_FLOAT, RANK_f223, f223_dims, &f223_id);
    check_err(stat,__LINE__,__FILE__);

    c231_dims[0] = D2_dim;
    c231_dims[1] = D3_dim;
    c231_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "c231", NC_CHAR, RANK_c231, c231_dims, &c231_id);
    check_err(stat,__LINE__,__FILE__);

    b232_dims[0] = D2_dim;
    b232_dims[1] = D3_dim;
    b232_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "b232", NC_BYTE, RANK_b232, b232_dims, &b232_id);
    check_err(stat,__LINE__,__FILE__);

    s233_dims[0] = D2_dim;
    s233_dims[1] = D3_dim;
    s233_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "s233", NC_SHORT, RANK_s233, s233_dims, &s233_id);
    check_err(stat,__LINE__,__FILE__);

    s311_dims[0] = D3_dim;
    s311_dims[1] = D1_dim;
    s311_dims[2] = D1_dim;
    stat = nc_def_var(ncid, "s311", NC_SHORT, RANK_s311, s311_dims, &s311_id);
    check_err(stat,__LINE__,__FILE__);

    i312_dims[0] = D3_dim;
    i312_dims[1] = D1_dim;
    i312_dims[2] = D2_dim;
    stat = nc_def_var(ncid, "i312", NC_INT, RANK_i312, i312_dims, &i312_id);
    check_err(stat,__LINE__,__FILE__);

    f313_dims[0] = D3_dim;
    f313_dims[1] = D1_dim;
    f313_dims[2] = D3_dim;
    stat = nc_def_var(ncid, "f313", NC_FLOAT, RANK_f313, f313_dims, &f313_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "var-name-dashes", NC_DOUBLE, RANK_var_MINUS_name_MINUS_dashes, 0, &var_MINUS_name_MINUS_dashes_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "var.name.dots", NC_DOUBLE, RANK_var_PERIOD_name_PERIOD_dots, 0, &var_PERIOD_name_PERIOD_dots_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "var+name+plusses", NC_DOUBLE, RANK_var_PLUS_name_PLUS_plusses, 0, &var_PLUS_name_PLUS_plusses_id);
    check_err(stat,__LINE__,__FILE__);

    stat = nc_def_var(ncid, "var@name@ats", NC_DOUBLE, RANK_var_ATSIGN_name_ATSIGN_ats, 0, &var_ATSIGN_name_ATSIGN_ats_id);
    check_err(stat,__LINE__,__FILE__);

    /* assign global attributes */

    {
    stat = nc_put_att_text(ncid, NC_GLOBAL, "Gc", 1, "");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const signed char c0_Gb_att[2] = {-128, 127} ;
    stat = nc_put_att_schar(ncid, NC_GLOBAL, "Gb", NC_BYTE, 2, c0_Gb_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const short c0_Gs_att[3] = {-32768, 0, 32767} ;
    stat = nc_put_att_short(ncid, NC_GLOBAL, "Gs", NC_SHORT, 3, c0_Gs_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_Gi_att[3] = {-2147483647, 0, 2147483647} ;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "Gi", NC_INT, 3, c0_Gi_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float c0_Gf_att[3] = {((float)-9.9999996e+35), ((float)0), ((float)9.9999996e+35)} ;
    stat = nc_put_att_float(ncid, NC_GLOBAL, "Gf", NC_FLOAT, 3, c0_Gf_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double c0_Gd_att[3] = {((double)-1e+308), ((double)0), ((double)1e+308)} ;
    stat = nc_put_att_double(ncid, NC_GLOBAL, "Gd", NC_DOUBLE, 3, c0_Gd_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_Gatt_MINUS_name_MINUS_dashes_att[1] = {-1} ;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "Gatt-name-dashes", NC_INT, 1, c0_Gatt_MINUS_name_MINUS_dashes_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_Gatt_DOT_name_DOT_dots_att[1] = {-2} ;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "Gatt.name.dots", NC_INT, 1, c0_Gatt_DOT_name_DOT_dots_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_Gatt_PLUS_name_PLUS_plusses_att[1] = {-3} ;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "Gatt+name+plusses", NC_INT, 1, c0_Gatt_PLUS_name_PLUS_plusses_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_Gatt_ATSIGN_name_ATSIGN_ats_att[1] = {-4} ;
    stat = nc_put_att_int(ncid, NC_GLOBAL, "Gatt@name@ats", NC_INT, 1, c0_Gatt_ATSIGN_name_ATSIGN_ats_att);
    check_err(stat,__LINE__,__FILE__);
    }


    /* assign per-variable attributes */

    {
    static const int c0_att_MINUS_name_MINUS_dashes_att[1] = {4} ;
    stat = nc_put_att_int(ncid, c_id, "att-name-dashes", NC_INT, 1, c0_att_MINUS_name_MINUS_dashes_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_att_DOT_name_DOT_dots_att[1] = {5} ;
    stat = nc_put_att_int(ncid, c_id, "att.name.dots", NC_INT, 1, c0_att_DOT_name_DOT_dots_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_att_PLUS_name_PLUS_plusses_att[1] = {6} ;
    stat = nc_put_att_int(ncid, c_id, "att+name+plusses", NC_INT, 1, c0_att_PLUS_name_PLUS_plusses_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_att_ATSIGN_name_ATSIGN_ats_att[1] = {7} ;
    stat = nc_put_att_int(ncid, c_id, "att@name@ats", NC_INT, 1, c0_att_ATSIGN_name_ATSIGN_ats_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, b_id, "c", 1, "");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const signed char c0_b_att[4] = {0, 127, -128, -1} ;
    stat = nc_put_att_schar(ncid, s_id, "b", NC_BYTE, 4, c0_b_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const short c0_s_att[3] = {-32768, 0, 32767} ;
    stat = nc_put_att_short(ncid, s_id, "s", NC_SHORT, 3, c0_s_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const int c0_i_att[3] = {-2147483647, 0, 2147483647} ;
    stat = nc_put_att_int(ncid, i_id, "i", NC_INT, 3, c0_i_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const float c0_f_att[3] = {((float)-9.9999996e+35), ((float)0), ((float)9.9999996e+35)} ;
    stat = nc_put_att_float(ncid, i_id, "f", NC_FLOAT, 3, c0_f_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    static const double c0_d_att[3] = {((double)-1e+308), ((double)0), ((double)1e+308)} ;
    stat = nc_put_att_double(ncid, i_id, "d", NC_DOUBLE, 3, c0_d_att);
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, f_id, "c", 1, "x");
    check_err(stat,__LINE__,__FILE__);
    }

    {
    stat = nc_put_att_text(ncid, d_id, "c", 8, "abcd\tZ$&");
    check_err(stat,__LINE__,__FILE__);
    }


    /* leave define mode */
    stat = nc_enddef (ncid);
    check_err(stat,__LINE__,__FILE__);

    /* assign variable data */

    {
    size_t count = 0;
    static char c_data[1] = {'2'};
    stat = nc_put_var1(ncid, c_id, &count, c_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static signed char b_data[1] = {-2};
    stat = nc_put_var1(ncid, b_id, &count, b_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static short s_data[1] = {-5};
    stat = nc_put_var1(ncid, s_id, &count, s_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static int i_data[1] = {-20};
    stat = nc_put_var1(ncid, i_id, &count, i_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static float f_data[1] = {((float)-9)};
    stat = nc_put_var1(ncid, f_id, &count, f_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static double d_data[1] = {((double)-10)};
    stat = nc_put_var1(ncid, d_id, &count, d_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* cr_data = "ab" ;
    size_t cr_startset[1] = {0} ;
    size_t cr_countset[1] = {2};
    stat = nc_put_vara(ncid, cr_id, cr_startset, cr_countset, cr_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char br_data[2] = {-128, 127} ;
    size_t br_startset[1] = {0} ;
    size_t br_countset[1] = {2};
    stat = nc_put_vara(ncid, br_id, br_startset, br_countset, br_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short sr_data[2] = {-32768, 32767} ;
    size_t sr_startset[1] = {0} ;
    size_t sr_countset[1] = {2};
    stat = nc_put_vara(ncid, sr_id, sr_startset, sr_countset, sr_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int ir_data[2] = {-2147483646, 2147483647} ;
    size_t ir_startset[1] = {0} ;
    size_t ir_countset[1] = {2};
    stat = nc_put_vara(ncid, ir_id, ir_startset, ir_countset, ir_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float fr_data[2] = {((float)-9.9999996e+35), ((float)9.9999996e+35)} ;
    size_t fr_startset[1] = {0} ;
    size_t fr_countset[1] = {2};
    stat = nc_put_vara(ncid, fr_id, fr_startset, fr_countset, fr_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double dr_data[2] = {((double)-1e+308), ((double)1e+308)} ;
    size_t dr_startset[1] = {0} ;
    size_t dr_countset[1] = {2};
    stat = nc_put_vara(ncid, dr_id, dr_startset, dr_countset, dr_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c1_data = "\000" ;
    size_t c1_startset[1] = {0} ;
    size_t c1_countset[1] = {1};
    stat = nc_put_vara(ncid, c1_id, c1_startset, c1_countset, c1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b1_data[1] = {-128} ;
    size_t b1_startset[1] = {0} ;
    size_t b1_countset[1] = {1};
    stat = nc_put_vara(ncid, b1_id, b1_startset, b1_countset, b1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s1_data[1] = {-32768} ;
    size_t s1_startset[1] = {0} ;
    size_t s1_countset[1] = {1};
    stat = nc_put_vara(ncid, s1_id, s1_startset, s1_countset, s1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i1_data[1] = {-2147483646} ;
    size_t i1_startset[1] = {0} ;
    size_t i1_countset[1] = {1};
    stat = nc_put_vara(ncid, i1_id, i1_startset, i1_countset, i1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f1_data[1] = {((float)-9.9999996e+35)} ;
    size_t f1_startset[1] = {0} ;
    size_t f1_countset[1] = {1};
    stat = nc_put_vara(ncid, f1_id, f1_startset, f1_countset, f1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d1_data[1] = {((double)-1e+308)} ;
    size_t d1_startset[1] = {0} ;
    size_t d1_countset[1] = {1};
    stat = nc_put_vara(ncid, d1_id, d1_startset, d1_countset, d1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c2_data = "ab" ;
    size_t c2_startset[1] = {0} ;
    size_t c2_countset[1] = {2};
    stat = nc_put_vara(ncid, c2_id, c2_startset, c2_countset, c2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b2_data[2] = {-128, 127} ;
    size_t b2_startset[1] = {0} ;
    size_t b2_countset[1] = {2};
    stat = nc_put_vara(ncid, b2_id, b2_startset, b2_countset, b2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s2_data[2] = {-32768, 32767} ;
    size_t s2_startset[1] = {0} ;
    size_t s2_countset[1] = {2};
    stat = nc_put_vara(ncid, s2_id, s2_startset, s2_countset, s2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i2_data[2] = {-2147483646, 2147483647} ;
    size_t i2_startset[1] = {0} ;
    size_t i2_countset[1] = {2};
    stat = nc_put_vara(ncid, i2_id, i2_startset, i2_countset, i2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f2_data[2] = {((float)-9.9999996e+35), ((float)9.9999996e+35)} ;
    size_t f2_startset[1] = {0} ;
    size_t f2_countset[1] = {2};
    stat = nc_put_vara(ncid, f2_id, f2_startset, f2_countset, f2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d2_data[2] = {((double)-1e+308), ((double)1e+308)} ;
    size_t d2_startset[1] = {0} ;
    size_t d2_countset[1] = {2};
    stat = nc_put_vara(ncid, d2_id, d2_startset, d2_countset, d2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c3_data = "\001\177." ;
    size_t c3_startset[1] = {0} ;
    size_t c3_countset[1] = {3};
    stat = nc_put_vara(ncid, c3_id, c3_startset, c3_countset, c3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b3_data[3] = {-128, 127, -1} ;
    size_t b3_startset[1] = {0} ;
    size_t b3_countset[1] = {3};
    stat = nc_put_vara(ncid, b3_id, b3_startset, b3_countset, b3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s3_data[3] = {-32768, 0, 32767} ;
    size_t s3_startset[1] = {0} ;
    size_t s3_countset[1] = {3};
    stat = nc_put_vara(ncid, s3_id, s3_startset, s3_countset, s3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i3_data[3] = {-2147483646, 0, 2147483647} ;
    size_t i3_startset[1] = {0} ;
    size_t i3_countset[1] = {3};
    stat = nc_put_vara(ncid, i3_id, i3_startset, i3_countset, i3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f3_data[3] = {((float)-9.9999996e+35), ((float)0), ((float)9.9999996e+35)} ;
    size_t f3_startset[1] = {0} ;
    size_t f3_countset[1] = {3};
    stat = nc_put_vara(ncid, f3_id, f3_startset, f3_countset, f3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d3_data[3] = {((double)-1e+308), ((double)0), ((double)1e+308)} ;
    size_t d3_startset[1] = {0} ;
    size_t d3_countset[1] = {3};
    stat = nc_put_vara(ncid, d3_id, d3_startset, d3_countset, d3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* cr1_data = "xy" ;
    size_t cr1_startset[2] = {0, 0} ;
    size_t cr1_countset[2] = {2, 1};
    stat = nc_put_vara(ncid, cr1_id, cr1_startset, cr1_countset, cr1_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char br2_data[4] = {-24, -26, -20, -22} ;
    size_t br2_startset[2] = {0, 0} ;
    size_t br2_countset[2] = {2, 2};
    stat = nc_put_vara(ncid, br2_id, br2_startset, br2_countset, br2_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short sr3_data[6] = {-375, -380, -385, -350, -355, -360} ;
    size_t sr3_startset[2] = {0, 0} ;
    size_t sr3_countset[2] = {2, 3};
    stat = nc_put_vara(ncid, sr3_id, sr3_startset, sr3_countset, sr3_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f11_data[1] = {((float)-2187)} ;
    size_t f11_startset[2] = {0, 0} ;
    size_t f11_countset[2] = {1, 1};
    stat = nc_put_vara(ncid, f11_id, f11_startset, f11_countset, f11_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d12_data[2] = {((double)-3000), ((double)-3010)} ;
    size_t d12_startset[2] = {0, 0} ;
    size_t d12_countset[2] = {1, 2};
    stat = nc_put_vara(ncid, d12_id, d12_startset, d12_countset, d12_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c13_data = "\tb\177" ;
    size_t c13_startset[2] = {0, 0} ;
    size_t c13_countset[2] = {1, 3};
    stat = nc_put_vara(ncid, c13_id, c13_startset, c13_countset, c13_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s21_data[2] = {-375, -350} ;
    size_t s21_startset[2] = {0, 0} ;
    size_t s21_countset[2] = {2, 1};
    stat = nc_put_vara(ncid, s21_id, s21_startset, s21_countset, s21_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i22_data[4] = {-24000, -24020, -23600, -23620} ;
    size_t i22_startset[2] = {0, 0} ;
    size_t i22_countset[2] = {2, 2};
    stat = nc_put_vara(ncid, i22_id, i22_startset, i22_countset, i22_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f23_data[6] = {((float)-2187), ((float)-2196), ((float)-2205), ((float)-2106), ((float)-2115), ((float)-2124)} ;
    size_t f23_startset[2] = {0, 0} ;
    size_t f23_countset[2] = {2, 3};
    stat = nc_put_vara(ncid, f23_id, f23_startset, f23_countset, f23_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c31_data = "+- " ;
    size_t c31_startset[2] = {0, 0} ;
    size_t c31_countset[2] = {3, 1};
    stat = nc_put_vara(ncid, c31_id, c31_startset, c31_countset, c31_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b32_data[6] = {-24, -26, -20, -22, -16, -18} ;
    size_t b32_startset[2] = {0, 0} ;
    size_t b32_countset[2] = {3, 2};
    stat = nc_put_vara(ncid, b32_id, b32_startset, b32_countset, b32_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s33_data[9] = {-375, -380, -385, -350, -355, -360, -325, -330, -335} ;
    size_t s33_startset[2] = {0, 0} ;
    size_t s33_countset[2] = {3, 3};
    stat = nc_put_vara(ncid, s33_id, s33_startset, s33_countset, s33_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short sr11_data[2] = {2500, 2375} ;
    size_t sr11_startset[3] = {0, 0, 0} ;
    size_t sr11_countset[3] = {2, 1, 1};
    stat = nc_put_vara(ncid, sr11_id, sr11_startset, sr11_countset, sr11_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int ir12_data[4] = {640000, 639980, 632000, 631980} ;
    size_t ir12_startset[3] = {0, 0, 0} ;
    size_t ir12_countset[3] = {2, 1, 2};
    stat = nc_put_vara(ncid, ir12_id, ir12_startset, ir12_countset, ir12_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float fr13_data[6] = {((float)26244), ((float)26235), ((float)26226), ((float)25515), ((float)25506), ((float)25497)} ;
    size_t fr13_startset[3] = {0, 0, 0} ;
    size_t fr13_countset[3] = {2, 1, 3};
    stat = nc_put_vara(ncid, fr13_id, fr13_startset, fr13_countset, fr13_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* cr21_data = "@DHL" ;
    size_t cr21_startset[3] = {0, 0, 0} ;
    size_t cr21_countset[3] = {2, 2, 1};
    stat = nc_put_vara(ncid, cr21_id, cr21_startset, cr21_countset, cr21_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char br22_data[8] = {64, 62, 68, 66, 56, 54, 60, 58} ;
    size_t br22_startset[3] = {0, 0, 0} ;
    size_t br22_countset[3] = {2, 2, 2};
    stat = nc_put_vara(ncid, br22_id, br22_startset, br22_countset, br22_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short sr23_data[12] = {2500, 2495, 2490, 2525, 2520, 2515, 2375, 2370, 2365, 2400, 2395, 2390} ;
    size_t sr23_startset[3] = {0, 0, 0} ;
    size_t sr23_countset[3] = {2, 2, 3};
    stat = nc_put_vara(ncid, sr23_id, sr23_startset, sr23_countset, sr23_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float fr31_data[6] = {((float)26244), ((float)26325), ((float)26406), ((float)25515), ((float)25596), ((float)25677)} ;
    size_t fr31_startset[3] = {0, 0, 0} ;
    size_t fr31_countset[3] = {2, 3, 1};
    stat = nc_put_vara(ncid, fr31_id, fr31_startset, fr31_countset, fr31_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double dr32_data[12] = {((double)40000), ((double)39990), ((double)40100), ((double)40090), ((double)40200), ((double)40190), ((double)39000), ((double)38990), ((double)39100), ((double)39090), ((double)39200), ((double)39190)} ;
    size_t dr32_startset[3] = {0, 0, 0} ;
    size_t dr32_countset[3] = {2, 3, 2};
    stat = nc_put_vara(ncid, dr32_id, dr32_startset, dr32_countset, dr32_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* cr33_data = "1\000\000two3\000\0004\000\0005\000\000six" ;
    size_t cr33_startset[3] = {0, 0, 0} ;
    size_t cr33_countset[3] = {2, 3, 3};
    stat = nc_put_vara(ncid, cr33_id, cr33_startset, cr33_countset, cr33_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c111_data = "@" ;
    size_t c111_startset[3] = {0, 0, 0} ;
    size_t c111_countset[3] = {1, 1, 1};
    stat = nc_put_vara(ncid, c111_id, c111_startset, c111_countset, c111_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b112_data[2] = {64, 62} ;
    size_t b112_startset[3] = {0, 0, 0} ;
    size_t b112_countset[3] = {1, 1, 2};
    stat = nc_put_vara(ncid, b112_id, b112_startset, b112_countset, b112_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s113_data[3] = {2500, 2495, 2490} ;
    size_t s113_startset[3] = {0, 0, 0} ;
    size_t s113_countset[3] = {1, 1, 3};
    stat = nc_put_vara(ncid, s113_id, s113_startset, s113_countset, s113_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f121_data[2] = {((float)26244), ((float)26325)} ;
    size_t f121_startset[3] = {0, 0, 0} ;
    size_t f121_countset[3] = {1, 2, 1};
    stat = nc_put_vara(ncid, f121_id, f121_startset, f121_countset, f121_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d122_data[4] = {((double)40000), ((double)39990), ((double)40100), ((double)40090)} ;
    size_t d122_startset[3] = {0, 0, 0} ;
    size_t d122_countset[3] = {1, 2, 2};
    stat = nc_put_vara(ncid, d122_id, d122_startset, d122_countset, d122_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c123_data = "one2\000\000" ;
    size_t c123_startset[3] = {0, 0, 0} ;
    size_t c123_countset[3] = {1, 2, 3};
    stat = nc_put_vara(ncid, c123_id, c123_startset, c123_countset, c123_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s131_data[3] = {2500, 2525, 2550} ;
    size_t s131_startset[3] = {0, 0, 0} ;
    size_t s131_countset[3] = {1, 3, 1};
    stat = nc_put_vara(ncid, s131_id, s131_startset, s131_countset, s131_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i132_data[6] = {640000, 639980, 640400, 640380, 640800, 640780} ;
    size_t i132_startset[3] = {0, 0, 0} ;
    size_t i132_countset[3] = {1, 3, 2};
    stat = nc_put_vara(ncid, i132_id, i132_startset, i132_countset, i132_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f133_data[9] = {((float)26244), ((float)26235), ((float)26226), ((float)26325), ((float)26316), ((float)26307), ((float)26406), ((float)26397), ((float)26388)} ;
    size_t f133_startset[3] = {0, 0, 0} ;
    size_t f133_countset[3] = {1, 3, 3};
    stat = nc_put_vara(ncid, f133_id, f133_startset, f133_countset, f133_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f211_data[2] = {((float)26244), ((float)25515)} ;
    size_t f211_startset[3] = {0, 0, 0} ;
    size_t f211_countset[3] = {2, 1, 1};
    stat = nc_put_vara(ncid, f211_id, f211_startset, f211_countset, f211_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    double d212_data[4] = {((double)40000), ((double)39990), ((double)39000), ((double)38990)} ;
    size_t d212_startset[3] = {0, 0, 0} ;
    size_t d212_countset[3] = {2, 1, 2};
    stat = nc_put_vara(ncid, d212_id, d212_startset, d212_countset, d212_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c213_data = "\000\000\000\000\000\000" ;
    size_t c213_startset[3] = {0, 0, 0} ;
    size_t c213_countset[3] = {2, 1, 3};
    stat = nc_put_vara(ncid, c213_id, c213_startset, c213_countset, c213_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s221_data[4] = {2500, 2525, 2375, 2400} ;
    size_t s221_startset[3] = {0, 0, 0} ;
    size_t s221_countset[3] = {2, 2, 1};
    stat = nc_put_vara(ncid, s221_id, s221_startset, s221_countset, s221_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i222_data[8] = {640000, 639980, 640400, 640380, 632000, 631980, 632400, 632380} ;
    size_t i222_startset[3] = {0, 0, 0} ;
    size_t i222_countset[3] = {2, 2, 2};
    stat = nc_put_vara(ncid, i222_id, i222_startset, i222_countset, i222_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f223_data[12] = {((float)26244), ((float)26235), ((float)26226), ((float)26325), ((float)26316), ((float)26307), ((float)25515), ((float)25506), ((float)25497), ((float)25596), ((float)25587), ((float)25578)} ;
    size_t f223_startset[3] = {0, 0, 0} ;
    size_t f223_countset[3] = {2, 2, 3};
    stat = nc_put_vara(ncid, f223_id, f223_startset, f223_countset, f223_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    char* c231_data = "@DHHLP" ;
    size_t c231_startset[3] = {0, 0, 0} ;
    size_t c231_countset[3] = {2, 3, 1};
    stat = nc_put_vara(ncid, c231_id, c231_startset, c231_countset, c231_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    signed char b232_data[12] = {64, 62, 68, 66, 72, 70, 56, 54, 60, 58, 64, 62} ;
    size_t b232_startset[3] = {0, 0, 0} ;
    size_t b232_countset[3] = {2, 3, 2};
    stat = nc_put_vara(ncid, b232_id, b232_startset, b232_countset, b232_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s233_data[18] = {2500, 2495, 2490, 2525, 2520, 2515, 2550, 2545, 2540, 2375, 2370, 2365, 2400, 2395, 2390, 2425, 2420, 2415} ;
    size_t s233_startset[3] = {0, 0, 0} ;
    size_t s233_countset[3] = {2, 3, 3};
    stat = nc_put_vara(ncid, s233_id, s233_startset, s233_countset, s233_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    short s311_data[3] = {2500, 2375, 2250} ;
    size_t s311_startset[3] = {0, 0, 0} ;
    size_t s311_countset[3] = {3, 1, 1};
    stat = nc_put_vara(ncid, s311_id, s311_startset, s311_countset, s311_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    int i312_data[6] = {640000, 639980, 632000, 631980, 624000, 623980} ;
    size_t i312_startset[3] = {0, 0, 0} ;
    size_t i312_countset[3] = {3, 1, 2};
    stat = nc_put_vara(ncid, i312_id, i312_startset, i312_countset, i312_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    float f313_data[9] = {((float)26244), ((float)26235), ((float)26226), ((float)25515), ((float)25506), ((float)25497), ((float)24786), ((float)24777), ((float)24768)} ;
    size_t f313_startset[3] = {0, 0, 0} ;
    size_t f313_countset[3] = {3, 1, 3};
    stat = nc_put_vara(ncid, f313_id, f313_startset, f313_countset, f313_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static double var_MINUS_name_MINUS_dashes_data[1] = {((double)-1)};
    stat = nc_put_var1(ncid, var_MINUS_name_MINUS_dashes_id, &count, var_MINUS_name_MINUS_dashes_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static double var_PERIOD_name_PERIOD_dots_data[1] = {((double)-2)};
    stat = nc_put_var1(ncid, var_PERIOD_name_PERIOD_dots_id, &count, var_PERIOD_name_PERIOD_dots_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static double var_PLUS_name_PLUS_plusses_data[1] = {((double)9.969209968386869e+36)};
    stat = nc_put_var1(ncid, var_PLUS_name_PLUS_plusses_id, &count, var_PLUS_name_PLUS_plusses_data);
    check_err(stat,__LINE__,__FILE__);
    }


    {
    size_t count = 0;
    static double var_ATSIGN_name_ATSIGN_ats_data[1] = {((double)9.969209968386869e+36)};
    stat = nc_put_var1(ncid, var_ATSIGN_name_ATSIGN_ats_id, &count, var_ATSIGN_name_ATSIGN_ats_data);
    check_err(stat,__LINE__,__FILE__);
    }


    stat = nc_close(ncid);
    check_err(stat,__LINE__,__FILE__);
    return 0;
}
