/* This is part of the netCDF package. Copyright 2005-2007 University
   Corporation for Atmospheric Research/Unidata. See COPYRIGHT file
   for conditions of use.

   Test attributes. 

   $Id$
*/

#include <nc_tests.h>
#include "netcdf.h"
#include <signal.h>

#define FILE_NAME "tst_atts.nc"
#define FILE_NAME2 "tst_atts_2.nc"
#define VAR1_NAME "Horace_Rumpole"
#define VAR2_NAME "Claude_Erskine-Brown"
#define VAR3_NAME "Phillida_Erskine-Brown_Q.C."
#define DIM1_NAME "Old_Bailey_case_number"
#define DIM1_LEN 10
#define DIM2_NAME "occupancy_in_chambers"
#define DIM2_LEN 15
#define ATT_INT_NAME "Old_Bailey_Room_Numbers"
#define ATT_DOUBLE_NAME "Equity_Court_Canteen_Charges"
#define ATT_SHORT_NAME "Ecclesiastical_Court_Appearences"
#define ATT_TEXT_NAME "Speech_to_Jury"
#define ATT_TEXT_NAME2 "Speech_to_She_Who_Must_be_Obeyed"
#define ATT_UCHAR_NAME "Number_of_current_briefs"
#define ATT_SCHAR_NAME "Slate_totals_at_Pomeroys_Wine_Bar"
#define ATT_USHORT_NAME "brief_no"
#define ATT_UINT_NAME "Orders_from_SWMBO"
#define ATT_INT64_NAME "judges_golf_score"
#define ATT_UINT64_NAME "Number_of_drinks_in_career_to_date"

/*
#define ATT_USHORT_NAME "Chamber_Gas_Electric_and_Telephone_Bill_Share"
*/
#define ATT_FLOAT_NAME "Average_Nanoseconds_for_Lose_Win_or_Appeal"
#define ATT_LEN 3

char speech[] = "Once more unto the breach, dear friends, once more;\n\
Or close the wall up with our English dead.\n\
In peace there's nothing so becomes a man\n\
As modest stillness and humility:\n\
But when the blast of war blows in our ears,\n\
Then imitate the action of the tiger;\n\
Stiffen the sinews, summon up the blood,\n\
Disguise fair nature with hard-favour'd rage;\n\
Then lend the eye a terrible aspect;\n\
Let pry through the portage of the head\n\
Like the brass cannon; let the brow o'erwhelm it\n\
As fearfully as doth a galled rock\n\
O'erhang and jutty his confounded base,\n\
Swill'd with the wild and wasteful ocean.\n\
Now set the teeth and stretch the nostril wide,\n\
Hold hard the breath and bend up every spirit\n\
To his full height. On, on, you noblest English.\n\
Whose blood is fet from fathers of war-proof!\n\
Fathers that, like so many Alexanders,\n\
Have in these parts from morn till even fought\n\
And sheathed their swords for lack of argument:\n\
Dishonour not your mothers; now attest\n\
That those whom you call'd fathers did beget you.\n\
Be copy now to men of grosser blood,\n\
And teach them how to war. And you, good yeoman,\n\
Whose limbs were made in England, show us here\n\
The mettle of your pasture; let us swear\n\
That you are worth your breeding; which I doubt not;\n\
For there is none of you so mean and base,\n\
That hath not noble lustre in your eyes.\n\
I see you stand like greyhounds in the slips,\n\
Straining upon the start. The game's afoot:\n\
Follow your spirit, and upon this charge\n\
Cry 'God for Harry, England, and Saint George!'";

/* Test the ordering of atts for a cmode. */
#define NUM_ATTS 8
#define ATT_MAX_NAME 25
int
tst_att_ordering(int cmode)
{
   int ncid;
   char name[NUM_ATTS][ATT_MAX_NAME + 1] = {"Gc", "Gb", "Gs", "Gi", "Gf", 
					    "Gd", "Gatt-name-dashes", "Gatt.name.dots"};
   int len[NUM_ATTS] = {0, 2, 3, 3, 3, 3, 1, 1};
   signed char b[2] = {-128, 127};
   short s[3] = {-32768, 0, 32767};
   int i[3] = {42, 0, -42};
   float f[3] = {42.0, -42.0, 42.0};
   double d[3] = {420.0, -420.0, 420.0};
   int att_name_dashes = -1, att_name_dots = -2;
   char name_in[NC_MAX_NAME];
   int j;

   /* Create a file with some global atts. */
   if (nc_create(FILE_NAME, cmode, &ncid)) ERR;
   if (nc_put_att_text(ncid, NC_GLOBAL, name[0], len[0], NULL)) ERR;      
   if (nc_put_att_schar(ncid, NC_GLOBAL, name[1], NC_BYTE, len[1], b)) ERR;      
   if (nc_put_att_short(ncid, NC_GLOBAL, name[2], NC_SHORT, len[2], s)) ERR;      
   if (nc_put_att_int(ncid, NC_GLOBAL, name[3], NC_INT, len[3], i)) ERR;       
   if (nc_put_att_float(ncid, NC_GLOBAL, name[4], NC_FLOAT, len[4], f)) ERR;       
   if (nc_put_att_double(ncid, NC_GLOBAL, name[5], NC_DOUBLE, len[5], d)) ERR;       
   if (nc_put_att_int(ncid, NC_GLOBAL, name[6], NC_INT, len[6], &att_name_dashes)) ERR;       
   if (nc_put_att_int(ncid, NC_GLOBAL, name[7], NC_INT, len[7], &att_name_dots)) ERR;       
   if (nc_close(ncid)) ERR;
      
   /* Reopen the file and check the order. */
   if (nc_open(FILE_NAME, 0, &ncid)) ERR;
   for (j = 0; j < NUM_ATTS; j++)
   {
      if (nc_inq_attname(ncid, NC_GLOBAL, j, name_in)) ERR;
      if (strcmp(name_in, name[j])) ERR;
   }

   /* Close up shop. */
   if (nc_close(ncid)) ERR;
   return err;
}

int
main(int argc, char **argv)
{

    signed char schar_in[ATT_LEN], schar_out[ATT_LEN] = {NC_MIN_BYTE, 1, NC_MAX_BYTE};
    unsigned char uchar_in[ATT_LEN], uchar_out[ATT_LEN] = {0, 128, NC_MAX_UBYTE};
    short short_in[ATT_LEN], short_out[ATT_LEN] = {NC_MIN_SHORT, -128, NC_MAX_SHORT};
    unsigned short ushort_in[ATT_LEN], ushort_out[ATT_LEN] = {0, 128, NC_MAX_USHORT};
    int int_in[ATT_LEN], int_out[ATT_LEN] = {-100000, 128, 100000};
    long long_in[ATT_LEN];
    unsigned int uint_in[ATT_LEN], uint_out[ATT_LEN] = {0, 128, NC_MAX_UINT};
    float float_in[ATT_LEN], float_out[ATT_LEN] = {-0.5, 0.25, 0.125};
    double double_in[ATT_LEN], double_out[ATT_LEN] = {-0.25, .5, 0.125};
    long long longlong_in[ATT_LEN], longlong_out[ATT_LEN] = {-3123456789LL, 128LL, 3123456789LL};
    unsigned long long ulonglong_in[ATT_LEN], ulonglong_out[ATT_LEN] = {0LL, 128LL, 3123456789LL};
	
   (void) signal(SIGFPE, SIG_IGN);
   printf("\n*** Testing netcdf-4 attribute functions.\n");
   printf("*** testing really simple global atts...");
#define NUM_SIMPLE_ATTS 9
   {
      int ncid;
      char name[NUM_SIMPLE_ATTS][ATT_MAX_NAME + 1] = {"Gc", "Gb", "Gs", "Gi", "Gf", 
						      "Gd", "G7", "G8", "G9"};
      char name_in[NC_MAX_NAME];
      int j;

      /* Create a file with some global atts. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLOBBER, &ncid)) ERR;
      for (j = 0; j < NUM_SIMPLE_ATTS; j++)
	 if (nc_put_att_int(ncid, NC_GLOBAL, name[j], NC_INT, 0, NULL)) ERR;      
      if (nc_close(ncid)) ERR;
      
      /* Reopen the file and check the order. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      for (j = 0; j < NUM_SIMPLE_ATTS; j++)
      {
	 if (nc_inq_attname(ncid, NC_GLOBAL, j, name_in)) ERR;
	 if (strcmp(name_in, name[j])) ERR;
      }

      /* Close up shop. */
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing simple global atts...");
   {      
      int ncid;
      nc_type att_type;
      size_t att_len;
      int i;

      char *speech_in;

      /* This won't work, because classic files can't create these types. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_put_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, NC_USHORT, ATT_LEN, 
			    ushort_out) != NC_ESTRICTNC3) ERR;      
      if (nc_put_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, NC_UINT, ATT_LEN, 
			  uint_out) != NC_ESTRICTNC3) ERR;      
      if (nc_put_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, NC_INT64, ATT_LEN, 
			      longlong_out) != NC_ESTRICTNC3) ERR;      
      if (nc_put_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, NC_UINT64, ATT_LEN, 
			       ulonglong_out) != NC_ESTRICTNC3) ERR;      
      /* But it's OK to put classic types like NC_INT converted from
       * supported C types, though there may be out-of-range errrors
       * for some values */
      if (nc_put_att_uint(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, 
			  uint_out) != NC_ERANGE) ERR;
      if (nc_put_att_longlong(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, 
			      longlong_out) != NC_ERANGE) ERR;
      if (nc_put_att_ulonglong(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, 
			       ulonglong_out) != NC_ERANGE) ERR;
      if (nc_put_att_ushort(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, 
			    ushort_out)) ERR;      
      /* restore to intended values for subsequent tests */
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, 
			    int_out)) ERR;      
      /* It should also be OK to read classic types converted into
       * supported C types. though the conversion may encounter
       * out-of-range values */
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_INT_NAME, uchar_in) != NC_ERANGE) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_INT_NAME, ushort_in) != NC_ERANGE) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_INT_NAME, uint_in) != NC_ERANGE) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_INT_NAME, longlong_in)) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_INT_NAME, ulonglong_in) != NC_ERANGE) ERR;

      if (nc_close(ncid)) ERR;

      /* Create a file with a global attribute of each type. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, strlen(speech)+1, speech)) ERR;      
      if (nc_put_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, NC_BYTE, ATT_LEN, schar_out)) ERR;      
      if (nc_put_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, NC_UBYTE, ATT_LEN, uchar_out)) ERR;
      if (nc_put_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, NC_SHORT, ATT_LEN, short_out)) ERR;      
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, int_out)) ERR;      
      if (nc_put_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, NC_FLOAT, ATT_LEN, float_out)) ERR;      
      if (nc_put_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, NC_DOUBLE, ATT_LEN, double_out)) ERR;      
      if (nc_put_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, NC_USHORT, ATT_LEN, ushort_out)) ERR;      
      if (nc_put_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, NC_UINT, ATT_LEN, uint_out)) ERR;      
      if (nc_put_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, NC_INT64, ATT_LEN, longlong_out)) ERR;      
      if (nc_put_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, NC_UINT64, ATT_LEN, ulonglong_out)) ERR;      
      if (nc_close(ncid)) ERR;

      /* Open the file and check attributes. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      /* Check text. */
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_TEXT_NAME, &att_type, &att_len))
	 ERR;
      if (att_type != NC_CHAR || att_len != strlen(speech) + 1) ERR;
      if (!(speech_in = malloc(att_len + 1))) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, speech_in)) ERR;      
      if (strcmp(speech, speech_in)) ERR;
      free(speech_in);
      /* Check numeric values. */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, schar_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != schar_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uchar_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, short_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (short_in[i] != short_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, int_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != int_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, float_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != float_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, double_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != double_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, ushort_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, uint_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != uint_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, longlong_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != longlong_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, ulonglong_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != ulonglong_out[i]) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing attribute data type conversions...");

   {
      int ncid;
      int i;

      /* Reopen the file and try different type conversions. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;

      /* No text conversions are allowed, and people who try them should
       * be locked up, away from decent folk! */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_TEXT_NAME, schar_in) != NC_ECHAR) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_TEXT_NAME, uchar_in) != NC_ECHAR) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_TEXT_NAME, short_in) != NC_ECHAR) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_TEXT_NAME, int_in) != NC_ECHAR) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_TEXT_NAME, float_in) != NC_ECHAR) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_TEXT_NAME, double_in) != NC_ECHAR) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_TEXT_NAME, ushort_in) != NC_ECHAR) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_TEXT_NAME, uint_in) != NC_ECHAR) ERR;
      if (nc_get_att_long(ncid, NC_GLOBAL, ATT_TEXT_NAME, long_in) != NC_ECHAR) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_TEXT_NAME, longlong_in) != NC_ECHAR) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_TEXT_NAME, ulonglong_in) != NC_ECHAR) ERR;

      /* Read all atts (except text) as double. */
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_SCHAR_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != schar_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_UCHAR_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_SHORT_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != short_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_USHORT_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_INT_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != int_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_UINT_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != uint_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_FLOAT_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != float_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != double_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_INT64_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (double_in[i] != longlong_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_UINT64_NAME, double_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (double_in[i] != ulonglong_out[i]) ERR;

      /* Read all atts (except text) as float. */
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_SCHAR_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != schar_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_UCHAR_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_SHORT_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != short_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_USHORT_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_INT_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != int_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_UINT_NAME, float_in)) ERR;
      /* Why is comparison failing in 32-bit? Looks like a gcc compiler error */
      /* for (i = 0; i < ATT_LEN; i++) */
      /* 	  if (float_in[i] != (float) uint_out[i]) ERR; */
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != float_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, float_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (float_in[i] != (float) double_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_INT64_NAME, float_in)) ERR;
      /* Why is comparison failing in 32-bit? */
      /* for (i = 0; i < ATT_LEN; i++) */
      /* 	  if (float_in[i] != (float) longlong_out[i]) ERR; */
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_UINT64_NAME, float_in)) ERR;
#if !defined(_WIN32) && !defined(_WIN64) 
	  for (i = 0; i < ATT_LEN; i++)
	  if (float_in[i] != (float) ulonglong_out[i]) ERR;
#endif
      /* Read all atts (except text) as int. */
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_SCHAR_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != schar_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_UCHAR_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) uchar_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_SHORT_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != short_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_USHORT_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != int_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_UINT_NAME, int_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) uint_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_FLOAT_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) float_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, int_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) double_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT64_NAME, int_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) longlong_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_UINT64_NAME, int_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (int_in[i] != (int) ulonglong_out[i]) ERR;

      /* Read all atts (except text) as uint. */
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_SCHAR_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != (unsigned int) schar_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uint_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_SHORT_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != (unsigned int) short_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_USHORT_NAME, uint_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_INT_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != (unsigned int) int_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, uint_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != uint_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_FLOAT_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != (unsigned int) float_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != (unsigned int) double_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_INT64_NAME, uint_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (uint_in[i] != (unsigned int) longlong_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_UINT64_NAME, uint_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (uint_in[i] != (unsigned int) ulonglong_out[i]) ERR;

      /* Read all atts (except text) as short. */
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SCHAR_NAME, short_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (short_in[i] != schar_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_UCHAR_NAME, short_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (short_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, short_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (short_in[i] != short_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_USHORT_NAME, short_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) ushort_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_INT_NAME, short_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) int_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_UINT_NAME, short_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) uint_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_FLOAT_NAME, short_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) float_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, short_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) double_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_INT64_NAME, short_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) longlong_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_UINT64_NAME, short_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (short_in[i] != (short) ulonglong_out[i]) ERR;

      /* Read all atts (except text) as ushort. */
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_SCHAR_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) schar_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_UCHAR_NAME, ushort_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_SHORT_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) short_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, ushort_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_INT_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) int_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_UINT_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) uint_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_FLOAT_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) float_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != (unsigned short) double_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_INT64_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ushort_in[i] != (unsigned short) longlong_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_UINT64_NAME, ushort_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ushort_in[i] != (unsigned short) ulonglong_out[i]) ERR;

      /* Read all atts (except text) as schar. */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, schar_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != schar_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (schar_in[i] != (signed char) uchar_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SHORT_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) short_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_USHORT_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) ushort_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_INT_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) int_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_UINT_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) uint_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_FLOAT_NAME, schar_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) float_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, schar_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != (signed char) double_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_INT64_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (schar_in[i] != (signed char) longlong_out[i]) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_UINT64_NAME, schar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (schar_in[i] != (signed char) ulonglong_out[i]) ERR;

      /* Read all atts (except text) as uchar. */
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) schar_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uchar_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_SHORT_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) short_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_USHORT_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) ushort_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_INT_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) int_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UINT_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) uint_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_FLOAT_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) float_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != (unsigned char) double_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_INT64_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (uchar_in[i] != (unsigned char) longlong_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UINT64_NAME, uchar_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (uchar_in[i] != (unsigned char) ulonglong_out[i]) ERR;

      /* Read all atts (except text) into long long variable. */
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_SCHAR_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != schar_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_UCHAR_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_SHORT_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != short_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_USHORT_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_INT_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != int_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_UINT_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != uint_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_FLOAT_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != (int)float_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != (int)double_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != longlong_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_UINT64_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (longlong_in[i] != (unsigned long long) ulonglong_out[i]) ERR;

      /* Read all atts (except text) into unsigned long long variable. */
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_SCHAR_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ulonglong_in[i] != (unsigned long long) schar_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_UCHAR_NAME, ulonglong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_SHORT_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ulonglong_in[i] != (unsigned long long) short_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_USHORT_NAME, ulonglong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_INT_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ulonglong_in[i] != (unsigned long long) int_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT_NAME, ulonglong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != uint_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_FLOAT_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != (unsigned long long) float_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != (unsigned long long) double_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_INT64_NAME, ulonglong_in) != NC_ERANGE) ERR;
      for (i = 0; i < ATT_LEN; i++)
	  if (ulonglong_in[i] != (unsigned long long) longlong_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, ulonglong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != ulonglong_out[i]) ERR;

      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing simple variable atts...");
   {
      int ncid, varid, dimids[2];
      nc_type att_type;
      size_t att_len;
      int i, v;

      char *speech_in;
      signed char schar_in[ATT_LEN], schar_out[ATT_LEN] = {NC_MIN_BYTE, 1, NC_MAX_BYTE};
      short short_in[ATT_LEN], short_out[ATT_LEN] = {NC_MIN_SHORT, -128, NC_MAX_SHORT};
      /*int int_in[ATT_LEN], int_out[ATT_LEN] = {NC_MIN_INT, 128, NC_MAX_INT};*/
      int int_in[ATT_LEN], int_out[ATT_LEN] = {-100000, 128, 100000};
      float float_in[ATT_LEN], float_out[ATT_LEN] = {.5, 0.25, 0.125};
      double double_in[ATT_LEN], double_out[ATT_LEN] = {0.25, .5, 0.125};

      /* Create a file with two vars, attaching to each an attribute of
       * each type. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimids[0])) ERR;
      if (nc_def_dim(ncid, DIM2_NAME, DIM2_LEN, &dimids[1])) ERR;
      if (nc_def_var(ncid, VAR1_NAME, NC_INT, 2, dimids, &varid)) ERR;
      if (nc_put_att_text(ncid, varid, ATT_TEXT_NAME, strlen(speech)+1, speech)) ERR;      
      if (nc_put_att_schar(ncid, varid, ATT_SCHAR_NAME, NC_BYTE, ATT_LEN, schar_out)) ERR;      
      if (nc_put_att_short(ncid, varid, ATT_SHORT_NAME, NC_SHORT, 3, short_out)) ERR;      
      if (nc_put_att_int(ncid, varid, ATT_INT_NAME, NC_INT, 3, int_out)) ERR;      
      if (nc_put_att_float(ncid, varid, ATT_FLOAT_NAME, NC_FLOAT, 3, float_out)) ERR;      
      if (nc_put_att_double(ncid, varid, ATT_DOUBLE_NAME, NC_DOUBLE, 3, double_out)) ERR;      
      if (nc_def_var(ncid, VAR2_NAME, NC_UINT, 2, dimids, &varid)) ERR;
      if (nc_put_att_text(ncid, varid, ATT_TEXT_NAME, strlen(speech)+1, speech)) ERR;      
      if (nc_put_att_schar(ncid, varid, ATT_SCHAR_NAME, NC_BYTE, ATT_LEN, schar_out)) ERR; 
      if (nc_put_att_short(ncid, varid, ATT_SHORT_NAME, NC_SHORT, 3, short_out)) ERR;           
      if (nc_put_att_int(ncid, varid, ATT_INT_NAME, NC_INT, 3, int_out)) ERR;      
      if (nc_put_att_float(ncid, varid, ATT_FLOAT_NAME, NC_FLOAT, 3, float_out)) ERR;      
      if (nc_put_att_double(ncid, varid, ATT_DOUBLE_NAME, NC_DOUBLE, 3, double_out)) ERR;      
      if (nc_close(ncid)) ERR;

      /* Open the file and check attributes. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      for (v=0; v<2; v++)
      {
	 if (nc_inq_att(ncid, v, ATT_TEXT_NAME, &att_type, &att_len)) ERR;
	 if (att_type != NC_CHAR || att_len != strlen(speech) + 1) ERR;
	 if (!(speech_in = malloc(att_len + 1))) ERR;
	 if (nc_get_att_text(ncid, v, ATT_TEXT_NAME, speech_in)) ERR;      
	 if (strcmp(speech, speech_in)) ERR;
	 free(speech_in);
	 if (nc_get_att_schar(ncid, v, ATT_SCHAR_NAME, schar_in)) ERR;      
	 for (i = 0; i < ATT_LEN; i++)
	    if (schar_in[i] != schar_out[i]) ERR;
	 if (nc_get_att_short(ncid, v, ATT_SHORT_NAME, short_in)) ERR;      
	 for (i = 0; i < ATT_LEN; i++)
	    if (short_in[i] != short_out[i]) ERR;
	 if (nc_get_att_int(ncid, v, ATT_INT_NAME, int_in)) ERR;      
	 for (i = 0; i < ATT_LEN; i++)
	    if (int_in[i] != int_out[i]) ERR;
	 if (nc_get_att_float(ncid, v, ATT_FLOAT_NAME, float_in)) ERR;      
	 for (i = 0; i < ATT_LEN; i++)
	    if (float_in[i] != float_out[i]) ERR;
	 if (nc_get_att_double(ncid, v, ATT_DOUBLE_NAME, double_in)) ERR;      
	 for (i = 0; i < ATT_LEN; i++)
	    if (double_in[i] != double_out[i]) ERR;
      }
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   printf("*** testing zero-length attributes...");
   {
      int ncid;

      /* Create a file with a global attribute of each type of zero length. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, 0, NULL)) ERR;
      if (nc_put_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, NC_BYTE, 0, NULL)) ERR;
      if (nc_put_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, NC_UBYTE, 0, uchar_out)) ERR;
      if (nc_put_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, NC_SHORT, 0, NULL)) ERR;
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, 0, NULL)) ERR;
      if (nc_put_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, NC_FLOAT, 0, NULL)) ERR;
      if (nc_put_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, NC_DOUBLE, 0, NULL)) ERR;
      if (nc_close(ncid)) ERR;
   }

   /* Make sure we can read all these zero-length atts. */
   {
      int ncid;
      signed char schar_in[ATT_LEN];
      unsigned char uchar_in[ATT_LEN];
      short short_in[ATT_LEN];
      int int_in[ATT_LEN];
      float float_in[ATT_LEN];
      double double_in[ATT_LEN];
      size_t len;
      nc_type xtype;

      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, NULL)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_TEXT_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_CHAR) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, schar_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_SCHAR_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_BYTE) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uchar_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_UCHAR_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_UBYTE) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, short_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_SHORT_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_SHORT) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, int_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_INT_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_INT) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, float_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_FLOAT_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_FLOAT) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, double_in)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, &xtype, &len)) ERR;
      if (len || xtype != NC_DOUBLE) ERR;
      /* Conversions no longer result in range errors, since there's no data. */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_FLOAT_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_INT_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SHORT_NAME, schar_in)) ERR;
      if (nc_close(ncid)) ERR;
   }

   SUMMARIZE_ERR;
   printf("*** testing zero-length attributes and redef...");
   {
      int ncid;
      signed char schar_in[ATT_LEN];
      unsigned char uchar_in[ATT_LEN];
      short short_in[ATT_LEN];
      int int_in[ATT_LEN];
      float float_in[ATT_LEN];
      double double_in[ATT_LEN];


      /* Create a file with a global attribute of each type of zero length. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_enddef(ncid)) ERR;
      if (nc_redef(ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, 0, NULL)) ERR;
      if (nc_put_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, NC_BYTE, 0, NULL)) ERR;
      if (nc_put_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, NC_UBYTE, 0, uchar_out)) ERR;
      if (nc_put_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, NC_SHORT, 0, NULL)) ERR;
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, 0, NULL)) ERR;
      if (nc_put_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, NC_FLOAT, 0, NULL)) ERR;
      if (nc_put_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, NC_DOUBLE, 0, NULL)) ERR;
      if (nc_close(ncid)) ERR;

      /* Make sure we can read all these zero-length atts added during a
       * redef. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, NULL)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, schar_in)) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uchar_in)) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, short_in)) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, int_in)) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, float_in)) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, double_in)) ERR;
      /* Conversions no longer result in range errors, since there's no data. */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_FLOAT_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_INT_NAME, schar_in)) ERR;
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SHORT_NAME, schar_in)) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;

   printf("*** testing attribute deletes and renames...");
   {
      int ncid, varid, dimids[2];
      nc_type att_type;
      size_t att_len;
      char *speech_in;
      char name_in[NC_MAX_NAME + 1];
      int attid_in, natts_in;
      int int_out[ATT_LEN] = {-100000, 128, 100000};

      /* Create a file with a global attribute. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, strlen(speech)+1, 
			  speech)) ERR;      
      if (nc_close(ncid)) ERR;
      
      /* Rename it. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_attid(ncid, NC_GLOBAL, ATT_TEXT_NAME, &attid_in)) ERR;
      if (attid_in != 0) ERR;
      if (nc_inq_attname(ncid, NC_GLOBAL, attid_in, name_in)) ERR;
      if (strcmp(name_in, ATT_TEXT_NAME)) ERR;
      if (nc_rename_att(ncid, NC_GLOBAL, ATT_TEXT_NAME, ATT_TEXT_NAME2)) ERR;      
      if (nc_inq_attname(ncid, NC_GLOBAL, attid_in, name_in)) ERR;
      if (strcmp(name_in, ATT_TEXT_NAME2)) ERR;
      if (nc_close(ncid)) ERR;

      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_TEXT_NAME2, &att_type, &att_len)) ERR;
      if (att_type != NC_CHAR || att_len != strlen(speech) + 1) ERR;
      if (!(speech_in = malloc(att_len + 1))) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME2, speech_in)) ERR;      
      if (strcmp(speech, speech_in)) ERR;
      free(speech_in);
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, speech_in) != NC_ENOTATT) ERR;      
      if (nc_close(ncid)) ERR;

      /* Now delete the att. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_del_att(ncid, NC_GLOBAL, ATT_TEXT_NAME2)) ERR;
      if (nc_close(ncid)) ERR;

      /* Now create a file with a variable, which has an att. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, strlen(speech)+1, speech)) ERR;      
      if (nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimids[0])) ERR;
      if (nc_def_dim(ncid, DIM2_NAME, DIM2_LEN, &dimids[1])) ERR;
      if (nc_def_var(ncid, VAR1_NAME, NC_INT, 2, dimids, &varid)) ERR;
      if (nc_put_att_int(ncid, varid, ATT_INT_NAME, NC_INT, 3, int_out)) ERR;      
      if (nc_close(ncid)) ERR;
      
      /* Reopen the file and delete it. Make sure it's gone. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_del_att(ncid, 0, ATT_INT_NAME)) ERR;
      if (nc_close(ncid)) ERR;

      /* Reopen the file and readd the attribute. Enddef and redef,
       * and delete it, then check to make sure it's gone. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_put_att_int(ncid, varid, ATT_INT_NAME, NC_INT, 3, int_out)) ERR;      
      if (nc_enddef(ncid)) ERR;
      if (nc_redef(ncid)) ERR;
      if (nc_del_att(ncid, 0, ATT_INT_NAME)) ERR;
      if (nc_inq_varnatts(ncid, 0, &natts_in)) ERR;
      if (natts_in != 0) ERR;
      if (nc_close(ncid)) ERR;
   }

   SUMMARIZE_ERR;
   printf("*** testing attribute create order...");

#define ATT0 "Maturin"
#define ATT1 "Aubery"
   {
      int ncid, varid, dimids[2];
      int attid_in;
      const int number = 42;

      /* Create a file with several global attributes. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT0, NC_INT, 1, &number)) ERR;
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT1, NC_INT, 1, &number)) ERR;
      if (nc_close(ncid)) ERR;
      
      /* Open it and check the order. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_attid(ncid, NC_GLOBAL, ATT0, &attid_in)) ERR;
      if (attid_in != 0) ERR;
      if (nc_inq_attid(ncid, NC_GLOBAL, ATT1, &attid_in)) ERR;
      if (attid_in != 1) ERR;
      if (nc_close(ncid)) ERR;

      /* Now create a file with a variable, which has two atts. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_def_dim(ncid, DIM1_NAME, DIM1_LEN, &dimids[0])) ERR;
      if (nc_def_dim(ncid, DIM2_NAME, DIM2_LEN, &dimids[1])) ERR;
      if (nc_def_var(ncid, VAR1_NAME, NC_INT, 2, dimids, &varid)) ERR;
      if (nc_put_att_int(ncid, varid, ATT0, NC_INT, 1, &number)) ERR;
      if (nc_put_att_int(ncid, varid, ATT1, NC_INT, 1, &number)) ERR;
      if (nc_close(ncid)) ERR;
      
      /* Reopen the file and check the order of the attributes on the var. */
      if (nc_open(FILE_NAME, NC_WRITE, &ncid)) ERR;
      if (nc_inq_attid(ncid, 0, ATT0, &attid_in)) ERR;
      if (attid_in != 0) ERR;
      if (nc_inq_attid(ncid, 0, ATT1, &attid_in)) ERR;
      if (attid_in != 1) ERR;
      if (nc_close(ncid)) ERR;
   }

   SUMMARIZE_ERR;
   printf("*** testing attribute ordering some more...");

#define VAR_NAME "i"
#define A1_NAME "i"      
#define A2_NAME "f"      
#define A3_NAME "d"      
#define A1_LEN 3
#define A2_LEN 4
#define A3_LEN 5
   {
      int ncid;
      int varid, natts, nvars;
      double dvalue[] = {999.99, 999.99, 999.99, 999.99, 999.99};
      int varids[1];
      char name_in[NC_MAX_NAME + 1];

      /* Create a file with one var, and attach three atts to it. */
      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL, &ncid)) ERR;
      if (nc_def_var(ncid, VAR_NAME, NC_INT, 0, NULL, &varid)) ERR;
      if (nc_put_att_double(ncid, varid, A1_NAME, NC_INT, A1_LEN, dvalue)) ERR;      
      if (nc_put_att_double(ncid, varid, A2_NAME, NC_INT, A2_LEN, dvalue)) ERR;      
      if (nc_put_att_double(ncid, varid, A3_NAME, NC_INT, A3_LEN, dvalue)) ERR;      
      if (nc_close(ncid)) ERR;
      
      /* Reopen the file and check. */
      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_inq_varids(ncid, &nvars, varids)) ERR;
      if (nvars != 1 || varids[0] != 0) ERR;
      if (nc_inq_varnatts(ncid, 0, &natts)) ERR;
      if (natts != 3) ERR;
      if (nc_inq_attname(ncid, 0, 0, name_in)) ERR;
      if (strcmp(name_in, A1_NAME)) ERR;
      if (nc_inq_attname(ncid, 0, 1, name_in)) ERR;
      if (strcmp(name_in, A2_NAME)) ERR;
      if (nc_inq_attname(ncid, 0, 2, name_in)) ERR;
      if (strcmp(name_in, A3_NAME)) ERR;

      /* Close up shop. */
      if (nc_close(ncid)) ERR;
   }

   SUMMARIZE_ERR;
   printf("*** testing attribute ordering even more...");

   /* Test the ordering of atts for each cmode. */
   if (tst_att_ordering(NC_CLOBBER)) ERR;
   if (tst_att_ordering(NC_CLOBBER|NC_64BIT_OFFSET)) ERR;
   if (tst_att_ordering(NC_CLOBBER|NC_NETCDF4)) ERR;
   if (tst_att_ordering(NC_CLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL)) ERR;

   SUMMARIZE_ERR;
   printf("*** testing attributes and enddef/redef...");

#define ATT_1 "a"
#define ATT_2 "b"
#define ATT_3 "c"
   {
      int ncid, att = 1;

      if (nc_create(FILE_NAME, NC_NETCDF4|NC_CLASSIC_MODEL|NC_CLOBBER, &ncid)) ERR;
      if (nc_enddef(ncid)) ERR;
      if (nc_redef(ncid)) ERR;
      if (nc_put_att(ncid, NC_GLOBAL, ATT_1, NC_INT, 1, &att)) ERR;
      if (nc_put_att(ncid, NC_GLOBAL, ATT_2, NC_INT, 1, &att)) ERR;
      if (nc_put_att(ncid, NC_GLOBAL, ATT_3, NC_INT, 1, &att)) ERR;

      if (nc_close(ncid)) ERR;

      if (nc_open(FILE_NAME, 0, &ncid)) ERR;
      if (nc_close(ncid)) ERR;
   }

   SUMMARIZE_ERR;
   printf("*** testing copy of simple global atts...");
   {      
      int ncid, ncid2;
      nc_type att_type;
      size_t att_len;
      int i;

      char *speech_in;
      signed char schar_in[ATT_LEN], schar_out[ATT_LEN] = {NC_MIN_BYTE, 1, NC_MAX_BYTE};
      unsigned char uchar_in[ATT_LEN], uchar_out[ATT_LEN] = {0, 128, NC_MAX_CHAR};
      short short_in[ATT_LEN], short_out[ATT_LEN] = {NC_MIN_SHORT, -128, NC_MAX_SHORT};
      int int_in[ATT_LEN], int_out[ATT_LEN] = {-100000, 128, 100000};
      float float_in[ATT_LEN], float_out[ATT_LEN] = {.5, 0.25, 0.125};
      double double_in[ATT_LEN], double_out[ATT_LEN] = {0.25, .5, 0.125};
      unsigned short ushort_in[ATT_LEN], ushort_out[ATT_LEN] = {0, 128, NC_MAX_USHORT};
      unsigned int uint_in[ATT_LEN], uint_out[ATT_LEN] = {0, 128, NC_MAX_UINT};
      unsigned long long ulonglong_in[ATT_LEN], ulonglong_out[ATT_LEN] = {0, 128, 18446744073709551612ULL};
      long long longlong_in[ATT_LEN], longlong_out[ATT_LEN] = {NC_MIN_INT64, 128, NC_MAX_INT64};

      /* Create a file with a global attribute of each type. */
      if (nc_create(FILE_NAME, NC_NETCDF4, &ncid)) ERR;
      if (nc_put_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, strlen(speech)+1, speech)) ERR;      
      if (nc_put_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, NC_BYTE, ATT_LEN, schar_out)) ERR;      
      if (nc_put_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, NC_UBYTE, ATT_LEN, uchar_out)) ERR;
      if (nc_put_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, NC_SHORT, ATT_LEN, short_out)) ERR;      
      if (nc_put_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, NC_INT, ATT_LEN, int_out)) ERR;      
      if (nc_put_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, NC_FLOAT, ATT_LEN, float_out)) ERR;      
      if (nc_put_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, NC_DOUBLE, ATT_LEN, double_out)) ERR;      
      if (nc_put_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, NC_USHORT, ATT_LEN, ushort_out)) ERR;      
      if (nc_put_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, NC_UINT, ATT_LEN, uint_out)) ERR;      
      if (nc_put_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, NC_INT64, ATT_LEN, longlong_out)) ERR;      
      if (nc_put_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, NC_UINT64, ATT_LEN, ulonglong_out)) ERR;      

      /* Create another file and copy all the attributes. */
      if (nc_create(FILE_NAME2, NC_NETCDF4, &ncid2)) ERR;      
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_TEXT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_SCHAR_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_UCHAR_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_SHORT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_INT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_FLOAT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_USHORT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_UINT_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_INT64_NAME, ncid2, NC_GLOBAL)) ERR;
      if (nc_copy_att(ncid, NC_GLOBAL, ATT_UINT64_NAME, ncid2, NC_GLOBAL)) ERR;

      /* Close both files. */
      if (nc_close(ncid)) ERR;
      if (nc_close(ncid2)) ERR;

      /* Open the file and check attributes. */
      if (nc_open(FILE_NAME2, 0, &ncid)) ERR;
      /* Check text. */
      if (nc_inq_att(ncid, NC_GLOBAL, ATT_TEXT_NAME, &att_type, &att_len)) ERR;
      if (att_type != NC_CHAR || att_len != strlen(speech) + 1) ERR;
      if (!(speech_in = malloc(att_len + 1))) ERR;
      if (nc_get_att_text(ncid, NC_GLOBAL, ATT_TEXT_NAME, speech_in)) ERR;      
      if (strcmp(speech, speech_in)) ERR;
      free(speech_in);
      /* Check numeric values. */
      if (nc_get_att_schar(ncid, NC_GLOBAL, ATT_SCHAR_NAME, schar_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (schar_in[i] != schar_out[i]) ERR;
      if (nc_get_att_uchar(ncid, NC_GLOBAL, ATT_UCHAR_NAME, uchar_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uchar_in[i] != uchar_out[i]) ERR;
      if (nc_get_att_short(ncid, NC_GLOBAL, ATT_SHORT_NAME, short_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (short_in[i] != short_out[i]) ERR;
      if (nc_get_att_int(ncid, NC_GLOBAL, ATT_INT_NAME, int_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (int_in[i] != int_out[i]) ERR;
      if (nc_get_att_float(ncid, NC_GLOBAL, ATT_FLOAT_NAME, float_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (float_in[i] != float_out[i]) ERR;
      if (nc_get_att_double(ncid, NC_GLOBAL, ATT_DOUBLE_NAME, double_in)) ERR;      
      for (i = 0; i < ATT_LEN; i++)
	 if (double_in[i] != double_out[i]) ERR;
      if (nc_get_att_ushort(ncid, NC_GLOBAL, ATT_USHORT_NAME, ushort_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ushort_in[i] != ushort_out[i]) ERR;
      if (nc_get_att_uint(ncid, NC_GLOBAL, ATT_UINT_NAME, uint_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (uint_in[i] != uint_out[i]) ERR;
      if (nc_get_att_longlong(ncid, NC_GLOBAL, ATT_INT64_NAME, longlong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (longlong_in[i] != longlong_out[i]) ERR;
      if (nc_get_att_ulonglong(ncid, NC_GLOBAL, ATT_UINT64_NAME, ulonglong_in)) ERR;
      for (i = 0; i < ATT_LEN; i++)
	 if (ulonglong_in[i] != ulonglong_out[i]) ERR;
      if (nc_close(ncid)) ERR;
   }
   SUMMARIZE_ERR;
   FINAL_RESULTS;
}


