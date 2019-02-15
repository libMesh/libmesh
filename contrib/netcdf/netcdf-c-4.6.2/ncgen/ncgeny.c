/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         ncgparse
#define yylex           ncglex
#define yyerror         ncgerror
#define yydebug         ncgdebug
#define yynerrs         ncgnerrs

#define yylval          ncglval
#define yychar          ncgchar

/* Copy the first part of user declarations.  */
#line 11 "ncgen.y" /* yacc.c:339  */

/*
static char SccsId[] = "$Id: ncgen.y,v 1.42 2010/05/18 21:32:46 dmh Exp $";
*/
#include        "includes.h"
#include        "netcdf_aux.h"
#include        "ncgeny.h"
#include        "ncgen.h"
#ifdef USE_NETCDF4
#include        "ncfilter.h"
#endif

/* Following are in ncdump (for now)*/
/* Need some (unused) definitions to get it to compile */
#define ncatt_t void*
#define ncvar_t void
#include "nctime.h"

/* parser controls */
#define YY_NO_INPUT 1

/* True if string a equals string b*/
#ifndef NCSTREQ
#define NCSTREQ(a, b)     (*(a) == *(b) && strcmp((a), (b)) == 0)
#endif
#define VLENSIZE  (sizeof(nc_vlen_t))
#define MAXFLOATDIM 4294967295.0

/* mnemonics */
typedef enum Attrkind {ATTRVAR, ATTRGLOBAL, DONTKNOW} Attrkind;

#define ISCONST 1
#define ISLIST 0

typedef nc_vlen_t vlen_t;

/* We retain the old representation of the symbol list
   as a linked list.
*/
List* symlist;

/* Track rootgroup separately*/
Symbol* rootgroup;

/* Track the group sequence */
static List* groupstack;

/* Provide a separate sequence for accumulating values
   during the parse.
*/
static List* stack;

/* track homogeneity of types for data lists*/
static nc_type consttype;

/* Misc. */
static int stackbase;
static int stacklen;
static int count;
static int opaqueid; /* counter for opaque constants*/
static int arrayuid; /* counter for pseudo-array types*/

char* primtypenames[PRIMNO] = {
"nat",
"byte", "char", "short",
"int", "float", "double",
"ubyte", "ushort", "uint",
"int64", "uint64",
"string"
};

static int GLOBAL_SPECIAL = _NCPROPS_FLAG
                            | _ISNETCDF4_FLAG
                            | _SUPERBLOCK_FLAG
                            | _FORMAT_FLAG ;

/*Defined in ncgen.l*/
extern int lineno;              /* line number for error messages */
extern Bytebuffer* lextext;           /* name or string with escapes removed */

extern double double_val;       /* last double value read */
extern float float_val;         /* last float value read */
extern long long int64_val;         /* last int64 value read */
extern int int32_val;             /* last int32 value read */
extern short int16_val;         /* last short value read */
extern unsigned long long uint64_val;         /* last int64 value read */
extern unsigned int uint32_val;             /* last int32 value read */
extern unsigned short uint16_val;         /* last short value read */
extern char char_val;           /* last char value read */
extern signed char byte_val;    /* last byte value read */
extern unsigned char ubyte_val;    /* last byte value read */

/* Track definitions of dims, types, attributes, and vars*/
List* grpdefs;
List* dimdefs;
List* attdefs; /* variable-specific attributes*/
List* gattdefs; /* global attributes only*/
List* xattdefs; /* unknown attributes*/
List* typdefs;
List* vardefs;
List* condefs; /* non-dimension constants used in type defs*/
List* tmp;

/* Forward */
static NCConstant* makeconstdata(nc_type);
static NCConstant* evaluate(Symbol* fcn, Datalist* arglist);
static NCConstant* makeenumconstref(Symbol*);
static void addtogroup(Symbol*);
static Symbol* currentgroup(void);
static Symbol* createrootgroup(const char*);
static Symbol* creategroup(Symbol*);
static int dupobjectcheck(nc_class,Symbol*);
static void setpathcurrent(Symbol* sym);
static Symbol* makeattribute(Symbol*,Symbol*,Symbol*,Datalist*,Attrkind);
static Symbol* makeprimitivetype(nc_type i);
static Symbol* makespecial(int tag, Symbol* vsym, Symbol* tsym, void* data, int isconst);
static int containsfills(Datalist* list);
static void vercheck(int ncid);
static long long extractint(NCConstant* con);
#ifdef USE_NETCDF4
static int parsefilterflag(const char* sdata0, Specialdata* special);
#endif

int yylex(void);

#ifndef NO_STDARG
static void yyerror(const char *fmt, ...);
#else
static void yyerror(fmt,va_alist) const char* fmt; va_dcl;
#endif

/* Extern */
extern int lex_init(void);


#line 210 "ncgeny.c" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 1
#endif

/* In a future release of Bison, this section will be replaced
   by #include "ncgeny.h".  */
#ifndef YY_NCG_NCGEN_TAB_H_INCLUDED
# define YY_NCG_NCGEN_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int ncgdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    NC_UNLIMITED_K = 258,
    CHAR_K = 259,
    BYTE_K = 260,
    SHORT_K = 261,
    INT_K = 262,
    FLOAT_K = 263,
    DOUBLE_K = 264,
    UBYTE_K = 265,
    USHORT_K = 266,
    UINT_K = 267,
    INT64_K = 268,
    UINT64_K = 269,
    STRING_K = 270,
    IDENT = 271,
    TERMSTRING = 272,
    CHAR_CONST = 273,
    BYTE_CONST = 274,
    SHORT_CONST = 275,
    INT_CONST = 276,
    INT64_CONST = 277,
    UBYTE_CONST = 278,
    USHORT_CONST = 279,
    UINT_CONST = 280,
    UINT64_CONST = 281,
    FLOAT_CONST = 282,
    DOUBLE_CONST = 283,
    DIMENSIONS = 284,
    VARIABLES = 285,
    NETCDF = 286,
    DATA = 287,
    TYPES = 288,
    COMPOUND = 289,
    ENUM = 290,
    OPAQUE_ = 291,
    OPAQUESTRING = 292,
    GROUP = 293,
    PATH = 294,
    FILLMARKER = 295,
    NIL = 296,
    _FILLVALUE = 297,
    _FORMAT = 298,
    _STORAGE = 299,
    _CHUNKSIZES = 300,
    _DEFLATELEVEL = 301,
    _SHUFFLE = 302,
    _ENDIANNESS = 303,
    _NOFILL = 304,
    _FLETCHER32 = 305,
    _NCPROPS = 306,
    _ISNETCDF4 = 307,
    _SUPERBLOCK = 308,
    _FILTER = 309,
    DATASETID = 310
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 149 "ncgen.y" /* yacc.c:355  */

Symbol* sym;
unsigned long  size; /* allow for zero size to indicate e.g. UNLIMITED*/
long           mark; /* track indices into the sequence*/
int            nctype; /* for tracking attribute list type*/
Datalist*      datalist;
NCConstant*    constant;

#line 315 "ncgeny.c" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE ncglval;

int ncgparse (void);

#endif /* !YY_NCG_NCGEN_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 332 "ncgeny.c" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  5
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   393

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  65
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  67
/* YYNRULES -- Number of rules.  */
#define YYNRULES  153
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  262

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   310

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      61,    62,    63,     2,    59,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    64,    58,
       2,    60,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    56,     2,    57,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   231,   231,   237,   239,   246,   253,   253,   256,   265,
     255,   270,   271,   272,   276,   276,   278,   288,   288,   291,
     292,   293,   294,   297,   297,   300,   330,   332,   349,   358,
     370,   384,   417,   418,   421,   435,   436,   437,   438,   439,
     440,   441,   442,   443,   444,   445,   446,   449,   450,   451,
     454,   455,   458,   458,   460,   461,   465,   473,   483,   495,
     496,   497,   500,   501,   504,   504,   506,   531,   535,   539,
     568,   569,   572,   573,   577,   591,   595,   600,   629,   630,
     634,   635,   640,   650,   670,   681,   692,   711,   718,   718,
     721,   723,   725,   727,   729,   738,   749,   751,   753,   755,
     757,   759,   761,   763,   765,   767,   769,   774,   781,   790,
     791,   792,   795,   796,   799,   803,   804,   808,   812,   813,
     818,   819,   823,   824,   825,   826,   827,   828,   832,   836,
     840,   842,   847,   848,   849,   850,   851,   852,   853,   854,
     855,   856,   857,   858,   862,   863,   867,   869,   871,   873,
     878,   882,   883,   889
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "NC_UNLIMITED_K", "CHAR_K", "BYTE_K",
  "SHORT_K", "INT_K", "FLOAT_K", "DOUBLE_K", "UBYTE_K", "USHORT_K",
  "UINT_K", "INT64_K", "UINT64_K", "STRING_K", "IDENT", "TERMSTRING",
  "CHAR_CONST", "BYTE_CONST", "SHORT_CONST", "INT_CONST", "INT64_CONST",
  "UBYTE_CONST", "USHORT_CONST", "UINT_CONST", "UINT64_CONST",
  "FLOAT_CONST", "DOUBLE_CONST", "DIMENSIONS", "VARIABLES", "NETCDF",
  "DATA", "TYPES", "COMPOUND", "ENUM", "OPAQUE_", "OPAQUESTRING", "GROUP",
  "PATH", "FILLMARKER", "NIL", "_FILLVALUE", "_FORMAT", "_STORAGE",
  "_CHUNKSIZES", "_DEFLATELEVEL", "_SHUFFLE", "_ENDIANNESS", "_NOFILL",
  "_FLETCHER32", "_NCPROPS", "_ISNETCDF4", "_SUPERBLOCK", "_FILTER",
  "DATASETID", "'{'", "'}'", "';'", "','", "'='", "'('", "')'", "'*'",
  "':'", "$accept", "ncdesc", "datasetid", "rootgroup", "groupbody",
  "subgrouplist", "namedgroup", "$@1", "$@2", "typesection", "typedecls",
  "typename", "type_or_attr_decl", "typedecl", "optsemicolon", "enumdecl",
  "enumidlist", "enumid", "opaquedecl", "vlendecl", "compounddecl",
  "fields", "field", "primtype", "dimsection", "dimdecls",
  "dim_or_attr_decl", "dimdeclist", "dimdecl", "dimd", "vasection",
  "vadecls", "vadecl_or_attr", "vardecl", "varlist", "varspec", "dimspec",
  "dimlist", "dimref", "fieldlist", "fieldspec", "fielddimspec",
  "fielddimlist", "fielddim", "varref", "typeref", "type_var_ref",
  "attrdecllist", "attrdecl", "path", "datasection", "datadecls",
  "datadecl", "datalist", "datalist0", "datalist1", "dataitem",
  "constdata", "econstref", "function", "arglist", "simpleconstant",
  "intlist", "constint", "conststring", "constbool", "ident", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   123,   125,    59,    44,
      61,    40,    41,    42,    58
};
# endif

#define YYPACT_NINF -145

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-145)))

#define YYTABLE_NINF -108

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -10,   -19,    54,  -145,    12,  -145,   219,  -145,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,    -8,  -145,  -145,   354,    -9,    26,    11,  -145,
    -145,    15,    28,    31,    32,    36,   -11,    34,   130,   183,
      70,   219,    83,    83,    45,    52,   301,    85,  -145,  -145,
      -1,    42,    46,    49,    50,    51,    57,    58,    59,    60,
      61,    85,    53,   183,  -145,  -145,    64,    64,    64,    64,
      88,   255,    66,   219,    95,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,   301,  -145,    68,  -145,  -145,  -145,  -145,  -145,
    -145,  -145,    71,    75,    73,    77,   301,    83,    52,    52,
      45,    83,    45,    45,    83,   301,    79,  -145,   117,  -145,
    -145,  -145,  -145,  -145,  -145,    85,    78,  -145,   219,    82,
      86,  -145,    84,  -145,    87,   219,   118,    30,   301,   355,
    -145,   301,   301,    68,  -145,    90,  -145,  -145,  -145,  -145,
    -145,  -145,  -145,    68,   354,    93,    96,    97,   100,  -145,
      85,    35,   219,   102,  -145,   354,  -145,   354,  -145,  -145,
    -145,   -28,  -145,   219,    68,    68,    52,   291,   103,    85,
    -145,    85,    85,    85,  -145,  -145,  -145,  -145,  -145,   104,
    -145,   105,  -145,   -32,   108,  -145,   354,   107,   355,  -145,
    -145,  -145,  -145,   115,  -145,   123,  -145,   125,  -145,    38,
    -145,   140,  -145,  -145,    85,     3,  -145,   301,   143,  -145,
    -145,   166,  -145,    85,    -7,  -145,  -145,    85,    52,  -145,
     145,    17,  -145,  -145,    68,  -145,   106,  -145,  -145,  -145,
      23,  -145,  -145,  -145,     3,  -145,   219,    -7,  -145,  -145,
    -145,  -145
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     3,     0,     1,    88,     2,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
     153,   108,     0,     6,    87,     0,    85,    11,     0,    86,
     107,     0,     0,     0,     0,     0,     0,     0,     0,    12,
      47,    88,     0,     0,     0,     0,   117,     0,     4,     7,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    13,    14,    17,    23,    23,    23,    23,
      87,     0,     0,    48,    59,    89,   150,   106,    90,   146,
     148,   147,   149,   152,   151,    91,    92,   143,   132,   133,
     134,   135,   136,   137,   138,   139,   140,   141,   142,   123,
     124,   125,   117,   128,    93,   115,   116,   118,   120,   126,
     127,   122,   107,     0,     0,     0,   117,     0,     0,     0,
       0,     0,     0,     0,     0,   117,     0,    16,     0,    15,
      24,    19,    22,    21,    20,     0,     0,    18,    49,     0,
      52,    54,     0,    53,   107,    60,   109,     0,     0,     0,
       8,   117,   117,    96,    98,    99,   144,   101,   102,   103,
     105,   100,   104,    95,     0,     0,     0,     0,     0,    50,
       0,     0,    61,     0,    64,     0,    65,   110,     5,   121,
     119,     0,   130,    88,    97,    94,     0,     0,     0,     0,
      85,     0,     0,     0,    51,    55,    58,    57,    56,     0,
      62,    66,    67,    70,     0,    84,   111,     0,     0,   129,
       6,   145,    31,     0,    32,    34,    75,    78,    29,     0,
      26,     0,    30,    63,     0,     0,    69,   117,     0,   112,
     131,     9,    33,     0,     0,    77,    25,     0,     0,    68,
      70,     0,    72,    74,   114,   113,     0,    76,    83,    82,
       0,    80,    27,    28,     0,    71,    88,     0,    79,    73,
      10,    81
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -145,  -145,  -145,  -145,    24,    -2,  -145,  -145,  -145,  -145,
    -145,  -128,   146,  -145,   -20,  -145,  -145,   -25,  -145,  -145,
    -145,  -145,    27,   -26,  -145,  -145,    80,  -145,    43,  -145,
    -145,  -145,    48,  -145,  -145,    -3,  -145,  -145,   -18,  -145,
       4,  -145,  -145,   -17,  -145,   -30,   -21,   -40,   -33,   -44,
    -145,  -145,    33,   -99,  -145,  -145,    94,  -145,  -145,  -145,
    -145,  -144,  -145,   -35,   -31,  -100,   -22
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     2,     4,     7,    23,    36,    49,   183,   246,    40,
      63,   126,    64,    65,   131,    66,   219,   220,    67,    68,
      69,   187,   188,    24,    74,   138,   139,   140,   141,   142,
     146,   172,   173,   174,   201,   202,   226,   241,   242,   215,
     216,   235,   250,   251,   204,    25,    26,    27,    28,    29,
     178,   206,   207,   104,   105,   106,   107,   108,   109,   110,
     181,   111,   155,    83,    84,    85,    30
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      35,    75,   103,   147,    37,   182,    72,   166,    20,    71,
      86,    77,    78,    70,   248,    20,    60,   153,   249,    20,
     158,     1,   160,   161,   112,   113,   163,    47,   115,   225,
      72,   208,  -107,    71,   209,    31,     3,    70,   197,   127,
     143,   114,    21,    32,    33,    34,    48,   132,   133,   134,
      37,   144,   184,   185,     5,    38,    79,    80,   103,    39,
      81,    82,    76,   218,   230,   222,    79,    80,     6,    41,
      81,    82,   103,    79,    80,    42,   254,    81,    82,   255,
     112,   103,   257,   156,   157,   258,   154,   179,    43,   148,
     159,    44,    45,   162,   112,   236,    46,   237,    50,    73,
      76,    20,   116,   112,   103,   143,   117,   103,   103,   118,
     119,   120,   176,   127,   128,   175,   144,   121,   122,   123,
     124,   125,   130,   135,   137,   145,   112,   148,   244,   112,
     112,   150,   149,   151,   189,   164,   198,   152,   165,   176,
     169,   167,   175,   190,   171,   170,    20,   -58,   196,   186,
     177,   211,   192,   203,    37,   191,   205,   189,   194,   193,
     200,   214,   223,   256,   224,   229,   190,   217,   227,   127,
     221,   127,    51,   232,    52,    53,    54,    55,    56,    57,
      58,   243,   233,   103,    59,   205,   234,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
     238,   245,   240,   253,    47,   112,   225,   210,   231,   129,
     243,   217,   252,   195,   213,   221,   260,    61,   168,    62,
     199,   239,    21,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,   259,   247,     0,   228,
     261,     0,   180,     0,     0,     0,     0,    22,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    21,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    22,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    21,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,     0,     0,
       0,     0,     0,     0,     0,     0,   136,    20,    87,    88,
      89,    90,    91,    92,    93,    94,    95,    96,    97,    98,
      21,     0,     0,     0,     0,     0,     0,     0,    99,     0,
      21,   100,   101,     0,     0,     0,     0,     0,   212,     0,
       0,     0,     0,     0,     0,     0,     0,   102,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,     0,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    21
};

static const yytype_int16 yycheck[] =
{
      22,    41,    46,   102,    25,   149,    39,   135,    16,    39,
      45,    42,    43,    39,    21,    16,    38,   116,    25,    16,
     120,    31,   122,   123,    46,    47,   125,    38,    50,    61,
      63,    59,    64,    63,    62,    43,    55,    63,     3,    61,
      73,    42,    39,    51,    52,    53,    57,    67,    68,    69,
      71,    73,   151,   152,     0,    64,    21,    22,   102,    33,
      25,    26,    17,   191,   208,   193,    21,    22,    56,    58,
      25,    26,   116,    21,    22,    60,    59,    25,    26,    62,
     102,   125,    59,   118,   119,    62,   117,    57,    60,    59,
     121,    60,    60,   124,   116,    57,    60,    59,    64,    29,
      17,    16,    60,   125,   148,   138,    60,   151,   152,    60,
      60,    60,   145,   135,    61,   145,   138,    60,    60,    60,
      60,    60,    58,    35,    58,    30,   148,    59,   227,   151,
     152,    56,    61,    60,   164,    56,   171,    60,    21,   172,
      58,    63,   172,   164,    60,    59,    16,    60,   170,    59,
      32,   186,    56,   175,   175,    62,   177,   187,    58,    62,
      58,    58,    58,    57,    59,    58,   187,   189,    60,   191,
     192,   193,    42,    58,    44,    45,    46,    47,    48,    49,
      50,   225,    59,   227,    54,   206,    61,     4,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      60,    58,   224,   238,    38,   227,    61,   183,   210,    63,
     254,   233,   237,   170,   187,   237,   256,    34,   138,    36,
     172,   224,    39,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,   254,   233,    -1,   206,
     257,    -1,   148,    -1,    -1,    -1,    -1,    64,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    39,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    64,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    39,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    61,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      39,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    37,    -1,
      39,    40,    41,    -1,    -1,    -1,    -1,    -1,    57,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    56,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    -1,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    39
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    31,    66,    55,    67,     0,    56,    68,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    39,    64,    69,    88,   110,   111,   112,   113,   114,
     131,    43,    51,    52,    53,   131,    70,   111,    64,    33,
      74,    58,    60,    60,    60,    60,    60,    38,    57,    71,
      64,    42,    44,    45,    46,    47,    48,    49,    50,    54,
     131,    34,    36,    75,    77,    78,    80,    83,    84,    85,
      88,   110,   113,    29,    89,   112,    17,   129,   129,    21,
      22,    25,    26,   128,   129,   130,   128,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    28,    37,
      40,    41,    56,   114,   118,   119,   120,   121,   122,   123,
     124,   126,   131,   131,    42,   131,    60,    60,    60,    60,
      60,    60,    60,    60,    60,    60,    76,   131,    61,    77,
      58,    79,    79,    79,    79,    35,    61,    58,    90,    91,
      92,    93,    94,   113,   131,    30,    95,   118,    59,    61,
      56,    60,    60,   118,   129,   127,   128,   128,   130,   129,
     130,   130,   129,   118,    56,    21,    76,    63,    91,    58,
      59,    60,    96,    97,    98,   110,   113,    32,   115,    57,
     121,   125,   126,    72,   118,   118,    59,    86,    87,   110,
     111,    62,    56,    62,    58,    93,   131,     3,   128,    97,
      58,    99,   100,   131,   109,   111,   116,   117,    59,    62,
      69,   128,    57,    87,    58,   104,   105,   131,    76,    81,
      82,   131,    76,    58,    59,    61,   101,    60,   117,    58,
     126,    70,    58,    59,    61,   106,    57,    59,    60,   100,
     131,   102,   103,   114,   118,    58,    73,   105,    21,    25,
     107,   108,    82,   128,    59,    62,    57,    59,    62,   103,
     112,   108
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    65,    66,    67,    68,    69,    70,    70,    72,    73,
      71,    74,    74,    74,    75,    75,    76,    77,    77,    78,
      78,    78,    78,    79,    79,    80,    81,    81,    82,    83,
      84,    85,    86,    86,    87,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    88,    88,    89,    89,    89,
      90,    90,    91,    91,    92,    92,    93,    93,    94,    95,
      95,    95,    96,    96,    97,    97,    98,    99,    99,   100,
     101,   101,   102,   102,   103,   104,   104,   105,   106,   106,
     107,   107,   108,   108,   109,   110,   111,   111,   112,   112,
     113,   113,   113,   113,   113,   113,   113,   113,   113,   113,
     113,   113,   113,   113,   113,   113,   113,   114,   114,   115,
     115,   115,   116,   116,   117,   118,   118,   119,   120,   120,
     121,   121,   122,   122,   122,   122,   122,   122,   123,   124,
     125,   125,   126,   126,   126,   126,   126,   126,   126,   126,
     126,   126,   126,   126,   127,   127,   128,   128,   128,   128,
     129,   130,   130,   131
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     3,     1,     4,     5,     0,     2,     0,     0,
       9,     0,     1,     2,     1,     2,     1,     1,     2,     2,
       2,     2,     2,     0,     1,     6,     1,     3,     3,     5,
       5,     5,     2,     3,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     1,     2,
       2,     3,     1,     1,     1,     3,     3,     3,     1,     0,
       1,     2,     2,     3,     1,     1,     2,     1,     3,     2,
       0,     3,     1,     3,     1,     1,     3,     2,     0,     3,
       1,     3,     1,     1,     1,     1,     1,     1,     0,     3,
       4,     4,     4,     4,     6,     5,     5,     6,     5,     5,
       5,     5,     5,     5,     5,     5,     4,     1,     1,     0,
       1,     2,     2,     3,     3,     1,     1,     0,     1,     3,
       1,     3,     1,     1,     1,     1,     1,     1,     1,     4,
       1,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     3,     1,     1,     1,     1,
       1,     1,     1,     1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 234 "ncgen.y" /* yacc.c:1646  */
    {if (error_count > 0) YYABORT;}
#line 1644 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 3:
#line 237 "ncgen.y" /* yacc.c:1646  */
    {createrootgroup(datasetname);}
#line 1650 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 8:
#line 256 "ncgen.y" /* yacc.c:1646  */
    {
		Symbol* id = (yyvsp[-1].sym);
                markcdf4("Group specification");
		if(creategroup(id) == NULL)
                    yyerror("duplicate group declaration within parent group for %s",
                                id->name);
            }
#line 1662 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 9:
#line 265 "ncgen.y" /* yacc.c:1646  */
    {listpop(groupstack);}
#line 1668 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 12:
#line 271 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1674 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 13:
#line 273 "ncgen.y" /* yacc.c:1646  */
    {markcdf4("Type specification");}
#line 1680 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 16:
#line 279 "ncgen.y" /* yacc.c:1646  */
    { /* Use when defining a type */
              (yyvsp[0].sym)->objectclass = NC_TYPE;
              if(dupobjectcheck(NC_TYPE,(yyvsp[0].sym)))
                    yyerror("duplicate type declaration for %s",
                            (yyvsp[0].sym)->name);
              listpush(typdefs,(void*)(yyvsp[0].sym));
	    }
#line 1692 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 17:
#line 288 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1698 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 18:
#line 288 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1704 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 25:
#line 302 "ncgen.y" /* yacc.c:1646  */
    {
		int i;
                addtogroup((yyvsp[-3].sym)); /* sets prefix*/
                (yyvsp[-3].sym)->objectclass=NC_TYPE;
                (yyvsp[-3].sym)->subclass=NC_ENUM;
                (yyvsp[-3].sym)->typ.basetype=(yyvsp[-5].sym);
                (yyvsp[-3].sym)->typ.size = (yyvsp[-5].sym)->typ.size;
                (yyvsp[-3].sym)->typ.alignment = (yyvsp[-5].sym)->typ.alignment;
                stackbase=(yyvsp[-1].mark);
                stacklen=listlength(stack);
                (yyvsp[-3].sym)->subnodes = listnew();
                /* Variety of field fixups*/
		/* 1. add in the enum values*/
		/* 2. make this type be their container*/
		/* 3. make constant names visible in the group*/
		/* 4. set field basetype to be same as enum basetype*/
                for(i=stackbase;i<stacklen;i++) {
                   Symbol* eid = (Symbol*)listget(stack,i);
		   assert(eid->subclass == NC_ECONST);
		   addtogroup(eid);
                   listpush((yyvsp[-3].sym)->subnodes,(void*)eid);
                   eid->container = (yyvsp[-3].sym);
		   eid->typ.basetype = (yyvsp[-3].sym)->typ.basetype;
                }
                listsetlength(stack,stackbase);/* remove stack nodes*/
              }
#line 1735 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 26:
#line 331 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack); listpush(stack,(void*)(yyvsp[0].sym));}
#line 1741 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 27:
#line 333 "ncgen.y" /* yacc.c:1646  */
    {
		    int i;
		    (yyval.mark)=(yyvsp[-2].mark);
		    /* check for duplicates*/
		    stackbase=(yyvsp[-2].mark);
		    stacklen=listlength(stack);
		    for(i=stackbase;i<stacklen;i++) {
		      Symbol* elem = (Symbol*)listget(stack,i);
		      if(strcmp((yyvsp[0].sym)->name,elem->name)==0)
  	                yyerror("duplicate enum declaration for %s",
        	                 elem->name);
		    }
		    listpush(stack,(void*)(yyvsp[0].sym));
		}
#line 1760 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 28:
#line 350 "ncgen.y" /* yacc.c:1646  */
    {
            (yyvsp[-2].sym)->objectclass=NC_TYPE;
            (yyvsp[-2].sym)->subclass=NC_ECONST;
            (yyvsp[-2].sym)->typ.econst=(yyvsp[0].constant);
	    (yyval.sym)=(yyvsp[-2].sym);
        }
#line 1771 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 29:
#line 359 "ncgen.y" /* yacc.c:1646  */
    {
		    vercheck(NC_OPAQUE);
                    addtogroup((yyvsp[0].sym)); /*sets prefix*/
                    (yyvsp[0].sym)->objectclass=NC_TYPE;
                    (yyvsp[0].sym)->subclass=NC_OPAQUE;
                    (yyvsp[0].sym)->typ.typecode=NC_OPAQUE;
                    (yyvsp[0].sym)->typ.size=int32_val;
                    (yyvsp[0].sym)->typ.alignment=ncaux_class_alignment(NC_OPAQUE);
                }
#line 1785 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 30:
#line 371 "ncgen.y" /* yacc.c:1646  */
    {
                    Symbol* basetype = (yyvsp[-4].sym);
		    vercheck(NC_VLEN);
                    addtogroup((yyvsp[0].sym)); /*sets prefix*/
                    (yyvsp[0].sym)->objectclass=NC_TYPE;
                    (yyvsp[0].sym)->subclass=NC_VLEN;
                    (yyvsp[0].sym)->typ.basetype=basetype;
                    (yyvsp[0].sym)->typ.typecode=NC_VLEN;
                    (yyvsp[0].sym)->typ.size=VLENSIZE;
                    (yyvsp[0].sym)->typ.alignment=ncaux_class_alignment(NC_VLEN);
                }
#line 1801 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 31:
#line 385 "ncgen.y" /* yacc.c:1646  */
    {
	    int i,j;
	    vercheck(NC_COMPOUND);
            addtogroup((yyvsp[-3].sym));
	    /* check for duplicate field names*/
	    stackbase=(yyvsp[-1].mark);
	    stacklen=listlength(stack);
	    for(i=stackbase;i<stacklen;i++) {
	      Symbol* elem1 = (Symbol*)listget(stack,i);
	      for(j=i+1;j<stacklen;j++) {
	          Symbol* elem2 = (Symbol*)listget(stack,j);
	          if(strcmp(elem1->name,elem2->name)==0) {
	            yyerror("duplicate field declaration for %s",elem1->name);
		  }
	      }
	    }
	    (yyvsp[-3].sym)->objectclass=NC_TYPE;
            (yyvsp[-3].sym)->subclass=NC_COMPOUND;
            (yyvsp[-3].sym)->typ.basetype=NULL;
            (yyvsp[-3].sym)->typ.typecode=NC_COMPOUND;
	    (yyvsp[-3].sym)->subnodes = listnew();
	    /* Add in the fields*/
	    for(i=stackbase;i<stacklen;i++) {
	        Symbol* fsym = (Symbol*)listget(stack,i);
		fsym->container = (yyvsp[-3].sym);
 	        listpush((yyvsp[-3].sym)->subnodes,(void*)fsym);
	    }
	    listsetlength(stack,stackbase);/* remove stack nodes*/
          }
#line 1835 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 32:
#line 417 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-1].mark);}
#line 1841 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 33:
#line 418 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-2].mark);}
#line 1847 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 34:
#line 422 "ncgen.y" /* yacc.c:1646  */
    {
	    int i;
	    (yyval.mark)=(yyvsp[0].mark);
	    stackbase=(yyvsp[0].mark);
	    stacklen=listlength(stack);
	    /* process each field in the fieldlist*/
            for(i=stackbase;i<stacklen;i++) {
                Symbol* f = (Symbol*)listget(stack,i);
		f->typ.basetype = (yyvsp[-1].sym);
            }
        }
#line 1863 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 35:
#line 435 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_CHAR]; }
#line 1869 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 36:
#line 436 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_BYTE]; }
#line 1875 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 37:
#line 437 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_SHORT]; }
#line 1881 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 38:
#line 438 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_INT]; }
#line 1887 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 39:
#line 439 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_FLOAT]; }
#line 1893 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 40:
#line 440 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym) = primsymbols[NC_DOUBLE]; }
#line 1899 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 41:
#line 441 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_UBYTE); (yyval.sym) = primsymbols[NC_UBYTE]; }
#line 1905 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 42:
#line 442 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_USHORT); (yyval.sym) = primsymbols[NC_USHORT]; }
#line 1911 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 43:
#line 443 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_UINT); (yyval.sym) = primsymbols[NC_UINT]; }
#line 1917 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 44:
#line 444 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_INT64); (yyval.sym) = primsymbols[NC_INT64]; }
#line 1923 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 45:
#line 445 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_UINT64); (yyval.sym) = primsymbols[NC_UINT64]; }
#line 1929 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 46:
#line 446 "ncgen.y" /* yacc.c:1646  */
    { vercheck(NC_STRING); (yyval.sym) = primsymbols[NC_STRING]; }
#line 1935 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 48:
#line 450 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1941 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 49:
#line 451 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1947 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 52:
#line 458 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1953 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 53:
#line 458 "ncgen.y" /* yacc.c:1646  */
    {}
#line 1959 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 56:
#line 466 "ncgen.y" /* yacc.c:1646  */
    {
		(yyvsp[-2].sym)->dim.declsize = (size_t)extractint((yyvsp[0].constant));
#ifdef GENDEBUG1
fprintf(stderr,"dimension: %s = %llu\n",(yyvsp[-2].sym)->name,(unsigned long long)(yyvsp[-2].sym)->dim.declsize);
#endif
		reclaimconstant((yyvsp[0].constant));
	      }
#line 1971 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 57:
#line 474 "ncgen.y" /* yacc.c:1646  */
    {
		        (yyvsp[-2].sym)->dim.declsize = NC_UNLIMITED;
		        (yyvsp[-2].sym)->dim.isunlimited = 1;
#ifdef GENDEBUG1
fprintf(stderr,"dimension: %s = UNLIMITED\n",(yyvsp[-2].sym)->name);
#endif
		   }
#line 1983 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 58:
#line 484 "ncgen.y" /* yacc.c:1646  */
    {
                     (yyvsp[0].sym)->objectclass=NC_DIM;
                     if(dupobjectcheck(NC_DIM,(yyvsp[0].sym)))
                        yyerror( "Duplicate dimension declaration for %s",
                                (yyvsp[0].sym)->name);
		     addtogroup((yyvsp[0].sym));
		     (yyval.sym)=(yyvsp[0].sym);
		     listpush(dimdefs,(void*)(yyvsp[0].sym));
                   }
#line 1997 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 60:
#line 496 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2003 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 61:
#line 497 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2009 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 64:
#line 504 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2015 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 65:
#line 504 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2021 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 66:
#line 507 "ncgen.y" /* yacc.c:1646  */
    {
		    int i;
		    stackbase=(yyvsp[0].mark);
		    stacklen=listlength(stack);
		    /* process each variable in the varlist*/
	            for(i=stackbase;i<stacklen;i++) {
	                Symbol* sym = (Symbol*)listget(stack,i);
			sym->objectclass = NC_VAR;
		        if(dupobjectcheck(NC_VAR,sym)) {
                            yyerror("Duplicate variable declaration for %s",
                                    sym->name);
			} else {
		  	    sym->typ.basetype = (yyvsp[-1].sym);
	                    addtogroup(sym);
		            listpush(vardefs,(void*)sym);
			    sym->var.special = ecalloc(sizeof(Specialdata));
			    if(sym->var.special == NULL)
			        derror("out of memory");
			}
		    }
		    listsetlength(stack,stackbase);/* remove stack nodes*/
		}
#line 2048 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 67:
#line 532 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack);
                 listpush(stack,(void*)(yyvsp[0].sym));
		}
#line 2056 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 68:
#line 536 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-2].mark); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2062 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 69:
#line 540 "ncgen.y" /* yacc.c:1646  */
    {
		    int i;
		    Dimset dimset;
		    Symbol* var = (yyvsp[-1].sym); /* for debugging */
		    stacklen=listlength(stack);
		    stackbase=(yyvsp[0].mark);
		    count = stacklen - stackbase;
		    if(count >= NC_MAX_VAR_DIMS) {
			yyerror("%s has too many dimensions",(yyvsp[-1].sym)->name);
			count = NC_MAX_VAR_DIMS - 1;
			stacklen = stackbase + count;
		    }
  	            dimset.ndims = count;
		    /* extract the actual dimensions*/
		    if(dimset.ndims > 0) {
		        for(i=0;i<count;i++) {
			    Symbol* dsym = (Symbol*)listget(stack,stackbase+i);
			    dimset.dimsyms[i] = dsym;
			}
			var->typ.dimset = dimset;
		    }
		    var->typ.basetype = NULL; /* not yet known*/
                    var->objectclass=NC_VAR;
		    listsetlength(stack,stackbase);/* remove stack nodes*/
		    (yyval.sym) = var;
		    }
#line 2093 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 70:
#line 568 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack);}
#line 2099 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 71:
#line 569 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-1].mark);}
#line 2105 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 72:
#line 572 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2111 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 73:
#line 574 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-2].mark); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2117 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 74:
#line 578 "ncgen.y" /* yacc.c:1646  */
    {Symbol* dimsym = (yyvsp[0].sym);
		dimsym->objectclass = NC_DIM;
		/* Find the actual dimension*/
		dimsym = locate(dimsym);
		if(dimsym == NULL) {
		    derror("Undefined or forward referenced dimension: %s",(yyvsp[0].sym)->name);
		    YYABORT;
		}
		(yyval.sym)=dimsym;
	    }
#line 2132 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 75:
#line 592 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack);
             listpush(stack,(void*)(yyvsp[0].sym));
	    }
#line 2140 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 76:
#line 596 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-2].mark); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2146 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 77:
#line 601 "ncgen.y" /* yacc.c:1646  */
    {
		int i;
		Dimset dimset;
		stackbase=(yyvsp[0].mark);
		stacklen=listlength(stack);
		count = stacklen - stackbase;
		if(count >= NC_MAX_VAR_DIMS) {
		    yyerror("%s has too many dimensions",(yyvsp[-1].sym)->name);
		    count = NC_MAX_VAR_DIMS - 1;
		    stacklen = stackbase + count;
		}
  	        dimset.ndims = count;
		if(count > 0) {
		    /* extract the actual dimensions*/
		    for(i=0;i<count;i++) {
		        Symbol* dsym = (Symbol*)listget(stack,stackbase+i);
		        dimset.dimsyms[i] = dsym;
		    }
		    (yyvsp[-1].sym)->typ.dimset = dimset;
		}
		(yyvsp[-1].sym)->typ.basetype = NULL; /* not yet known*/
                (yyvsp[-1].sym)->objectclass=NC_TYPE;
                (yyvsp[-1].sym)->subclass=NC_FIELD;
		listsetlength(stack,stackbase);/* remove stack nodes*/
		(yyval.sym) = (yyvsp[-1].sym);
	    }
#line 2177 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 78:
#line 629 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack);}
#line 2183 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 79:
#line 630 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-1].mark);}
#line 2189 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 80:
#line 634 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=listlength(stack); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2195 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 81:
#line 636 "ncgen.y" /* yacc.c:1646  */
    {(yyval.mark)=(yyvsp[-2].mark); listpush(stack,(void*)(yyvsp[0].sym));}
#line 2201 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 82:
#line 641 "ncgen.y" /* yacc.c:1646  */
    {  /* Anonymous integer dimension.
	         Can only occur in type definitions*/
	     char anon[32];
	     sprintf(anon,"const%u",uint32_val);
	     (yyval.sym) = install(anon);
	     (yyval.sym)->objectclass = NC_DIM;
	     (yyval.sym)->dim.isconstant = 1;
	     (yyval.sym)->dim.declsize = uint32_val;
	    }
#line 2215 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 83:
#line 651 "ncgen.y" /* yacc.c:1646  */
    {  /* Anonymous integer dimension.
	         Can only occur in type definitions*/
	     char anon[32];
	     if(int32_val <= 0) {
		derror("field dimension must be positive");
		YYABORT;
	     }
	     sprintf(anon,"const%d",int32_val);
	     (yyval.sym) = install(anon);
	     (yyval.sym)->objectclass = NC_DIM;
	     (yyval.sym)->dim.isconstant = 1;
	     (yyval.sym)->dim.declsize = int32_val;
	    }
#line 2233 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 84:
#line 671 "ncgen.y" /* yacc.c:1646  */
    {Symbol* vsym = (yyvsp[0].sym);
		if(vsym->objectclass != NC_VAR) {
		    derror("Undefined or forward referenced variable: %s",vsym->name);
		    YYABORT;
		}
		(yyval.sym)=vsym;
	    }
#line 2245 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 85:
#line 682 "ncgen.y" /* yacc.c:1646  */
    {Symbol* tsym = (yyvsp[0].sym);
		if(tsym->objectclass != NC_TYPE) {
		    derror("Undefined or forward referenced type: %s",tsym->name);
		    YYABORT;
		}
		(yyval.sym)=tsym;
	    }
#line 2257 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 86:
#line 693 "ncgen.y" /* yacc.c:1646  */
    {Symbol* tvsym = (yyvsp[0].sym); Symbol* sym;
		/* disambiguate*/
		tvsym->objectclass = NC_VAR;
		sym = locate(tvsym);
		if(sym == NULL) {
		    tvsym->objectclass = NC_TYPE;
		    sym = locate(tvsym);
		    if(tvsym == NULL) {
		        derror("Undefined or forward referenced name: %s",(yyvsp[0].sym)->name);
		        YYABORT;
		    } else tvsym = sym;
		} else tvsym = sym;
		if(tvsym == NULL) {
		    derror("Undefined name (line %d): %s",(yyvsp[0].sym)->lineno,(yyvsp[0].sym)->name);
		    YYABORT;
		}
		(yyval.sym)=tvsym;
	    }
#line 2280 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 87:
#line 711 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym)=(yyvsp[0].sym);}
#line 2286 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 88:
#line 718 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2292 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 89:
#line 718 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2298 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 90:
#line 722 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_NCPROPS_FLAG,NULL,NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2304 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 91:
#line 724 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_ISNETCDF4_FLAG,NULL,NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2310 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 92:
#line 726 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_SUPERBLOCK_FLAG,NULL,NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2316 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 93:
#line 728 "ncgen.y" /* yacc.c:1646  */
    { (yyval.sym)=makeattribute((yyvsp[-2].sym),NULL,NULL,(yyvsp[0].datalist),ATTRGLOBAL);}
#line 2322 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 94:
#line 730 "ncgen.y" /* yacc.c:1646  */
    {Symbol* tsym = (yyvsp[-5].sym); Symbol* vsym = (yyvsp[-4].sym); Symbol* asym = (yyvsp[-2].sym);
		if(vsym->objectclass == NC_VAR) {
		    (yyval.sym)=makeattribute(asym,vsym,tsym,(yyvsp[0].datalist),ATTRVAR);
		} else {
		    derror("Doubly typed attribute: %s",asym->name);
		    YYABORT;
		}
	    }
#line 2335 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 95:
#line 739 "ncgen.y" /* yacc.c:1646  */
    {Symbol* sym = (yyvsp[-4].sym); Symbol* asym = (yyvsp[-2].sym);
		if(sym->objectclass == NC_VAR) {
		    (yyval.sym)=makeattribute(asym,sym,NULL,(yyvsp[0].datalist),ATTRVAR);
		} else if(sym->objectclass == NC_TYPE) {
		    (yyval.sym)=makeattribute(asym,NULL,sym,(yyvsp[0].datalist),ATTRGLOBAL);
		} else {
		    derror("Attribute prefix not a variable or type: %s",asym->name);
		    YYABORT;
		}
	    }
#line 2350 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 96:
#line 750 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_FILLVALUE_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].datalist),ISLIST);}
#line 2356 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 97:
#line 752 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_FILLVALUE_FLAG,(yyvsp[-4].sym),(yyvsp[-5].sym),(void*)(yyvsp[0].datalist),ISLIST);}
#line 2362 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 98:
#line 754 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_STORAGE_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2368 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 99:
#line 756 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_CHUNKSIZES_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].datalist),ISLIST);}
#line 2374 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 100:
#line 758 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_FLETCHER32_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2380 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 101:
#line 760 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_DEFLATE_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2386 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 102:
#line 762 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_SHUFFLE_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2392 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 103:
#line 764 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_ENDIAN_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2398 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 104:
#line 766 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_FILTER_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2404 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 105:
#line 768 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_NOFILL_FLAG,(yyvsp[-4].sym),NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2410 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 106:
#line 770 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym) = makespecial(_FORMAT_FLAG,NULL,NULL,(void*)(yyvsp[0].constant),ISCONST);}
#line 2416 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 107:
#line 775 "ncgen.y" /* yacc.c:1646  */
    {
	        (yyval.sym)=(yyvsp[0].sym);
                (yyvsp[0].sym)->ref.is_ref=1;
                (yyvsp[0].sym)->is_prefixed=0;
                setpathcurrent((yyvsp[0].sym));
	    }
#line 2427 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 108:
#line 782 "ncgen.y" /* yacc.c:1646  */
    {
	        (yyval.sym)=(yyvsp[0].sym);
                (yyvsp[0].sym)->ref.is_ref=1;
                (yyvsp[0].sym)->is_prefixed=1;
	        /* path is set in ncgen.l*/
	    }
#line 2438 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 110:
#line 791 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2444 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 111:
#line 792 "ncgen.y" /* yacc.c:1646  */
    {}
#line 2450 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 114:
#line 800 "ncgen.y" /* yacc.c:1646  */
    {(yyvsp[-2].sym)->data = (yyvsp[0].datalist);}
#line 2456 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 115:
#line 803 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = (yyvsp[0].datalist);}
#line 2462 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 116:
#line 804 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = (yyvsp[0].datalist);}
#line 2468 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 117:
#line 808 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = builddatalist(0);}
#line 2474 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 118:
#line 812 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = const2list((yyvsp[0].constant));}
#line 2480 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 119:
#line 814 "ncgen.y" /* yacc.c:1646  */
    {dlappend((yyvsp[-2].datalist),((yyvsp[0].constant))); (yyval.datalist)=(yyvsp[-2].datalist); }
#line 2486 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 120:
#line 818 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=(yyvsp[0].constant);}
#line 2492 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 121:
#line 819 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=builddatasublist((yyvsp[-1].datalist));}
#line 2498 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 122:
#line 823 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=(yyvsp[0].constant);}
#line 2504 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 123:
#line 824 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_OPAQUE);}
#line 2510 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 124:
#line 825 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_FILLVALUE);}
#line 2516 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 125:
#line 826 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_NIL);}
#line 2522 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 126:
#line 827 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=(yyvsp[0].constant);}
#line 2528 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 128:
#line 832 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant) = makeenumconstref((yyvsp[0].sym));}
#line 2534 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 129:
#line 836 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=evaluate((yyvsp[-3].sym),(yyvsp[-1].datalist));}
#line 2540 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 130:
#line 841 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = const2list((yyvsp[0].constant));}
#line 2546 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 131:
#line 843 "ncgen.y" /* yacc.c:1646  */
    {dlappend((yyvsp[-2].datalist),((yyvsp[0].constant))); (yyval.datalist)=(yyvsp[-2].datalist);}
#line 2552 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 132:
#line 847 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_CHAR);}
#line 2558 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 133:
#line 848 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_BYTE);}
#line 2564 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 134:
#line 849 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_SHORT);}
#line 2570 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 135:
#line 850 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_INT);}
#line 2576 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 136:
#line 851 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_INT64);}
#line 2582 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 137:
#line 852 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_UBYTE);}
#line 2588 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 138:
#line 853 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_USHORT);}
#line 2594 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 139:
#line 854 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_UINT);}
#line 2600 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 140:
#line 855 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_UINT64);}
#line 2606 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 141:
#line 856 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_FLOAT);}
#line 2612 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 142:
#line 857 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_DOUBLE);}
#line 2618 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 143:
#line 858 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_STRING);}
#line 2624 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 144:
#line 862 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist) = const2list((yyvsp[0].constant));}
#line 2630 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 145:
#line 863 "ncgen.y" /* yacc.c:1646  */
    {(yyval.datalist)=(yyvsp[-2].datalist); dlappend((yyvsp[-2].datalist),((yyvsp[0].constant)));}
#line 2636 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 146:
#line 868 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_INT);}
#line 2642 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 147:
#line 870 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_UINT);}
#line 2648 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 148:
#line 872 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_INT64);}
#line 2654 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 149:
#line 874 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_UINT64);}
#line 2660 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 150:
#line 878 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=makeconstdata(NC_STRING);}
#line 2666 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 151:
#line 882 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=(yyvsp[0].constant);}
#line 2672 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 152:
#line 883 "ncgen.y" /* yacc.c:1646  */
    {(yyval.constant)=(yyvsp[0].constant);}
#line 2678 "ncgeny.c" /* yacc.c:1646  */
    break;

  case 153:
#line 889 "ncgen.y" /* yacc.c:1646  */
    {(yyval.sym)=(yyvsp[0].sym);}
#line 2684 "ncgeny.c" /* yacc.c:1646  */
    break;


#line 2688 "ncgeny.c" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 892 "ncgen.y" /* yacc.c:1906  */


#ifndef NO_STDARG
static void
yyerror(const char *fmt, ...)
#else
static void
yyerror(fmt,va_alist) const char* fmt; va_dcl
#endif
{
    va_list argv;
    va_start(argv,fmt);
    (void)fprintf(stderr,"%s: %s line %d: ", progname, cdlname, lineno);
    vderror(fmt,argv);
    va_end(argv);
}

/* undefine yywrap macro, in case we are using bison instead of yacc */
#ifdef yywrap
#undef yywrap
#endif

static int
ncgwrap(void)                    /* returns 1 on EOF if no more input */
{
    return  1;
}

/* get lexical input routine generated by lex  */
#include "ncgenl.c"

/* Really should init our data within this file */
void
parse_init(void)
{
    int i;
    opaqueid = 0;
    arrayuid = 0;
    symlist = listnew();
    stack = listnew();
    groupstack = listnew();
    consttype = NC_NAT;
    grpdefs = listnew();
    dimdefs = listnew();
    attdefs = listnew();
    gattdefs = listnew();
    xattdefs = listnew();
    typdefs = listnew();
    vardefs = listnew();
    condefs = listnew();
    tmp = listnew();
    /* Create the primitive types */
    for(i=NC_NAT+1;i<=NC_STRING;i++) {
        primsymbols[i] = makeprimitivetype(i);
    }
    lex_init();
}

static Symbol*
makeprimitivetype(nc_type nctype)
{
    Symbol* sym = install(primtypenames[nctype]);
    sym->objectclass=NC_TYPE;
    sym->subclass=NC_PRIM;
    sym->nc_id = nctype;
    sym->typ.typecode = nctype;
    sym->typ.size = ncsize(nctype);
    sym->typ.nelems = 1;
    sym->typ.alignment = ncaux_class_alignment(nctype); 
    /* Make the basetype circular so we can always ask for it */
    sym->typ.basetype = sym;
    sym->prefix = listnew();
    return sym;
}

/* Symbol table operations for ncgen tool */
/* install sname in symbol table even if it is already there */
Symbol*
install(const char *sname)
{
    Symbol* sp;
    sp = (Symbol*) ecalloc (sizeof (struct Symbol));
    sp->name = nulldup(sname);
    sp->lineno = lineno;
    sp->location = currentgroup();
    sp->container = currentgroup();
    listpush(symlist,sp);
    return sp;
}

static Symbol*
currentgroup(void)
{
    if(listlength(groupstack) == 0) return rootgroup;
    return (Symbol*)listtop(groupstack);
}

static Symbol*
createrootgroup(const char* dataset)
{
    Symbol* gsym = install(dataset);
    gsym->objectclass = NC_GRP;
    gsym->container = NULL;
    gsym->subnodes = listnew();
    gsym->grp.is_root = 1;
    gsym->prefix = listnew();
    listpush(grpdefs,(void*)gsym);
    rootgroup = gsym;
    return gsym;
}

static Symbol*
creategroup(Symbol * gsym)
{
    /* See if this group already exists in currentgroup */
    gsym->objectclass = NC_GRP;
    if(dupobjectcheck(NC_GRP,gsym)) {
        derror("Duplicate group name in same scope: %s",gsym->name);
	return NULL;
    }
    addtogroup(gsym);
    gsym->subnodes = listnew();
    listpush(groupstack,(void*)gsym);
    listpush(grpdefs,(void*)gsym);
    return gsym;
}

static NCConstant*
makeconstdata(nc_type nctype)
{
    NCConstant* con = nullconst();
    consttype = nctype;
    con->nctype = nctype;
    con->lineno = lineno;
    con->filled = 0;
    switch (nctype) {
	case NC_CHAR: con->value.charv = char_val; break;
        case NC_BYTE: con->value.int8v = byte_val; break;
        case NC_SHORT: con->value.int16v = int16_val; break;
        case NC_INT: con->value.int32v = int32_val; break;
        case NC_FLOAT:
	    con->value.floatv = float_val;
	    break;
        case NC_DOUBLE:
	    con->value.doublev = double_val;
	    break;
        case NC_STRING: { /* convert to a set of chars*/
	    size_t len;
	    len = bbLength(lextext);
	    con->value.stringv.len = len;
	    con->value.stringv.stringv = bbExtract(lextext);
	    }
	    break;

	/* Allow these constants even in netcdf-3 */
        case NC_UBYTE: con->value.uint8v = ubyte_val; break;
        case NC_USHORT: con->value.uint16v = uint16_val; break;
        case NC_UINT: con->value.uint32v = uint32_val; break;
        case NC_INT64: con->value.int64v = int64_val; break;
        case NC_UINT64: con->value.uint64v = uint64_val; break;

#ifdef USE_NETCDF4
	case NC_OPAQUE: {
	    char* s;
	    int len;
	    len = bbLength(lextext);
	    s = (char*)ecalloc(len+1);
	    strncpy(s,bbContents(lextext),len);
	    s[len] = '\0';
	    con->value.opaquev.stringv = s;
	    con->value.opaquev.len = len;
	    } break;

	case NC_NIL:
	    break; /* no associated value*/
#endif

 	case NC_FILLVALUE:
	    break; /* no associated value*/

	default:
	    yyerror("Data constant: unexpected NC type: %s",
		    nctypename(nctype));
	    con->value.stringv.stringv = NULL;
	    con->value.stringv.len = 0;
    }
    return con;
}

static NCConstant*
makeenumconstref(Symbol* refsym)
{
    NCConstant* con = nullconst();

    markcdf4("Enum type");
    consttype = NC_ENUM;
    con->nctype = NC_ECONST;
    con->lineno = lineno;
    con->filled = 0;
    refsym->objectclass = NC_TYPE;
    refsym->subclass = NC_ECONST;
    con->value.enumv = refsym;
    return con;
}

static void
addtogroup(Symbol* sym)
{
    Symbol* grp = currentgroup();
    sym->container = grp;
    listpush(grp->subnodes,(void*)sym);
    setpathcurrent(sym);
}

/* Check for duplicate name of given type within current group*/
static int
dupobjectcheck(nc_class objectclass, Symbol* pattern)
{
    int i;
    Symbol* grp;
    if(pattern == NULL) return 0;
    grp = pattern->container;
    if(grp == NULL || grp->subnodes == NULL) return 0;
    for(i=0;i<listlength(grp->subnodes);i++) {
	Symbol* sym = (Symbol*)listget(grp->subnodes,i);
	if(!sym->ref.is_ref && sym->objectclass == objectclass
	   && strcmp(sym->name,pattern->name)==0) return 1;
    }
    return 0;
}

static void
setpathcurrent(Symbol* sym)
{
    sym->is_prefixed = 0;
    sym->prefix = prefixdup(groupstack);
}

/* Convert an nc_type code to the corresponding Symbol*/
Symbol*
basetypefor(nc_type nctype)
{
    return primsymbols[nctype];
}

static int
truefalse(NCConstant* con, int tag)
{
    if(con->nctype == NC_STRING) {
	char* sdata = con->value.stringv.stringv;
	if(strncmp(sdata,"false",NC_MAX_NAME) == 0
           || strncmp(sdata,"0",NC_MAX_NAME) == 0)
	    return 0;
	else if(strncmp(sdata,"true",NC_MAX_NAME) == 0
           || strncmp(sdata,"1",NC_MAX_NAME) == 0)
	    return 1;
	else goto fail;
    } else if(con->value.int32v < 0 || con->value.int32v > 1)
	goto fail;
    return con->value.int32v;

fail:
    derror("%s: illegal value",specialname(tag));
    return 0;
}

/* Since this may be affected by the _Format attribute, which
   may come last, capture all the special info and sort it out
   in semantics.
*/
static Symbol*
makespecial(int tag, Symbol* vsym, Symbol* tsym, void* data, int isconst)
{
    Symbol* attr = NULL;
    Datalist* list = NULL;
    NCConstant* con = NULL;
    NCConstant* tmp = NULL;
    int tf = 0;
    char* sdata = NULL;
    int idata =  -1;

    if((GLOBAL_SPECIAL & tag) != 0) {
        if(vsym != NULL) {
            derror("_Format: must be global attribute");
            vsym = NULL;
        }
    } else {
        if(vsym == NULL) {
	    derror("%s: must have non-NULL vsym", specialname(tag));
	    return NULL;
        }
    }

    if(tag != _FILLVALUE_FLAG && tag != _FORMAT_FLAG)
        /*Main.*/specials_flag++;

    if(isconst)
	con = (NCConstant*)data;
    else
        list = (Datalist*)data;

    switch (tag) {
    case _FLETCHER32_FLAG:
    case _SHUFFLE_FLAG:
    case _ISNETCDF4_FLAG:
    case _NOFILL_FLAG:
	tmp = nullconst();
	tmp->nctype = (con->nctype == NC_STRING?NC_STRING:NC_INT);
	convert1(con,tmp);
	tf = truefalse(tmp,tag);
	reclaimconstant(tmp);
	break;
    case _FORMAT_FLAG:
    case _STORAGE_FLAG:
    case _NCPROPS_FLAG:
    case _ENDIAN_FLAG:
    case _FILTER_FLAG:
	tmp = nullconst();
        tmp->nctype = NC_STRING;
	convert1(con,tmp);
	if(tmp->nctype == NC_STRING) {
	    sdata = tmp->value.stringv.stringv;
	    tmp->value.stringv.stringv = NULL;
	    tmp->value.stringv.len = 0;
	} else
	    derror("%s: illegal value",specialname(tag));
	reclaimconstant(tmp);
	break;
    case _SUPERBLOCK_FLAG:
    case _DEFLATE_FLAG:
	tmp = nullconst();
        tmp->nctype = NC_INT;
	convert1(con,tmp);
	if(tmp->nctype == NC_INT)
	    idata = tmp->value.int32v;
	else
	    derror("%s: illegal value",specialname(tag));
	reclaimconstant(tmp);
	break;
    case _CHUNKSIZES_FLAG:
    case _FILLVALUE_FLAG:
	/* Handle below */
	break;
    default: PANIC1("unexpected special tag: %d",tag);
    }

    if(tag == _FORMAT_FLAG) {
	/* Watch out: this is a global attribute */
	struct Kvalues* kvalue;
	int found = 0;
	/* Use the table in main.c */
        for(kvalue = legalkinds; kvalue->name; kvalue++) {
          if(sdata) {
            if(strcmp(sdata, kvalue->name) == 0) {
              globalspecials._Format = kvalue->k_flag;
	      /*Main.*/format_attribute = 1;
              found = 1;
              break;
            }
          }
	}
	if(!found)
	    derror("_Format: illegal value: %s",sdata);
    } else if((GLOBAL_SPECIAL & tag) != 0) {
	if(tag == _ISNETCDF4_FLAG)
	    globalspecials._IsNetcdf4 = tf;
	else if(tag == _SUPERBLOCK_FLAG)
	    globalspecials._Superblock = idata;
	else if(tag == _NCPROPS_FLAG) {
	    globalspecials._NCProperties = sdata;
	    sdata = NULL;
	}	    
    } else {
        Specialdata* special;
        /* Set up special info */
	if(vsym->var.special == NULL) {
            vsym->var.special = ecalloc(sizeof(Specialdata));
	    if(vsym->var.special == NULL)
	        derror("Out of memory");
	}
        special = vsym->var.special;
        if(tag == _FILLVALUE_FLAG) {
            /* fillvalue must be a single value*/
	    if(!isconst && datalistlen(list) != 1)
                derror("_FillValue: must be a single (possibly compound) value",
                            vsym->name);
	    if(isconst) {
	        list = const2list(con);
		con = NULL;
	    }
            /* check that the attribute value contains no fill values*/
            if(containsfills(list)) {
                derror("Attribute data may not contain fill values (i.e. _ )");
            }
            /* _FillValue is also a real attribute*/
            if(vsym->objectclass != NC_VAR) {
                derror("_FillValue attribute not associated with variable: %s",vsym->name);
            }
            if(tsym  == NULL) tsym = vsym->typ.basetype;
            else if(vsym->typ.basetype != tsym) {
                derror("_FillValue attribute type does not match variable type: %s",vsym->name);
            }
            special->_Fillvalue = clonedatalist(list);
	    /* Create the corresponding attribute */
            attr = makeattribute(install("_FillValue"),vsym,tsym,list,ATTRVAR);
	    list = NULL;
        } else switch (tag) {
	    /* These will be output as attributes later */
            case _STORAGE_FLAG:
              if(!sdata)
                derror("_Storage: illegal NULL value");
              else if(strcmp(sdata,"contiguous") == 0)
                special->_Storage = NC_CONTIGUOUS;
              else if(strcmp(sdata,"chunked") == 0)
                special->_Storage = NC_CHUNKED;
              else
                derror("_Storage: illegal value: %s",sdata);
              special->flags |= _STORAGE_FLAG;
              break;
          case _FLETCHER32_FLAG:
                special->_Fletcher32 = tf;
                special->flags |= _FLETCHER32_FLAG;
                break;
            case _DEFLATE_FLAG:
                special->_DeflateLevel = idata;
                special->flags |= _DEFLATE_FLAG;
                break;
            case _SHUFFLE_FLAG:
                special->_Shuffle = tf;
                special->flags |= _SHUFFLE_FLAG;
                break;
            case _ENDIAN_FLAG:
              if(!sdata)
                derror("_Endianness: illegal NULL value");
              else if(strcmp(sdata,"little") == 0)
                special->_Endianness = 1;
              else if(strcmp(sdata,"big") == 0)
                special->_Endianness = 2;
              else
                derror("_Endianness: illegal value: %s",sdata);
              special->flags |= _ENDIAN_FLAG;
              break;
          case _NOFILL_FLAG:
                special->_Fill = (1 - tf); /* negate */
                special->flags |= _NOFILL_FLAG;
                break;
          case _CHUNKSIZES_FLAG: {
                int i;
		list = (isconst ? const2list(con) : list);		
                special->nchunks = list->length;
                special->_ChunkSizes = (size_t*)ecalloc(sizeof(size_t)*special->nchunks);
                for(i=0;i<special->nchunks;i++) {
		    tmp = nullconst();
                    tmp->nctype = NC_INT;
                    convert1(list->data[i],tmp);
                    if(tmp->nctype == NC_INT) {
                        special->_ChunkSizes[i] = (size_t)tmp->value.int32v;
                    } else {
                        efree(special->_ChunkSizes);
                        derror("%s: illegal value",specialname(tag));
                    }
		    reclaimconstant(tmp);
                }
                special->flags |= _CHUNKSIZES_FLAG;
                /* Chunksizes => storage == chunked */
                special->flags |= _STORAGE_FLAG;
                special->_Storage = NC_CHUNKED;
                } break;
          case _FILTER_FLAG:
#ifdef USE_NETCDF4
		/* Parse the filter spec */
		if(parsefilterflag(sdata,special) == NC_NOERR)
                    special->flags |= _FILTER_FLAG;
		else {
		    efree(special->_FilterParams);
		    derror("_Filter: unparseable filter spec: %s",sdata);
		}
#else
	        derror("%s: the filter attribute requires netcdf-4 to be enabled",specialname(tag));
#endif
                break;
            default: PANIC1("makespecial: illegal token: %d",tag);
         }
    }
    if(sdata) efree(sdata);
    if(con) reclaimconstant(con);
    if(list) reclaimdatalist(list);
    return attr;
}

static Symbol*
makeattribute(Symbol* asym,
		Symbol* vsym,
		Symbol* tsym,
		Datalist* data,
		Attrkind kind) /* global var or unknown*/
{
    asym->objectclass = NC_ATT;
    asym->data = data;
    switch (kind) {
    case ATTRVAR:
        asym->att.var = vsym;
        asym->typ.basetype = tsym;
        listpush(attdefs,(void*)asym);
        addtogroup(asym);
	break;
    case ATTRGLOBAL:
        asym->att.var = NULL; /* NULL => NC_GLOBAL*/
        asym->typ.basetype = tsym;
        listpush(gattdefs,(void*)asym);
        addtogroup(asym);
	break;
    default: PANIC1("unexpected attribute type: %d",kind);
    }
    /* finally; check that the attribute value contains no fill values*/
    if(containsfills(data)) {
	derror("Attribute data may not contain fill values (i.e. _ ): %s",asym->name);
    }
    return asym;
}

static long long
extractint(NCConstant* con)
{
    switch (con->nctype) {
    case NC_BYTE: return (long long)(con->value.int8v);
    case NC_SHORT: return (long long)(con->value.int16v);
    case NC_INT: return (long long)(con->value.int32v);
    case NC_UBYTE: return (long long)(con->value.uint8v);
    case NC_USHORT: return (long long)(con->value.uint16v);
    case NC_UINT: return (long long)(con->value.uint32v);
    case NC_INT64: return (long long)(con->value.int64v);
    default:
	derror("Not a signed integer type: %d",con->nctype);
	break;
    }
    return 0;
}

static int
containsfills(Datalist* list)
{
    if(list != NULL) {
        int i;
        NCConstant** cons = list->data;
        for(i=0;i<list->length;i++) {
	    if(cons[i]->nctype == NC_COMPOUND) {
	        if(containsfills(cons[i]->value.compoundv)) return 1;
	    } else if(cons[i]->nctype == NC_FILLVALUE)
		return 1;
	}
    }
    return 0;
}

/*
Try to infer the file type from the
kinds of constructs used in the cdl file.
*/
static void
vercheck(int tid)
{
    switch (tid) {
    case NC_UBYTE: markcdf4("netCDF4/5 type: UBYTE"); break;
    case NC_USHORT: markcdf4("netCDF4/5 type: USHORT"); break;
    case NC_UINT: markcdf4("netCDF4/5 type: UINT"); break;
    case NC_INT64: markcdf4("netCDF4/5 type: INT64"); break;
    case NC_UINT64: markcdf4("netCDF4/5 type: UINT64"); break;
    case NC_STRING: markcdf4("netCDF4 type: STRING"); break;
    case NC_VLEN: markcdf4("netCDF4 type: VLEN"); break;
    case NC_OPAQUE: markcdf4("netCDF4 type: OPAQUE"); break;
    case NC_ENUM: markcdf4("netCDF4 type: ENUM"); break;
    case NC_COMPOUND: markcdf4("netCDF4 type: COMPOUND"); break;
    default: break;
    }
}

const char*
specialname(int tag)
{
    struct Specialtoken* spp = specials;
    for(;spp->name;spp++) {
	if(spp->tag == tag)
	    return spp->name;
    }
    return "<unknown>";
}

#ifdef USE_NETCDF4
/*
Parse a filter spec string and store it in special
*/
static int
parsefilterflag(const char* sdata, Specialdata* special)
{
    int stat = NC_NOERR;

    if(sdata == NULL || strlen(sdata) == 0) return NC_EINVAL;

    stat = NC_parsefilterspec(sdata, &special->_FilterID, &special->nparams, &special->_FilterParams);
    if(stat)
        derror("Malformed filter spec: %s",sdata);
    return stat;
}
#endif

/*
Since the arguments are all simple constants,
we can evaluate the function immediately
and return its value.
Note that currently, only a single value can
be returned.
*/

static NCConstant*
evaluate(Symbol* fcn, Datalist* arglist)
{
    NCConstant* result = nullconst();

    /* prepare the result */
    result->lineno = fcn->lineno;

    if(strcasecmp(fcn->name,"time") == 0) {
        char* timekind = NULL;
        char* timevalue = NULL;
        result->nctype = NC_DOUBLE;
        result->value.doublev = 0;
	/* int time([string],string) */
	switch (arglist->length) {
	case 2:
	    if(arglist->data[1]->nctype != NC_STRING) {
	        derror("Expected function signature: time([string,]string)");
	        goto done;
	    }
	    /* fall thru */
	case 1:
	    if(arglist->data[0]->nctype != NC_STRING) {
	        derror("Expected function signature: time([string,]string)");
	        goto done;
	    }
	    break;
	case 0:
	default:
	    derror("Expected function signature: time([string,]string)");
	    goto done;
	}
	if(arglist->length == 2) {
	    timekind = arglist->data[0]->value.stringv.stringv;
            timevalue = arglist->data[1]->value.stringv.stringv;
	} else
            timevalue = arglist->data[0]->value.stringv.stringv;
	if(timekind == NULL) { /* use cd time as the default */
            cdCompTime comptime;
	    CdTime cdtime;
	    cdCalenType timetype = cdStandard;
	    cdChar2Comp(timetype,timevalue,&comptime);
	    /* convert comptime to cdTime */
	    cdtime.year = comptime.year;
	    cdtime.month = comptime.month;
	    cdtime.day = comptime.day;
	    cdtime.hour = comptime.hour;
	    cdtime.baseYear = 1970;
	    cdtime.timeType = CdChron;
	    /* convert to double value */
	    Cdh2e(&cdtime,&result->value.doublev);
        } else {
	    derror("Time conversion '%s' not supported",timekind);
	    goto done;
	}
    } else {	/* Unknown function */
	derror("Unknown function name: %s",fcn->name);
	goto done;
    }

done:
    return result;
}
