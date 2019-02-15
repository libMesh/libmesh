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
#define YYPURE 1

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1


/* Substitute the variable and function names.  */
#define yyparse         dapparse
#define yylex           daplex
#define yyerror         daperror
#define yydebug         dapdebug
#define yynerrs         dapnerrs


/* Copy the first part of user declarations.  */
#line 11 "dap.y" /* yacc.c:339  */

#include "config.h"
#include "dapparselex.h"
#include "dapy.h"
int dapdebug = 0;

#line 79 "dap.tab.c" /* yacc.c:339  */

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
   by #include "dap.tab.h".  */
#ifndef YY_DAP_DAP_TAB_H_INCLUDED
# define YY_DAP_DAP_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif
#if YYDEBUG
extern int dapdebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    SCAN_ALIAS = 258,
    SCAN_ARRAY = 259,
    SCAN_ATTR = 260,
    SCAN_BYTE = 261,
    SCAN_CODE = 262,
    SCAN_DATASET = 263,
    SCAN_DATA = 264,
    SCAN_ERROR = 265,
    SCAN_FLOAT32 = 266,
    SCAN_FLOAT64 = 267,
    SCAN_GRID = 268,
    SCAN_INT16 = 269,
    SCAN_INT32 = 270,
    SCAN_MAPS = 271,
    SCAN_MESSAGE = 272,
    SCAN_SEQUENCE = 273,
    SCAN_STRING = 274,
    SCAN_STRUCTURE = 275,
    SCAN_UINT16 = 276,
    SCAN_UINT32 = 277,
    SCAN_URL = 278,
    SCAN_PTYPE = 279,
    SCAN_PROG = 280,
    WORD_WORD = 281,
    WORD_STRING = 282
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif



int dapparse (DAPparsestate* parsestate);

#endif /* !YY_DAP_DAP_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 157 "dap.tab.c" /* yacc.c:358  */

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
#define YYFINAL  9
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   369

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  36
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  34
/* YYNRULES -- Number of rules.  */
#define YYNRULES  106
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  201

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   282

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
       2,     2,     2,     2,    35,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    31,    30,
       2,    34,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    32,     2,    33,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    28,     2,    29,     2,     2,     2,     2,
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
      25,    26,    27
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    54,    54,    55,    56,    57,    58,    62,    66,    70,
      75,    81,    82,    88,    90,    92,    94,    97,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   115,   116,   120,
     121,   122,   123,   128,   129,   133,   136,   137,   142,   143,
     147,   148,   150,   152,   154,   156,   158,   160,   162,   164,
     166,   167,   172,   173,   177,   178,   182,   183,   187,   188,
     192,   193,   196,   197,   200,   201,   204,   205,   209,   210,
     214,   218,   219,   230,   234,   238,   238,   239,   239,   240,
     240,   241,   241,   247,   248,   249,   250,   251,   252,   253,
     254,   255,   256,   257,   258,   259,   260,   261,   262,   263,
     264,   265,   266,   267,   268,   269,   270
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 1
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "SCAN_ALIAS", "SCAN_ARRAY", "SCAN_ATTR",
  "SCAN_BYTE", "SCAN_CODE", "SCAN_DATASET", "SCAN_DATA", "SCAN_ERROR",
  "SCAN_FLOAT32", "SCAN_FLOAT64", "SCAN_GRID", "SCAN_INT16", "SCAN_INT32",
  "SCAN_MAPS", "SCAN_MESSAGE", "SCAN_SEQUENCE", "SCAN_STRING",
  "SCAN_STRUCTURE", "SCAN_UINT16", "SCAN_UINT32", "SCAN_URL", "SCAN_PTYPE",
  "SCAN_PROG", "WORD_WORD", "WORD_STRING", "'{'", "'}'", "';'", "':'",
  "'['", "']'", "'='", "','", "$accept", "start", "dataset", "attr", "err",
  "datasetbody", "declarations", "declaration", "base_type", "array_decls",
  "array_decl", "datasetname", "var_name", "attributebody", "attr_list",
  "attribute", "bytes", "int16", "uint16", "int32", "uint32", "float32",
  "float64", "strs", "urls", "url", "str_or_id", "alias", "errorbody",
  "errorcode", "errormsg", "errorptype", "errorprog", "name", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   123,   125,
      59,    58,    91,    93,    61,    44
};
# endif

#define YYPACT_NINF -91

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-91)))

#define YYTABLE_NINF -1

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
       6,   -91,   -91,   -91,   -91,     9,   -22,     7,   -16,   -91,
     -91,    10,   -91,   -91,   -91,    20,   -91,    37,   -91,   191,
      -6,    14,   -91,   -91,   -91,   -91,    17,   -91,   -91,    18,
     -91,    19,   -91,   -91,   -91,   271,   -91,   320,   -91,    27,
     -91,   -91,   320,   -91,   -91,   -91,   -91,   320,   320,   -91,
     320,   320,   -91,   -91,   -91,   320,   -91,   320,   320,   320,
     -91,   -91,   -91,   -91,   -91,    24,    43,    35,    39,    50,
      74,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   -91,    55,   -91,   -91,   -91,    60,    67,
      68,    70,    71,    73,   295,    77,    78,   295,   -91,   -91,
      65,    79,    66,    80,    76,    69,   127,   -91,     4,   -91,
     -91,   -20,   -91,   -13,   -91,   -12,   -91,   -10,   -91,    -9,
     -91,    32,   -91,   -91,   -91,    33,   -91,    34,    42,   -91,
     -91,   218,   -91,    81,    82,    75,    83,   346,   320,   320,
     -91,   -91,   159,   -91,   -91,    84,   -91,    88,   -91,    89,
     -91,    90,   -91,    91,   -91,   295,   -91,    92,   -91,    93,
     -91,   295,   -91,   -91,    95,    94,    96,   105,    97,   -91,
      98,   103,   100,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   102,   -91,    99,   -91,    12,   -91,   111,
     109,   -91,   -91,   -91,   -91,   118,   244,   -91,   320,   106,
     -91
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     6,     8,     7,     9,     0,     0,     0,     0,     1,
      11,     2,    37,    38,     4,    75,     5,     0,     3,     0,
       0,    77,    17,    18,    23,    24,     0,    19,    21,     0,
      26,     0,    20,    22,    25,     0,    12,     0,    51,    84,
      85,    86,    87,   103,    88,    89,    90,    91,    92,    93,
      94,    95,    96,   104,    97,    98,    99,   100,   101,   102,
     106,   105,    83,    36,    39,     0,     0,     0,     0,    79,
       0,    11,    11,    34,    84,    87,    91,    92,    94,    95,
      98,   100,   101,   102,     0,    33,    35,    27,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    40,    38,
       0,     0,     0,    81,     0,     0,     0,    10,     0,    73,
      52,     0,    62,     0,    64,     0,    54,     0,    58,     0,
      72,     0,    66,    71,    56,     0,    60,     0,     0,    68,
      70,     0,    76,     0,     0,     0,     0,     0,     0,     0,
      32,    13,     0,    28,    41,     0,    46,     0,    47,     0,
      42,     0,    44,     0,    48,     0,    43,     0,    45,     0,
      49,     0,    50,    78,     0,     0,     0,     0,     0,    27,
      83,     0,     0,    53,    63,    65,    55,    59,    67,    57,
      61,    69,    80,     0,    74,     0,    15,     0,    29,     0,
       0,    82,    11,    14,    30,     0,     0,    31,     0,     0,
      16
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -91,   -91,   -91,   -91,   -91,   -91,   -69,   -15,   -91,   -17,
     -91,   -91,   -37,   -91,    54,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   -91,   -91,    -7,   -90,   -91,   -91,   -91,
     -91,   -91,   -91,   -18
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     5,     6,     7,     8,    11,    17,    36,    37,   108,
     143,    84,    85,    14,    19,    64,   111,   117,   125,   119,
     127,   113,   115,   121,   128,   129,   130,    65,    16,    21,
      69,   103,   136,    86
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_uint8 yytable[] =
{
      87,    66,   105,   106,   122,   140,    10,     1,    12,     9,
     144,     2,    15,   140,     3,   145,     4,   146,   148,    18,
     150,   152,   147,   149,    89,   151,   153,    20,    67,    90,
      91,    68,    92,    93,   141,    13,   142,    94,    22,    95,
      96,    97,   193,    23,   142,    70,    71,    72,    24,    25,
      26,    27,    28,    88,    98,    29,    30,    31,    32,    33,
      34,   100,   154,   156,   158,   178,    35,   155,   157,   159,
      22,    99,   160,   101,   102,    23,   123,   161,   104,   123,
      24,    25,    26,    27,    28,   107,   109,    29,    30,    31,
      32,    33,    34,   110,   112,   132,   114,   116,   138,   118,
     134,   168,   169,   124,   126,   135,   133,   137,   164,   165,
     173,   163,   166,    66,   174,   175,   176,   177,   179,   180,
     183,   185,   167,   196,   172,   182,   184,   186,    22,   189,
     192,   188,   191,    23,   190,   195,   200,   123,    24,    25,
      26,    27,    28,   123,   194,    29,    30,    31,    32,    33,
      34,   197,   187,   131,   181,     0,   139,     0,     0,     0,
       0,   199,    74,    40,    41,    75,    43,    44,    45,    46,
      76,    77,    49,    78,    79,    52,    53,    54,    80,    56,
      81,    82,    83,    60,    61,   170,     0,     0,     0,     0,
       0,     0,    38,   171,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,     0,    38,
      63,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    54,    55,    56,    57,
      58,    59,    60,    61,    62,    22,     0,   162,     0,     0,
      23,     0,     0,     0,     0,    24,    25,    26,    27,    28,
       0,     0,    29,    30,    31,    32,    33,    34,     0,     0,
       0,     0,    73,   198,    74,    40,    41,    75,    43,    44,
      45,    46,    76,    77,    49,    78,    79,    52,    53,    54,
      80,    56,    81,    82,    83,    60,    61,    62,    74,    40,
      41,    75,    43,    44,    45,    46,    76,    77,    49,    78,
      79,    52,    53,    54,    80,    56,    81,    82,    83,    60,
      61,    62,   120,    74,    40,    41,    75,    43,    44,    45,
      46,    76,    77,    49,    78,    79,    52,    53,    54,    80,
      56,    81,    82,    83,    60,    61,    62,    22,     0,     0,
       0,     0,    23,     0,     0,     0,     0,    24,    25,    26,
      27,    28,     0,     0,    29,    30,    31,    32,    33,    34
};

static const yytype_int16 yycheck[] =
{
      37,    19,    71,    72,    94,     1,    28,     1,     1,     0,
      30,     5,    28,     1,     8,    35,    10,    30,    30,     9,
      30,    30,    35,    35,    42,    35,    35,     7,    34,    47,
      48,    17,    50,    51,    30,    28,    32,    55,     1,    57,
      58,    59,    30,     6,    32,    28,    28,    28,    11,    12,
      13,    14,    15,    26,    30,    18,    19,    20,    21,    22,
      23,    26,    30,    30,    30,   155,    29,    35,    35,    35,
       1,    28,    30,    34,    24,     6,    94,    35,     4,    97,
      11,    12,    13,    14,    15,    30,    26,    18,    19,    20,
      21,    22,    23,    26,    26,    30,    26,    26,    29,    26,
      34,   138,   139,    26,    26,    25,    27,    31,    26,    34,
      26,    30,    29,   131,    26,    26,    26,    26,    26,    26,
      26,    16,   137,   192,   142,    30,    30,    30,     1,    26,
      31,    33,    30,     6,    34,    26,    30,   155,    11,    12,
      13,    14,    15,   161,    33,    18,    19,    20,    21,    22,
      23,    33,   169,    99,   161,    -1,    29,    -1,    -1,    -1,
      -1,   198,     3,     4,     5,     6,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    -1,    -1,    -1,    -1,
      -1,    -1,     1,    34,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    -1,     1,
      29,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,     1,    -1,    29,    -1,    -1,
       6,    -1,    -1,    -1,    -1,    11,    12,    13,    14,    15,
      -1,    -1,    18,    19,    20,    21,    22,    23,    -1,    -1,
      -1,    -1,     1,    29,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,     3,     4,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,     1,    -1,    -1,
      -1,    -1,     6,    -1,    -1,    -1,    -1,    11,    12,    13,
      14,    15,    -1,    -1,    18,    19,    20,    21,    22,    23
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     1,     5,     8,    10,    37,    38,    39,    40,     0,
      28,    41,     1,    28,    49,    28,    64,    42,     9,    50,
       7,    65,     1,     6,    11,    12,    13,    14,    15,    18,
      19,    20,    21,    22,    23,    29,    43,    44,     1,     3,
       4,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    29,    51,    63,    69,    34,    17,    66,
      28,    28,    28,     1,     3,     6,    11,    12,    14,    15,
      19,    21,    22,    23,    47,    48,    69,    48,    26,    69,
      69,    69,    69,    69,    69,    69,    69,    69,    30,    28,
      26,    34,    24,    67,     4,    42,    42,    30,    45,    26,
      26,    52,    26,    57,    26,    58,    26,    53,    26,    55,
      27,    59,    62,    69,    26,    54,    26,    56,    60,    61,
      62,    50,    30,    27,    34,    25,    68,    31,    29,    29,
       1,    30,    32,    46,    30,    35,    30,    35,    30,    35,
      30,    35,    30,    35,    30,    35,    30,    35,    30,    35,
      30,    35,    29,    30,    26,    34,    29,    43,    48,    48,
      26,    34,    69,    26,    26,    26,    26,    26,    62,    26,
      26,    61,    30,    26,    30,    16,    30,    45,    33,    26,
      34,    30,    31,    30,    33,    26,    42,    33,    29,    48,
      30
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    36,    37,    37,    37,    37,    37,    38,    39,    40,
      41,    42,    42,    43,    43,    43,    43,    43,    44,    44,
      44,    44,    44,    44,    44,    44,    44,    45,    45,    46,
      46,    46,    46,    47,    47,    48,    49,    49,    50,    50,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      51,    51,    52,    52,    53,    53,    54,    54,    55,    55,
      56,    56,    57,    57,    58,    58,    59,    59,    60,    60,
      61,    62,    62,    63,    64,    65,    65,    66,    66,    67,
      67,    68,    68,    69,    69,    69,    69,    69,    69,    69,
      69,    69,    69,    69,    69,    69,    69,    69,    69,    69,
      69,    69,    69,    69,    69,    69,    69
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     3,     2,     2,     1,     1,     1,     1,
       5,     0,     2,     4,     7,     6,    11,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     2,     3,
       4,     5,     1,     1,     1,     1,     3,     1,     0,     2,
       2,     4,     4,     4,     4,     4,     4,     4,     4,     4,
       4,     1,     1,     3,     1,     3,     1,     3,     1,     3,
       1,     3,     1,     3,     1,     3,     1,     3,     1,     3,
       1,     1,     1,     3,     7,     0,     4,     0,     4,     0,
       4,     0,     4,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1
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
      yyerror (parsestate, YY_("syntax error: cannot back up")); \
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
                  Type, Value, parsestate); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, DAPparsestate* parsestate)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  YYUSE (parsestate);
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
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep, DAPparsestate* parsestate)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep, parsestate);
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
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule, DAPparsestate* parsestate)
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
                                              , parsestate);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule, parsestate); \
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
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep, DAPparsestate* parsestate)
{
  YYUSE (yyvaluep);
  YYUSE (parsestate);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/*----------.
| yyparse.  |
`----------*/

int
yyparse (DAPparsestate* parsestate)
{
/* The lookahead symbol.  */
int yychar;


/* The semantic value of the lookahead symbol.  */
/* Default value used for initialization, for pacifying older GCCs
   or non-GCC compilers.  */
YY_INITIAL_VALUE (static YYSTYPE yyval_default;)
YYSTYPE yylval YY_INITIAL_VALUE (= yyval_default);

    /* Number of syntax errors so far.  */
    int yynerrs;

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
      yychar = yylex (&yylval, parsestate);
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
        case 6:
#line 58 "dap.y" /* yacc.c:1646  */
    {dap_unrecognizedresponse(parsestate); YYABORT;}
#line 1412 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 7:
#line 63 "dap.y" /* yacc.c:1646  */
    {dap_tagparse(parsestate,SCAN_DATASET);}
#line 1418 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 8:
#line 67 "dap.y" /* yacc.c:1646  */
    {dap_tagparse(parsestate,SCAN_ATTR);}
#line 1424 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 9:
#line 71 "dap.y" /* yacc.c:1646  */
    {dap_tagparse(parsestate,SCAN_ERROR);}
#line 1430 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 10:
#line 76 "dap.y" /* yacc.c:1646  */
    {dap_datasetbody(parsestate,(yyvsp[-1]),(yyvsp[-3]));}
#line 1436 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 11:
#line 81 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_declarations(parsestate,null,null);}
#line 1442 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 12:
#line 82 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_declarations(parsestate,(yyvsp[-1]),(yyvsp[0]));}
#line 1448 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 13:
#line 89 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_makebase(parsestate,(yyvsp[-2]),(yyvsp[-3]),(yyvsp[-1]));}
#line 1454 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 14:
#line 91 "dap.y" /* yacc.c:1646  */
    {if(((yyval)=dap_makestructure(parsestate,(yyvsp[-2]),(yyvsp[-1]),(yyvsp[-4])))==null) {YYABORT;}}
#line 1460 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 15:
#line 93 "dap.y" /* yacc.c:1646  */
    {if(((yyval)=dap_makesequence(parsestate,(yyvsp[-1]),(yyvsp[-3])))==null) {YYABORT;}}
#line 1466 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 16:
#line 96 "dap.y" /* yacc.c:1646  */
    {if(((yyval)=dap_makegrid(parsestate,(yyvsp[-1]),(yyvsp[-6]),(yyvsp[-3])))==null) {YYABORT;}}
#line 1472 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 17:
#line 98 "dap.y" /* yacc.c:1646  */
    {dapsemanticerror(parsestate,OC_EBADTYPE,"Unrecognized type"); YYABORT;}
#line 1478 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 18:
#line 103 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_BYTE;}
#line 1484 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 19:
#line 104 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_INT16;}
#line 1490 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 20:
#line 105 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_UINT16;}
#line 1496 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 21:
#line 106 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_INT32;}
#line 1502 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 22:
#line 107 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_UINT32;}
#line 1508 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 23:
#line 108 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_FLOAT32;}
#line 1514 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 24:
#line 109 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_FLOAT64;}
#line 1520 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 25:
#line 110 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_URL;}
#line 1526 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 26:
#line 111 "dap.y" /* yacc.c:1646  */
    {(yyval)=(Object)SCAN_STRING;}
#line 1532 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 27:
#line 115 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_arraydecls(parsestate,null,null);}
#line 1538 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 28:
#line 116 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_arraydecls(parsestate,(yyvsp[-1]),(yyvsp[0]));}
#line 1544 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 29:
#line 120 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_arraydecl(parsestate,null,(yyvsp[-1]));}
#line 1550 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 30:
#line 121 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_arraydecl(parsestate,null,(yyvsp[-1]));}
#line 1556 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 31:
#line 122 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_arraydecl(parsestate,(yyvsp[-3]),(yyvsp[-1]));}
#line 1562 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 32:
#line 124 "dap.y" /* yacc.c:1646  */
    {dapsemanticerror(parsestate,OC_EDIMSIZE,"Illegal dimension declaration"); YYABORT;}
#line 1568 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 33:
#line 128 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[0]);}
#line 1574 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 34:
#line 130 "dap.y" /* yacc.c:1646  */
    {dapsemanticerror(parsestate,OC_EDDS,"Illegal dataset declaration"); YYABORT;}
#line 1580 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 35:
#line 133 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[0]);}
#line 1586 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 36:
#line 136 "dap.y" /* yacc.c:1646  */
    {dap_attributebody(parsestate,(yyvsp[-1]));}
#line 1592 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 37:
#line 138 "dap.y" /* yacc.c:1646  */
    {dapsemanticerror(parsestate,OC_EDAS,"Illegal DAS body"); YYABORT;}
#line 1598 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 38:
#line 142 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrlist(parsestate,null,null);}
#line 1604 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 39:
#line 143 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrlist(parsestate,(yyvsp[-1]),(yyvsp[0]));}
#line 1610 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 40:
#line 147 "dap.y" /* yacc.c:1646  */
    {(yyval)=null;}
#line 1616 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 41:
#line 149 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_BYTE);}
#line 1622 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 42:
#line 151 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_INT16);}
#line 1628 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 43:
#line 153 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_UINT16);}
#line 1634 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 44:
#line 155 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_INT32);}
#line 1640 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 45:
#line 157 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_UINT32);}
#line 1646 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 46:
#line 159 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_FLOAT32);}
#line 1652 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 47:
#line 161 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_FLOAT64);}
#line 1658 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 48:
#line 163 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_STRING);}
#line 1664 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 49:
#line 165 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attribute(parsestate,(yyvsp[-2]),(yyvsp[-1]),(Object)SCAN_URL);}
#line 1670 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 50:
#line 166 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrset(parsestate,(yyvsp[-3]),(yyvsp[-1]));}
#line 1676 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 51:
#line 168 "dap.y" /* yacc.c:1646  */
    {dapsemanticerror(parsestate,OC_EDAS,"Illegal attribute"); YYABORT;}
#line 1682 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 52:
#line 172 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_BYTE);}
#line 1688 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 53:
#line 174 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_BYTE);}
#line 1694 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 54:
#line 177 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_INT16);}
#line 1700 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 55:
#line 179 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_INT16);}
#line 1706 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 56:
#line 182 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_UINT16);}
#line 1712 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 57:
#line 184 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_UINT16);}
#line 1718 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 58:
#line 187 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_INT32);}
#line 1724 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 59:
#line 189 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_INT32);}
#line 1730 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 60:
#line 192 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_UINT32);}
#line 1736 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 61:
#line 193 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_UINT32);}
#line 1742 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 62:
#line 196 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_FLOAT32);}
#line 1748 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 63:
#line 197 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_FLOAT32);}
#line 1754 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 64:
#line 200 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_FLOAT64);}
#line 1760 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 65:
#line 201 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_FLOAT64);}
#line 1766 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 66:
#line 204 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_STRING);}
#line 1772 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 67:
#line 205 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_STRING);}
#line 1778 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 68:
#line 209 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,null,(yyvsp[0]),(Object)SCAN_URL);}
#line 1784 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 69:
#line 210 "dap.y" /* yacc.c:1646  */
    {(yyval)=dap_attrvalue(parsestate,(yyvsp[-2]),(yyvsp[0]),(Object)SCAN_URL);}
#line 1790 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 70:
#line 214 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[0]);}
#line 1796 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 71:
#line 218 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[0]);}
#line 1802 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 72:
#line 219 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[0]);}
#line 1808 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 73:
#line 230 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[-1]); (yyval)=(yyvsp[0]); (yyval)=null;}
#line 1814 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 74:
#line 235 "dap.y" /* yacc.c:1646  */
    {dap_errorbody(parsestate,(yyvsp[-5]),(yyvsp[-4]),(yyvsp[-3]),(yyvsp[-2]));}
#line 1820 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 75:
#line 238 "dap.y" /* yacc.c:1646  */
    {(yyval)=null;}
#line 1826 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 76:
#line 238 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[-1]);}
#line 1832 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 77:
#line 239 "dap.y" /* yacc.c:1646  */
    {(yyval)=null;}
#line 1838 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 78:
#line 239 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[-1]);}
#line 1844 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 79:
#line 240 "dap.y" /* yacc.c:1646  */
    {(yyval)=null;}
#line 1850 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 80:
#line 240 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[-1]);}
#line 1856 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 81:
#line 241 "dap.y" /* yacc.c:1646  */
    {(yyval)=null;}
#line 1862 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 82:
#line 241 "dap.y" /* yacc.c:1646  */
    {(yyval)=(yyvsp[-1]);}
#line 1868 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 83:
#line 247 "dap.y" /* yacc.c:1646  */
    {(yyval)=dapdecode(parsestate->lexstate,(yyvsp[0]));}
#line 1874 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 84:
#line 248 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1880 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 85:
#line 249 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1886 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 86:
#line 250 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1892 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 87:
#line 251 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1898 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 88:
#line 252 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1904 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 89:
#line 253 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1910 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 90:
#line 254 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1916 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 91:
#line 255 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1922 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 92:
#line 256 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1928 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 93:
#line 257 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1934 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 94:
#line 258 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1940 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 95:
#line 259 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1946 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 96:
#line 260 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1952 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 97:
#line 261 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1958 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 98:
#line 262 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1964 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 99:
#line 263 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1970 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 100:
#line 264 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1976 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 101:
#line 265 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1982 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 102:
#line 266 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1988 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 103:
#line 267 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 1994 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 104:
#line 268 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 2000 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 105:
#line 269 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 2006 "dap.tab.c" /* yacc.c:1646  */
    break;

  case 106:
#line 270 "dap.y" /* yacc.c:1646  */
    {(yyval)=strdup((yyvsp[0]));}
#line 2012 "dap.tab.c" /* yacc.c:1646  */
    break;


#line 2016 "dap.tab.c" /* yacc.c:1646  */
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
      yyerror (parsestate, YY_("syntax error"));
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
        yyerror (parsestate, yymsgp);
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
                      yytoken, &yylval, parsestate);
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
                  yystos[yystate], yyvsp, parsestate);
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
  yyerror (parsestate, YY_("memory exhausted"));
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
                  yytoken, &yylval, parsestate);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp, parsestate);
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
#line 273 "dap.y" /* yacc.c:1906  */

