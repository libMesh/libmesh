/* A Bison parser, made by GNU Bison 3.0.  */

/* Bison interface for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2013 Free Software Foundation, Inc.

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
