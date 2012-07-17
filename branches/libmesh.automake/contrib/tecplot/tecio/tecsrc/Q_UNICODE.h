/*
 * NOTICE and LICENSE for Tecplot Input/Output Library (TecIO) - OpenFOAM
 *
 * Copyright (C) 1988-2009 Tecplot, Inc.  All rights reserved worldwide.
 *
 * Tecplot hereby grants OpenCFD limited authority to distribute without
 * alteration the source code to the Tecplot Input/Output library, known 
 * as TecIO, as part of its distribution of OpenFOAM and the 
 * OpenFOAM_to_Tecplot converter.  Users of this converter are also hereby
 * granted access to the TecIO source code, and may redistribute it for the
 * purpose of maintaining the converter.  However, no authority is granted
 * to alter the TecIO source code in any form or manner.
 *
 * This limited grant of distribution does not supersede Tecplot, Inc.'s 
 * copyright in TecIO.  Contact Tecplot, Inc. for further information.
 * 
 * Tecplot, Inc.
 * 3535 Factoria Blvd, Ste. 550
 * Bellevue, WA 98006, USA
 * Phone: +1 425 653 1200
 * http://www.tecplot.com/
 *
 */
/*
******************************************************************
******************************************************************
*******                                                   ********
******  (C) 1988-2008 Tecplot, Inc.                        *******
*******                                                   ********
******************************************************************
******************************************************************
*/


#if !defined Q_UNICODE_H_
# define Q_UNICODE_H_

#if defined EXTERN
#undef EXTERN
#endif
#if defined Q_UNICODEMODULE
#define EXTERN
#else
#define EXTERN extern
#endif

namespace tecplot
{
namespace strutil
{

// functions
Boolean_t IsValidUtf8LeadByte(Byte_t ch);
Boolean_t IsValidUtf8ContinuingByte(Byte_t ch);
Boolean_t IsValidUtf8Byte(Byte_t ch);

Boolean_t IsValidUtf8String(const char *str);
Boolean_t ShouldConvertWideStringToUtf8String(const wchar_t *str);
void InitTranslatedStrings();
void CleanUpTranslatedStrings();

Boolean_t IsNullOrZeroLengthString(const char *S);
Boolean_t IsNullOrZeroLengthString(tecplot::strutil::TranslatedString TS);

Boolean_t IsEmptyString(const char *S);
Boolean_t IsEmptyString(tecplot::strutil::TranslatedString S);
Boolean_t IsEmptyString(const wchar_t* S);

#if defined MSWIN

std::string  LookUpTranslation(std::string& strEnglish);
void MsWinInitTranslatedStrings();

std::string    WStringToString(std::wstring str);
std::wstring   StringToWString(std::string str);

std::wstring   MultiByteToWideChar(const char *Utf8Str,
                                   unsigned int    CodePage);

std::string    WideCharToMultiByte(const wchar_t *WideStr,
                                   unsigned int    CodePage);

// Conversion
std::string    WideCharToUtf8(const wchar_t* str);
std::wstring   Utf8ToWideChar(const char *str);
char *getenv(const char *str);

#endif

}
}

#endif
