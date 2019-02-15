/*
 * Copyright 1998-2015 University Corporation for Atmospheric Research/Unidata
 *  See the LICENSE file for more information.
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netcdf.h"
#include "ncutf8.h"
/*

This test is taken from the UTF-8 decoder
capability and stress test file created by

Markus Kuhn <http://www.cl.cam.ac.uk/~mgk25/> - 2015-08-28 - CC BY 4.0

This test file can help you examine, how your UTF-8 decoder handles
various types of correct, malformed, or otherwise interesting UTF-8
sequences. This file is not meant to be a conformance test. It does
not prescribe any particular outcome. Therefore, there is no way to
"pass" or "fail" this test file, even though the text does suggest a
preferable decoder behaviour at some places. Its aim is, instead, to
help you think about, and test, the behaviour of your UTF-8 decoder on a
systematic collection of unusual inputs. Experience so far suggests
that most first-time authors of UTF-8 decoders find at least one
serious problem in their decoder using this file.

The test lines below cover boundary conditions, malformed UTF-8
sequences, as well as correctly encoded UTF-8 sequences of Unicode code
points that should never occur in a correct UTF-8 file.

According to ISO 10646-1:2000, sections D.7 and 2.3c, a device
receiving UTF-8 shall interpret a "malformed sequence in the same way
that it interprets a character that is outside the adopted subset" and
"characters that are not within the adopted subset shall be indicated
to the user" by a receiving device. One commonly used approach in
UTF-8 decoders is to replace any malformed UTF-8 sequence by a
replacement character (U+FFFD), which looks a bit like an inverted
question mark, or a similar symbol. It might be a good idea to
visually distinguish a malformed UTF-8 sequence from a correctly
encoded Unicode character that is just not available in the current
font but otherwise fully legal, even though ISO 10646-1 doesn't
mandate this. In any case, just ignoring malformed sequences or
unavailable characters does not conform to ISO 10646, will make
debugging more difficult, and can lead to user confusion.

Please check, whether a malformed UTF-8 sequence is (1) represented at
all, (2) represented by exactly one single replacement character (or
equivalent signal), and (3) the following quotation mark after an
illegal UTF-8 sequence is correctly displayed, i.e. proper
resynchronization takes place immediately after any malformed
sequence. This file says "THE END" in the last line, so if you don't
see that, your decoder crashed somehow before, which should always be
cause for concern.

All lines in this file are exactly 79 characters long (plus the line
feed). In addition, all lines end with "|", except for the two test
lines 2.1.1 and 2.2.1, which contain non-printable ASCII controls
U+0000 and U+007F. If you display this file with a fixed-width font,
these "|" characters should all line up in column 79 (right margin).
This allows you to test quickly, whether your UTF-8 decoder finds the
correct number of characters in every line, that is whether each
malformed sequences is replaced by a single replacement character.

Note that, as an alternative to the notion of malformed sequence used
here, it is also a perfectly acceptable (and in some situations even
preferable) solution to represent each individual byte of a malformed
sequence with a replacement character. If you follow this strategy in
your decoder, then please ignore the "|" column.
*/


struct Test {
    int xfail;
    const char* id;
    const char* description;
    const char* data;
};
#define NULLTEST {0,NULL,NULL,NULL}

/* The following tests are in envv form */

/*1  Some correct UTF-8 text
     You should see the Greek word 'kosme':
*/
static const struct Test utf8ok[] = {
{0, "1.1.1", "Greek word 'kosme'",
"Îºá½¹ÏƒÎ¼Îµ"},
NULLTEST
};

static const struct Test utf8boundary[] = {

/*2  Boundary condition test */
/*2.1  First possible sequence of a certain length */
{0,"2.1.1", "1 byte  (U-00000000)",        "\000"},
{0,"2.1.2", "2 bytes (U-00000080)",        "Â€"},
{0,"2.1.3", "3 bytes (U-00000800)",        "à €"},
{0,"2.1.4", "4 bytes (U-00010000)",        "ð€€"},
{1,"2.1.5", "5 bytes (U-00200000)",        "øˆ€€€"},
{1,"2.1.6", "6 bytes (U-04000000)",        "ü„€€€€"},

/*2.2  Last possible sequence of a certain length*/
{0,"2.2.1", "1 byte  (U-0000007F)",        ""},
{0,"2.2.2", "2 bytes (U-000007FF)",        "ß¿"},
{0,"2.2.3", "3 bytes (U-0000FFFF)",        "ï¿¿"}, /*See 5.3.2 */
{1,"2.2.4", "4 bytes (U-001FFFFF)",        "÷¿¿¿"},
{1,"2.2.5", "5 bytes (U-03FFFFFF)",        "û¿¿¿¿"},
{1,"2.2.6", "6 bytes (U-7FFFFFFF)",        "ý¿¿¿¿¿"},

/*2.3  Other boundary conditions*/

{0,"2.3.1", "U-0000D7FF = ed 9f bf", "íŸ¿"},
{0,"2.3.2", "U-0000E000 = ee 80 80", "î€€"},
{0,"2.3.3", "U-0000FFFD = ef bf bd", "ï¿½"},
{0,"2.3.4", "U-0010FFFF = f4 8f bf bf", "ô¿¿"},
{1,"2.3.5", "U-00110000 = f4 90 80 80", "ô€€"},
NULLTEST
};

static const struct Test utf8bad[] = {

/*3  Malformed sequences*/

/*3.1  Unexpected continuation bytes
       Each unexpected continuation byte should be separately signalled
       as a malformed sequence of its own.
*/
{1,"3.1.1", "First continuation byte 0x80", "€"},
{1,"3.1.2", "Last  continuation byte 0xbf", "¿"},

{1,"3.1.3", "2 continuation bytes", "€¿"},
{1,"3.1.4", "3 continuation bytes", "€¿€"},
{1,"3.1.5", "4 continuation bytes", "€¿€¿"},
{1,"3.1.6", "5 continuation bytes", "€¿€¿€"},
{1,"3.1.7", "6 continuation bytes", "€¿€¿€¿"},
{1,"3.1.8", "7 continuation bytes", "€¿€¿€¿€"},
{1,"3.1.9", "Sequence of all 64 possible continuation bytes (0x80-0xbf)",
   "€‚ƒ„…†‡ˆ‰Š‹ŒŽ‘’“”•–—˜™š›œžŸ ¡¢£¤¥¦§¨©ª«¬­®¯°±²³´µ¶·¸¹º»¼½¾¿"
},

/*3.2  Lonely start characters*/

/*3.2.1  All 32 first bytes of 2-byte sequences (0xc0-0xdf),
       each followed by a space character*/

{1,"3.2.1", "All 32 first bytes of 2-byte sequences",
"À Á Â Ã Ä Å Æ Ç È É Ê Ë Ì Í Î Ï Ð Ñ Ò Ó Ô Õ Ö × Ø Ù Ú Û Ü Ý Þ ß "
},

/*3.2.2  All 16 first bytes of 3-byte sequences (0xe0-0xef),
       each followed by a space character:*/
{1,"3.2.2", "All 16 first bytes of 3-byte sequences",
"à á â ã ä å æ ç è é ê ë ì í î ï "
},

/*3.2.3  All 8 first bytes of 4-byte sequences (0xf0-0xf7),
       each followed by a space character:*/
{1,"3.2.3", "All 8 first bytes of 4-byte sequences",
   "ð ñ ò ó ô õ ö ÷ "
},

/*3.2.4  All 4 first bytes of 5-byte sequences (0xf8-0xfb),
       each followed by a space character:*/
{1,"3.2.4", "All 4 first bytes of 5-byte sequences",
   "ø ù ú û "
},

/*3.2.5  All 2 first bytes of 6-byte sequences (0xfc-0xfd),
       each followed by a space character:*/
{1,"3.2.5", "All 2 first bytes of 6-byte sequences",
   "ü ý "
},

/*3.3  Sequences with last continuation byte missing
All bytes of an incomplete sequence should be signalled as a single
malformed sequence, i.e., you should see only a single replacement
character in each of the next 10 tests. (Characters as in section 2)
*/
{1,"3.3.1", "2-byte sequence with last byte missing (U+0000)",     "À"},
{1,"3.3.2", "3-byte sequence with last byte missing (U+0000)",     "à€"},
{1,"3.3.3", "4-byte sequence with last byte missing (U+0000)",     "ð€€"},
{1,"3.3.4", "5-byte sequence with last byte missing (U+0000)",     "ø€€€"},
{1,"3.3.5", "6-byte sequence with last byte missing (U+0000)",     "ü€€€€"},
{1,"3.3.6", "2-byte sequence with last byte missing (U-000007FF)", "ß"},
{1,"3.3.7", "3-byte sequence with last byte missing (U-0000FFFF)", "ï¿"},
{1,"3.3.8", "4-byte sequence with last byte missing (U-001FFFFF)", "÷¿¿"},
{1,"3.3.9", "5-byte sequence with last byte missing (U-03FFFFFF)", "û¿¿¿"},
{1,"3.3.10", "6-byte sequence with last byte missing (U-7FFFFFFF)", "ý¿¿¿¿"},

/*3.4 Concatenation of incomplete sequences
All the 10 sequences of 3.3 concatenated; you should see 10 malformed
sequences being signalled:
*/
{1, "3.4.1", "All the 10 sequences of 3.3 concatenated",
   "Àà€ð€€ø€€€ü€€€€ßï¿÷¿¿û¿¿¿ý¿¿¿¿"
},

/*3.5  Impossible bytes
The following two bytes cannot appear in a correct UTF-8 string
*/

{1,"3.5.1", "fe", "þ"},
{1,"3.5.2", "ff", "ÿ"},
{1,"3.5.3", "fe fe ff ff", "þþÿÿ"},

/*
4  Overlong sequences

The following sequences are not malformed according to the letter of
the Unicode 2.0 standard. However, they are longer then necessary and
a correct UTF-8 encoder is not allowed to produce them. A "safe UTF-8
decoder" should reject them just like malformed sequences for two
reasons: (1) It helps to debug applications if overlong sequences are
not treated as valid representations of characters, because this helps
to spot problems more quickly. (2) Overlong sequences provide
alternative representations of characters, that could maliciously be
used to bypass filters that check only for ASCII characters. For
instance, a 2-byte encoded line feed (LF) would not be caught by a
line counter that counts only 0x0a bytes, but it would still be
processed as a line feed by an unsafe UTF-8 decoder later in the
pipeline. From a security point of view, ASCII compatibility of UTF-8
sequences means also, that ASCII characters are *only* allowed to be
represented by ASCII bytes in the range 0x00-0x7f. To ensure this
aspect of ASCII compatibility, use only "safe UTF-8 decoders" that
reject overlong UTF-8 sequences for which a shorter encoding exists.
*/

/*4.1  Examples of an overlong ASCII character

With a safe UTF-8 decoder, all of the following five overlong
representations of the ASCII character slash ("/") should be rejected
like a malformed UTF-8 sequence, for instance by substituting it with
a replacement character. If you see a slash below, you do not have a
safe UTF-8 decoder!
*/

{1,"4.1.1", "U+002F = c0 af             ", "À¯"},
{1,"4.1.2", "U+002F = e0 80 af          ", "à€¯"},
{1,"4.1.3", "U+002F = f0 80 80 af       ", "ð€€¯"},
{1,"4.1.4", "U+002F = f8 80 80 80 af    ", "ø€€€¯"},
{1,"4.1.5", "U+002F = fc 80 80 80 80 af ", "ü€€€€¯"},

/*4.2  Maximum overlong sequences

Below you see the highest Unicode value that is still resulting in an
overlong sequence if represented with the given number of bytes. This
is a boundary test for safe UTF-8 decoders. All five characters should
be rejected like malformed UTF-8 sequences.
*/

{1,"4.2.1", "U-0000007F = c1 bf             ", "Á¿"},
{1,"4.2.2", "U-000007FF = e0 9f bf          ", "àŸ¿"},
{1,"4.2.3", "U-0000FFFF = f0 8f bf bf       ", "ð¿¿"},
{1,"4.2.4", "U-001FFFFF = f8 87 bf bf bf    ", "ø‡¿¿¿"},
{1,"4.2.5", "U-03FFFFFF = fc 83 bf bf bf bf ", "üƒ¿¿¿¿"},

/*
4.3  Overlong representation of the NUL character

The following five sequences should also be rejected like malformed
UTF-8 sequences and should not be treated like the ASCII NUL
character.
*/

{1,"4.3.1", "U+0000 = c0 80             ", "À€"},
{1,"4.3.2", "U+0000 = e0 80 80          ", "à€€"},
{1,"4.3.3", "U+0000 = f0 80 80 80       ", "ð€€€"},
{1,"4.3.4", "U+0000 = f8 80 80 80 80    ", "ø€€€€"},
{1,"4.3.5", "U+0000 = fc 80 80 80 80 80 ", "ü€€€€€"},

/*
5  Illegal code positions

The following UTF-8 sequences should be rejected like malformed
sequences, because they never represent valid ISO 10646 characters and
a UTF-8 decoder that accepts them might introduce security problems
comparable to overlong UTF-8 sequences.
*/
/*5.1 Single UTF-16 surrogates*/

{1,"5.1.1", "U+D800 = ed a0 80 ", "í €"},
{1,"5.1.2", "U+DB7F = ed ad bf ", "í­¿"},
{1,"5.1.3", "U+DB80 = ed ae 80 ", "í®€"},
{1,"5.1.4", "U+DBFF = ed af bf ", "í¯¿"},
{1,"5.1.5", "U+DC00 = ed b0 80 ", "í°€"},
{1,"5.1.6", "U+DF80 = ed be 80 ", "í¾€"},
{1,"5.1.7", "U+DFFF = ed bf bf ", "í¿¿"},

/*5.2 Paired UTF-16 surrogates */

{1,"5.2.1", "U+D800 U+DC00 = ed a0 80 ed b0 80 ", "í €í°€"},
{1,"5.2.2", "U+D800 U+DFFF = ed a0 80 ed bf bf ", "í €í¿¿"},
{1,"5.2.3", "U+DB7F U+DC00 = ed ad bf ed b0 80 ", "í­¿í°€"},
{1,"5.2.4", "U+DB7F U+DFFF = ed ad bf ed bf bf ", "í­¿í¿¿"},
{1,"5.2.5", "U+DB80 U+DC00 = ed ae 80 ed b0 80 ", "í®€í°€"},
{1,"5.2.6", "U+DB80 U+DFFF = ed ae 80 ed bf bf ", "í®€í¿¿"},
{1,"5.2.7", "U+DBFF U+DC00 = ed af bf ed b0 80 ", "í¯¿í°€"},
{1,"5.2.8", "U+DBFF U+DFFF = ed af bf ed bf bf ", "í¯¿í¿¿"},
NULLTEST
};

/*5.3 Noncharacter code positions

The following "noncharacters" are "reserved for internal use" by
applications, and according to older versions of the Unicode Standard
"should never be interchanged". Unicode Corrigendum #9 dropped the
latter restriction. Nevertheless, their presence in incoming UTF-8 data
can remain a potential security risk, depending on what use is made of
these codes subsequently. Examples of such internal use:

 - Some file APIs with 16-bit characters may use the integer value -1
   = U+FFFF to signal an end-of-file (EOF) or error condition.

 - In some UTF-16 receivers, code point U+FFFE might trigger a
   byte-swap operation (to convert between UTF-16LE and UTF-16BE).

With such internal use of noncharacters, it may be desirable and safer
to block those code points in UTF-8 decoders, as they should never
occur legitimately in incoming UTF-8 data, and could trigger unsafe
behaviour in subsequent processing.

Particularly problematic noncharacters in 16-bit applications:
*/

static const struct Test utf8problematic[] = {
{0,"5.3.1", "U+FFFE = ef bf be ", "ï¿¾"},
{0,"5.3.2", "U+FFFF = ef bf bf ", "ï¿¿"},
NULLTEST
};

/* Other (utf16) noncharacters: */
static const struct Test utf8nonchars[] = {
{0,"5.3.3", "U+FDD0 .. U+FDEF ",
"ï·ï·‘ï·’ï·“ï·”ï·•ï·–ï·—ï·˜ï·™ï·šï·›ï·œï·ï·žï·Ÿï· ï·¡ï·¢ï·£ï·¤ï·¥ï·¦ï·§ï·¨ï·©ï·ªï·«ï·¬ï·­ï·®ï·¯"
},

/* Do not understand this test; it passes, but should it? */
{0,"5.3.4", "U+nFFFE U+nFFFF (for n = 1..10)",
       "ðŸ¿¾ðŸ¿¿ð¯¿¾ð¯¿¿ð¿¿¾ð¿¿¿ñ¿¾ñ¿¿ñŸ¿¾ñŸ¿¿ñ¯¿¾ñ¯¿¿ñ¿¿¾ñ¿¿¿ò¿¾ò¿¿òŸ¿¾òŸ¿¿ò¯¿¾ò¯¿¿ò¿¿¾ò¿¿¿ó¿¾ó¿¿óŸ¿¾óŸ¿¿ó¯¿¾ó¯¿¿ó¿¿¾ó¿¿¿ô¿¾ô¿¿"
},
NULLTEST
};

static char*
trim(const char* s)
{
    int i;
    size_t l = strlen(s);
    char* t = strdup(s);
    for(i=l-1;i >= 0; i--) {
        if(t[i] != ' ') break;
    }
    t[i+1] = '\0';
    return t;
}

static int
test(const struct Test* tests, const char* title)
{
    int status = NC_NOERR;
    int failures = 0;
    const struct Test* p;

    fprintf(stderr,"Testing %s...\n",title);
    for(p=tests;p->id;p++) {
	char* id;
        char* description;
        const char* pf;
        id = trim(p->id);
        description = trim(p->description);
        status = nc_utf8_validate((const unsigned char*)p->data);
        if(status == NC_NOERR && p->xfail) {pf = "Fail"; failures++;}
        else if(status != NC_NOERR && p->xfail) pf = "Pass";
        else if(status == NC_NOERR && !p->xfail) pf = "Pass";
        else if(status != NC_NOERR && !p->xfail) {pf = "Fail"; failures++;}
        fprintf(stderr,"%s: %s %s\n",pf,id,description);
        fflush(stderr);
	free(id);
	free(description);
    }
    return failures;
}

int
main(int argc, char** argv)
{
    int failures = 0;

    printf("\n Testing UTF-8 sequences.\n");
    failures += test(utf8ok,"Correct Sequences");
    failures += test(utf8boundary,"Boundary Tests");
    failures += test(utf8bad,"Invalid strings");
    failures += test(utf8problematic,"Problematic strings");
    failures += test(utf8nonchars,"Other non-characters");
    fprintf(stderr,"No. of failures = %d\n",failures);
    exit(failures == 0 ? 0 : 1);
}
