/* Copyright 2009, UCAR/Unidata and OPeNDAP, Inc.
   See the COPYRIGHT file for more information. */

#ifndef OCURI_H
#define OCURI_H

/*! This is an open structure meaning
	it is ok to directly access its fields*/
typedef struct OCURI {
    char* uri;        /* as passed by the caller */
    char* params;     /* all params */
    char** paramlist;    /*!<null terminated list */
    char* constraint; /*!< projection+selection */
    char* projection; /*!< without leading '?'*/
    char* selection;  /*!< with leading '&'*/
    char* strings;    /* first char of strings is always '\0' */
    /* Following all point into the strings field */
    char* protocol;
    char* user; /* from user:password@ */
    char* password; /* from user:password@ */
    char* host;	      /*!< host*/
    char* port;	      /*!< host */
    char* file;	      /*!< file */
} OCURI;

extern int ocuriparse(const char* s, OCURI** ocuri);
extern void ocurifree(OCURI* ocuri);

/* Replace the constraints */
extern void ocurisetconstraints(OCURI*,const char* constraints);

/* Construct a complete OC URI; caller frees returned string */

/* Define flags to control what is included */
#define OCURICONSTRAINTS	 1
#define OCURIUSERPWD	  	 2
#define OCURIPREFIXPARAMS  	 4
#define OCURISUFFIXPARAMS	 8
#define OCURIPARAMS	  	OCURIPREFIXPARAMS
#define OCURIENCODE		16 /* If output should be encoded */
#define OCURISTD	  	(OCURICONSTRAINTS|OCURIUSERPWD)

extern char* ocuribuild(OCURI*,const char* prefix, const char* suffix, int flags);


/* Param Management */
extern int ocuridecodeparams(OCURI* ocuri);
extern int ocurisetparams(OCURI* ocuri,const char*);

/*! 0 result => entry not found; 1=>found; result holds value (may be null).
    In any case, the result is imutable and should not be free'd.
*/
extern int ocurilookup(OCURI*, const char* param, const char** result);

extern char* ocuriencode(char* s, char* allowable);
extern char* ocuridecode(char* s);
extern char* ocuridecodeonly(char* s, char*);

#endif /*OCURI_H*/
