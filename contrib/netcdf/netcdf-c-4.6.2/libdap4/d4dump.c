/* Copyright 2012, UCAR/Unidata.
   See the LICENSE file for more information.
*/

#include "d4includes.h"

/*
Provide a simple dump of binary data
*/

/**************************************************/

void
NCD4_dumpbytes(size_t size, const void* data0, int swap)
{
    size_t extended;
    void* data;
    char* pos;
    size_t i;

    extended = size + 8;
    data = calloc(1,extended); /* provide some space to simplify the code */
    memcpy(data,data0,size);
    for(i=0,pos=data;i<size; pos++,i++) {
	struct {
	    unsigned char u8[1];
	      signed char i8[1];
	    unsigned short u16[1];
	             short i16[1];
	    unsigned int u32[1];
	             int i32[1];
	    unsigned long long u64[1];
	             long long i64[1];
	    float f32[1];
	    double f64[1];
	    char s[8];
	} v;
	v.s[0] = *((char*)pos);
	v.s[1] = '\0';
	v.u8[0] = *((unsigned char*)pos);
	v.i8[0] = *((signed char*)pos);
        v.u16[0] = *((unsigned short*)pos);
        v.i16[0] = *((short*)pos);
        v.u32[0] = *((unsigned int*)pos);
        v.i32[0] = *((int*)pos);
        v.u64[0] = *((unsigned long long*)pos);
        v.i64[0] = *((long long*)pos);
	if(swap) {
	    swapinline16(v.u16);
	    swapinline32(v.u32);
	    swapinline32(v.u64);
	    swapinline16(v.i16);
	    swapinline32(v.i32);
	    swapinline32(v.i64);
	    swapinline32(v.f32);
	    swapinline32(v.f64);
        }
        if(v.s[0] == '\r') strcpy(v.s,"\\r");
        else if(v.s[0] == '\n') strcpy(v.s,"\\n");
        else if(v.s[0] < ' ' || v.s[0] >= 0x7f) v.s[0] = '?';
        fprintf(stderr,"[%03lu] %02x %03u %4d", (unsigned long)i, v.u8[0], v.u8[0], v.i8[0]);
        fprintf(stderr," 0x%08x %12u %13d", v.u32[0], v.u32[0], v.i32[0]);
        fprintf(stderr," 0x%04x %06u %7d", v.u16[0], v.u16[0], v.i16[0]);
        fprintf(stderr," '%s'\n",v.s);
	fflush(stderr);
    }
}

void
NCD4_tagdump(size_t size, const void* data0, int swap, const char* tag)
{
    fprintf(stderr,"++++++++++ %s ++++++++++\n",tag);
    NCD4_dumpbytes(size,data0,swap);
    fprintf(stderr,"++++++++++ %s ++++++++++\n",tag);
    fflush(stderr);
}

/* Dump the variables in a group */
void
NCD4_dumpvars(NCD4node* group)
{
    int i;
    fprintf(stderr,"%s.vars:\n",group->name);
    for(i=0;i<nclistlength(group->vars);i++) {
	NCD4node* var = (NCD4node*)nclistget(group->vars,i);
	NCD4node* type;

	switch (var->subsort) {
	default:
	    type = var->basetype;
	    fprintf(stderr,"<%s name=\"%s\"/>\n",type->name,var->name);
	    break;
	case NC_STRUCT:
	    fprintf(stderr,"<%s name=\"%s\"/>\n","Struct",var->name);
	    break;
	case NC_SEQ:
	    fprintf(stderr,"<%s name=\"%s\"/>\n","Sequence",var->name);
	    break;
        }
    }
    fflush(stderr);
}

union ATOMICS*
NCD4_dumpatomic(NCD4node* var, void* data)
{
    union ATOMICS* p = (union ATOMICS*)data;
    return p;
}

