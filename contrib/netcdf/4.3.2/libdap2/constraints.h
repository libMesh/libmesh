 /*********************************************************************
  *   Copyright 1993, UCAR/Unidata
  *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
  *********************************************************************/
#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H 1

extern NCerror parsedapconstraints(NCDAPCOMMON*, char*, DCEconstraint*);
extern NCerror mapconstraints(DCEconstraint*,CDFnode*);
extern NCerror qualifyconstraints(DCEconstraint* constraint);
extern NCerror computeprojectedvars(NCDAPCOMMON*,DCEconstraint*);

extern char* simplepathstring(NClist* segments, char* separator);
extern void makesegmentstring(NClist* segments, NCbytes* buf, char* separator);

extern int iswholeslice(DCEslice*, struct CDFnode* dim);
extern int iswholesegment(DCEsegment*);

extern int iswholeconstraint(DCEconstraint* con);

extern char* buildprojectionstring(NClist* projections);
extern char* buildselectionstring(NClist* selections);
extern char* buildconstraintstring(DCEconstraint* constraints);

extern void makewholesegment(DCEsegment*,struct CDFnode*);
extern void makewholeslice(DCEslice* slice, struct CDFnode* dim);

extern NCerror fixprojections(NClist* list);

extern int dapvar2projection(CDFnode* var, DCEprojection** projectionp);
extern int daprestrictprojection(NClist* projections, DCEprojection* var, DCEprojection** resultp);
extern int dapshiftprojection(DCEprojection*);

#endif /*CONSTRAINTS_H*/
