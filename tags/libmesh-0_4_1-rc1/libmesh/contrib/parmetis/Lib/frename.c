/*
 * frename.c
 *
 * This file contains some renaming routines to deal with different
 * Fortran compilers.
 *
 * Started 6/1/98
 * George
 *
 * $Id: frename.c,v 1.1 2003-06-24 05:33:51 benkirk Exp $
 *
 */

#include <parmetis.h>



void PARMETIS_PARTKWAY(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, comm);
}
void parmetis_partkway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, comm);
}
void parmetis_partkway_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, comm);
}
void parmetis_partkway__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part, comm);
}


void PARMETIS_REFINEKWAY(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RefineKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_refinekway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RefineKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_refinekway_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RefineKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_refinekway__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RefineKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}


void PARMETIS_PARTGEOMKWAY(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, nparts, options, edgecut, part, comm);
}
void parmetis_partgeomkway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, nparts, options, edgecut, part, comm);
}
void parmetis_partgeomkway_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, nparts, options, edgecut, part, comm);
}
void parmetis_partgeomkway__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, nparts, options, edgecut, part, comm);
}


void PARMETIS_PARTGEOMREFINE(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomRefine(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, options, edgecut, part, comm);
}
void parmetis_partgeomrefine(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomRefine(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, options, edgecut, part, comm);
}
void parmetis_partgeomrefine_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomRefine(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, options, edgecut, part, comm);
}
void parmetis_partgeomrefine__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeomRefine(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, ndims, xyz, options, edgecut, part, comm);
}


void PARMETIS_PARTGEOM(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeom(vtxdist, ndims, xyz, part, comm);
}
void parmetis_partgeom(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeom(vtxdist, ndims, xyz, part, comm);
}
void parmetis_partgeom_(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeom(vtxdist, ndims, xyz, part, comm);
}
void parmetis_partgeom__(idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_PartGeom(vtxdist, ndims, xyz, part, comm);
}


void PARMETIS_REPARTLDIFFUSION(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartLDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartldiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartLDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartldiffusion_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartLDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartldiffusion__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartLDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}

void PARMETIS_REPARTGDIFFUSION(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartgdiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartgdiffusion_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartgdiffusion__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartGDiffusion(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}

void PARMETIS_REPARTREMAP(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartremap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartremap_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartremap__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}


void PARMETIS_REPARTMLREMAP(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartMLRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartmlremap(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartMLRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartmlremap_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartMLRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}
void parmetis_repartmlremap__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  ParMETIS_RepartMLRemap(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, options, edgecut, part, comm);
}


void PARMETIS_NODEND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_NodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_nodend(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_NodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_nodend_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_NodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_nodend__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_NodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}


void PARMETIS_SERIALNODEND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_serialnodend(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_serialnodend_(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}
void parmetis_serialnodend__(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)
{
  ParMETIS_SerialNodeND(vtxdist, xadj, adjncy, numflag, options, order, sizes, comm);
}

