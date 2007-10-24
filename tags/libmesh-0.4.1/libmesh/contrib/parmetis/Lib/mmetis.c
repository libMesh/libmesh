/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmetis.c
 *
 * This is the entry point of ParMETIS_V3_PartMeshKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: mmetis.c,v 1.1 2003-06-24 05:33:51 benkirk Exp $
 *
 */

#include <parmetis.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel mesh partitionioner. 
* This function assumes nothing about the mesh distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_V3_PartMeshKway(idxtype *elmdist, idxtype *elements, idxtype *elmwgt,
  int *etype, int *mgcnum, int *wgtflag, int *numflag, int *ncon, int *nparts,
  float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part,
  MPI_Comm *comm)
{
  int i;
  int npes, mype;
  int mywgtflag, mynumflag, gnedges;
  idxtype *xadj, *adjncy;
  GraphType *graph;
  timer TotalTmr, Mesh2DualTmr, ParMETISTmr;
  CtrlType ctrl;
  MeshType *mesh;
  int esize, esizes[5] = {-1, 3, 4, 8, 4};
  int dbglvl = 0;
  int iwgtflag, inumflag, incon, inparts, ioptions[10], ietype, imgcnum;
  float *itpwgts, iubvec[MAXNCON];


  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];
  CheckInputs(MESH_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag,
  ncon, &incon, nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, etype, &ietype,
  mgcnum, &imgcnum, NULL, NULL, options, ioptions, part, comm);

  cleartimer(TotalTmr);
  cleartimer(Mesh2DualTmr);
  cleartimer(ParMETISTmr);

  MPI_Barrier(*comm);
  starttimer(TotalTmr);
  starttimer(Mesh2DualTmr);

  esize = esizes[ietype];
  if (inumflag == 1)
    ChangeNumberingMesh(elmdist, elements, NULL, NULL, NULL, npes, mype, esize, 1);

  ctrl.npes = npes;
  ctrl.mype = mype;
  ctrl.comm = *comm;
  if (ioptions[0] == 1) 
    ctrl.dbglvl = ioptions[PMV3_OPTION_DBGLVL];
  else 
    ctrl.dbglvl = GLOBAL_DBGLVL;
  mesh = SetUpMesh(&ietype, &incon, elmdist, elements, elmwgt, &iwgtflag, comm);

  mynumflag = 0;
  ParMETIS_V3_Mesh2Dual(elmdist, elements, &ietype, &imgcnum, &mynumflag, &xadj, &adjncy, comm);
  graph = CreateGraph();
  graph->vtxdist = elmdist;
  graph->xadj = xadj;
  graph->adjncy = adjncy;
  graph->vwgt = mesh->elmwgt;
  graph->ncon = mesh->ncon;
  graph->nvtxs = elmdist[mype+1]-elmdist[mype];
  graph->nedges = xadj[graph->nvtxs];
  graph->gnvtxs = elmdist[npes];

  MPI_Allreduce((void *)&(graph->nedges), (void *)&gnedges, 1, MPI_INT, MPI_SUM, *comm);
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "Completed Dual Graph -- Nvtxs: %d, Nedges: %d\n", graph->gnvtxs, gnedges/2));

  stoptimer(Mesh2DualTmr);
  MPI_Barrier(*comm);

  /* Put the elements array back the way you found it */
  for (i=0; i<mesh->nelms*mesh->esize; i++)
    elements[i] += mesh->gminnode;

  /***********************/
  /* Partition the graph */
  /***********************/
  starttimer(ParMETISTmr);

  mynumflag = 0;
  mywgtflag = 2;
  ParMETIS_V3_PartKway(graph->vtxdist, graph->xadj, graph->adjncy, graph->vwgt, NULL,
  &mywgtflag, &mynumflag, &incon, &inparts, itpwgts, iubvec, ioptions, edgecut, part, comm);

  MPI_Barrier(*comm);
  stoptimer(ParMETISTmr);
  stoptimer(TotalTmr);

  if (inumflag == 1)
    ChangeNumberingMesh(elmdist, elements, graph->xadj, graph->adjncy, part, npes, mype, esize, 0);

  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, Mesh2DualTmr,	"   Mesh2Dual"));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, ParMETISTmr,	"    ParMETIS"));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimer(&ctrl, TotalTmr,		"       Total"));

  FreeGraph(graph);
  GKfree((void **)&mesh, (void **)&itpwgts, LTERM);
  return;
}

