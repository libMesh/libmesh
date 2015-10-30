/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of Moc_PARMETIS_PartGraphKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c,v 1.1 2003-06-24 05:33:51 benkirk Exp $
 *
 */

#include <parmetis.h>


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void ParMETIS_V3_PartKway(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
  idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, float *tpwgts,
  float *ubvec, int *options, int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int h, i;
  int nvtxs = -1, npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;
  float avg, maximb, *mytpwgts;
  int moptions[10];
  int seed, dbglvl = 0;
  int iwgtflag, inumflag, incon, inparts, ioptions[10];
  float *itpwgts, iubvec[MAXNCON];

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

/*
dbglvl = 7;
printf("[%d]    wgtflag: %d, numflag: %d, ncon: %d, nparts: %d, tpwgts: %.3f %.3f %.3f, ubvec: %.3f, options: %d %d %d\n", mype, *wgtflag, *numflag, *ncon, *nparts, tpwgts[0], tpwgts[1], tpwgts[2], ubvec[0], options[0], options[1], options[2]);
*/

  /********************************/
  /* Try and take care bad inputs */
  /********************************/
  if (options != NULL && options[0] == 1)
    dbglvl = options[PMV3_OPTION_DBGLVL];
  CheckInputs(STATIC_PARTITION, npes, dbglvl, wgtflag, &iwgtflag, numflag, &inumflag, ncon, &incon,
  nparts, &inparts, tpwgts, &itpwgts, ubvec, iubvec, NULL, NULL, NULL, NULL, NULL, NULL,
  options, ioptions, part, comm);

/*
printf("[%d]    wgtflag: %d, numflag: %d, ncon: %d, nparts: %d, tpwgts: %.3f %.3f %.3f, ubvec: %.3f, options: %d %d %d\n", mype, iwgtflag, inumflag, incon, inparts, itpwgts[0], itpwgts[1], itpwgts[2], iubvec[0], ioptions[0], ioptions[1], ioptions[2]);
*/
  /*********************************/
  /* Take care the nparts = 1 case */
  /*********************************/
  if (inparts <= 1) {
    idxset(vtxdist[mype+1]-vtxdist[mype], 0, part); 
    *edgecut = 0;
    return;
  }

  /******************************/
  /* Take care of npes = 1 case */
  /******************************/
  if (npes == 1 && inparts > 1) {
    moptions[0] = 0;
    nvtxs = vtxdist[1];

    if (incon == 1) {
      METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt,
      &iwgtflag, &inumflag, &inparts, moptions, edgecut, part);
    }
    else {
      /* ADD: this is because METIS does not support tpwgts for all constraints */
      mytpwgts = fmalloc(inparts, "mytpwgts");
      for (i=0; i<inparts; i++)
        mytpwgts[i] = itpwgts[i*incon];

      METIS_mCPartGraphRecursive2(&nvtxs, &incon, xadj, adjncy, vwgt, adjwgt,
      &iwgtflag, &inumflag, &inparts, mytpwgts, moptions, edgecut, part);

      free(mytpwgts);

/*
      METIS_mCPartGraphKway(&nvtxs, &incon, xadj, adjncy, vwgt, adjwgt,
      &iwgtflag, &inumflag, &inparts, iubvec, moptions, edgecut, part);
*/
    }
 
    return;
  }


  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  /*****************************/
  /* Set up control structures */
  /*****************************/
  if (ioptions[0] == 1) {
    dbglvl = ioptions[PMV3_OPTION_DBGLVL];
    seed = ioptions[PMV3_OPTION_SEED];
  }
  else {
    dbglvl = GLOBAL_DBGLVL;
    seed = GLOBAL_SEED;
  }
  SetUpCtrl(&ctrl, inparts, &dbglvl, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*incon*amax(npes, inparts));
  ctrl.seed = (seed == 0) ? mype : seed*mype;
  ctrl.sync = GlobalSEMax(&ctrl, seed);
  ctrl.partType = STATIC_PARTITION;
  ctrl.ps_relation = -1;
  ctrl.tpwgts = itpwgts;
  scopy(incon, iubvec, ctrl.ubvec);

  graph = Moc_SetUpGraph(&ctrl, incon, vtxdist, xadj, vwgt, adjncy, adjwgt, &iwgtflag);

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  /*******************************************/
  /* Check for funny cases                   */
  /* 	-graph with no edges                 */
  /* 	-graph with self edges               */
  /* 	-graph with poor vertex distribution */
  /* 	-graph with less than 2*npe nodes    */
  /*******************************************/
  if (vtxdist[npes] < SMALLGRAPH || vtxdist[npes] < npes*20) {
    IFSET(ctrl.dbglvl, DBG_INFO,
    rprintf(&ctrl, "Partitioning a graph of size %d serially\n", vtxdist[npes]));
    PartitionSmallGraph(&ctrl, graph, &wspace);
  }
  else {
    /***********************/
    /* Partition the graph */
    /***********************/
    Moc_Global_Partition(&ctrl, graph, &wspace);
    ParallelReMapGraph(&ctrl, graph, &wspace);
  }

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  /*******************/
  /* Print out stats */
  /*******************/
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_INFO,
     rprintf(&ctrl, "Final %d-way CUT: %6d \tBalance: ", inparts, graph->mincut));
  if (ctrl.dbglvl&DBG_INFO) {
    avg = 0.0;
    for (h=0; h<incon; h++) {
      maximb = 0.0;
      for (i=0; i<inparts; i++)
        maximb = amax(maximb, graph->gnpwgts[i*incon+h]/itpwgts[i*incon+h]);
      avg += maximb;
      rprintf(&ctrl, "%.3f ", maximb);
    }
  }
  IFSET(ctrl.dbglvl, DBG_INFO, rprintf(&ctrl, "  avg: %.3f\n", avg/(float)incon));

  GKfree((void **)&itpwgts, (void **)&graph->lnpwgts, (void **)&graph->gnpwgts, (void **)&graph->nvwgt, LTERM);
  FreeInitialGraphAndRemap(graph, iwgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (inumflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);

}

