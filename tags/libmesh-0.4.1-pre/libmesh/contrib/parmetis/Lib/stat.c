/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file computes various statistics
 *
 * Started 7/25/97
 * George
 *
 * $Id: stat.c,v 1.1 2003-06-24 05:33:51 benkirk Exp $
 *
 */

#include <parmetis.h>



/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void Moc_ComputeSerialBalance(CtrlType *ctrl, GraphType *graph, idxtype *where,
float *ubvec)
{
  int i, j, nvtxs, ncon, nparts;
  idxtype *pwgts, *tvwgts, *vwgt;
  float *tpwgts, maximb;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  vwgt = graph->vwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;

  pwgts = idxsmalloc(nparts*ncon, 0, "pwgts");
  tvwgts = idxsmalloc(ncon, 0, "tvwgts");

  for (i=0; i<graph->nvtxs; i++)
    for (j=0; j<ncon; j++) {
      pwgts[where[i]*ncon+j] += vwgt[i*ncon+j];
      tvwgts[j] += vwgt[i*ncon+j];
    }

  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb = amax(maximb, (float)pwgts[i*ncon+j]/(tpwgts[i*ncon+j]*(float)tvwgts[j]));
    ubvec[j] = maximb;
  }

  GKfree((void **)&pwgts, (void **)&tvwgts, LTERM);

}

/*************************************************************************
* This function computes the balance of the partitioning
**************************************************************************/
void Moc_ComputeParallelBalance(CtrlType *ctrl, GraphType *graph, idxtype *where,
float *ubvec)
{
  int i, j, nvtxs, ncon, nparts;
  float *nvwgt, *lnpwgts, *gnpwgts;
  float *tpwgts, maximb;
  MPI_Comm comm;

  ncon = graph->ncon;
  nvtxs = graph->nvtxs;
  nvwgt = graph->nvwgt;
  nparts = ctrl->nparts;
  tpwgts = ctrl->tpwgts;
  comm = ctrl->comm;

  lnpwgts = fmalloc(nparts*ncon, "CPB: lnpwgts");
  gnpwgts = fmalloc(nparts*ncon, "CPB: gnpwgts");
  sset(nparts*ncon, 0.0, lnpwgts);

  for (i=0; i<nvtxs; i++)
    for (j=0; j<ncon; j++)
      lnpwgts[where[i]*ncon+j] += nvwgt[i*ncon+j];

  MPI_Allreduce((void *)(lnpwgts), (void *)(gnpwgts), nparts*ncon,
      MPI_FLOAT, MPI_SUM, comm);

  for (j=0; j<ncon; j++) {
    maximb = 0.0;
    for (i=0; i<nparts; i++)
      maximb = amax(maximb, gnpwgts[i*ncon+j]/tpwgts[i*ncon+j]);
    ubvec[j] = maximb;
  }

  GKfree((void **)&lnpwgts, (void **)&gnpwgts, LTERM);
  return;
}


/*************************************************************************
* This function prints a matrix
**************************************************************************/
void Moc_PrintThrottleMatrix(CtrlType *ctrl, GraphType *graph, float *matrix)
{
  int i, j;

  for (i=0; i<ctrl->npes; i++) {
    if (i == ctrl->mype) {
      for (j=0; j<ctrl->npes; j++)
        printf("%.3f ", matrix[j]);
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(ctrl->comm);
  }

  if (ctrl->mype == 0) {
    printf("****************************\n");
    fflush(stdout);
  }
  MPI_Barrier(ctrl->comm);

  return;
}


/*************************************************************************
*  This function computes stats for refinement
**************************************************************************/
void Moc_ComputeRefineStats(CtrlType *ctrl, GraphType *graph, float *ubvec)
{
  int h, i, j, k;
  int nvtxs, ncon;
  idxtype *xadj, *adjncy, *adjwgt, *where;
  float *nvwgt, *lnpwgts, *gnpwgts;
  RInfoType *rinfo;
  int mype = ctrl->mype, nparts = ctrl->nparts;
  idxtype *gborder, *border, *gfrom, *from, *gto, *to, *connect, *gconnect;
  idxtype gain[20] = {0}, ggain[20];
  int lnborders, gnborders;
  int bestgain, pmoves, gpmoves, other;
  float tpwgts[MAXNCON], badmaxpwgt[MAXNCON];
  int HIST_FACTOR = graph->level + 1;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  lnpwgts = graph->lnpwgts;
  gnpwgts = graph->gnpwgts;
  rinfo = graph->rinfo;

  connect = idxsmalloc(nparts*nparts, 0, "CRS: connect");
  gconnect = idxmalloc(nparts*nparts, "CRS: gconnect");
  border = idxsmalloc(nparts, 0, "CRS: border");
  gborder = idxmalloc(nparts, "CRS: gborder");
  from = idxsmalloc(nparts, 0, "CRS: from");
  gfrom = idxmalloc(nparts, "CRS: gfrom");
  to = idxsmalloc(nparts, 0, "CRS: to");
  gto = idxmalloc(nparts, "CRS: gto");

  for (h=0; h<ncon; h++) {
    tpwgts[h] = ssum_strd(nparts, gnpwgts+h, ncon)/(float)(nparts);
    badmaxpwgt[h] = ubvec[h]*tpwgts[h];
  }

  if (mype == 0) printf("******************************\n");
  if (mype == 0) printf("******************************\n");

  /***************************************/
  if (mype == 0) {
    printf("subdomain weights:\n");
    for (h=0; h<ncon; h++) {
      for (i=0; i<nparts; i++)
        printf("%9.3f ", gnpwgts[i*ncon+h]);
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  if (mype == 0) {
    printf("subdomain imbalance:\n");
    for (h=0; h<ncon; h++) {
      for (i=0; i<nparts; i++)
        printf("%9.3f ", gnpwgts[i*ncon+h] * (float)(nparts));
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  for (i=0; i<nparts; i++)
    connect[i*nparts+i] = -1;

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        connect[where[i]*nparts+where[adjncy[j]]] = 1;
        connect[where[adjncy[j]]*nparts+where[i]] = 1;
      }
    }
  }

  MPI_Reduce((void *)connect, (void *)gconnect, nparts*nparts, IDX_DATATYPE, MPI_MAX, 0, ctrl->comm);
  if (mype == 0) { 
    printf("connectivity\n");
    for (i=0; i<nparts; i++) {
      printf("%d: ", i);
      for (j=0; j<nparts; j++)
        printf("%9d ", gconnect[i*nparts+j]);
      printf("\n");
    }
    printf("\n");
  }

  /***************************************/
  lnborders = 0;
  for (i=0; i<nvtxs; i++)
    if (rinfo[i].ndegrees > 0) {
      lnborders++;
      border[where[i]]++;
    } 

  MPI_Reduce((void *)border, (void *)gborder, nparts, IDX_DATATYPE, MPI_SUM, 0, ctrl->comm);
  gnborders = GlobalSESum(ctrl, lnborders);
  if (mype == 0) {
    printf("number of borders: %d\n", gnborders);
    for (i=0; i<nparts; i++)
      printf("%9d ", gborder[i]);
    printf("\n\n");
  }

  /***************************************/
  pmoves = 0;
  for (i=0; i<nvtxs; i++) {
    nvwgt = graph->nvwgt+i*ncon;

    for (j=0; j<rinfo[i].ndegrees; j++) {
      other = rinfo[i].degrees[j].edge;
      for (h=0; h<ncon; h++)
        if (gnpwgts[other*ncon+h]+nvwgt[h] > badmaxpwgt[h])
          break;

      if (h == ncon)
        break;
    }

    if (j < rinfo[i].ndegrees) {
      pmoves++;
      from[where[i]]++;
      to[other]++;
      for (k=j+1; k<rinfo[i].ndegrees; k++) {
        other = rinfo[i].degrees[k].edge;
        for (h=0; h<ncon; h++)
          if (gnpwgts[other*ncon+h]+nvwgt[h] > badmaxpwgt[h])
            break;

        if (h == ncon) {
          pmoves++;
          from[where[i]]++;
          to[other]++;
        }
      }
    }
  }

  gpmoves = GlobalSESum(ctrl, pmoves);
  MPI_Reduce((void *)from, (void *)gfrom, nparts, IDX_DATATYPE, MPI_SUM, 0, ctrl->comm);
  MPI_Reduce((void *)to, (void *)gto, nparts, IDX_DATATYPE, MPI_SUM, 0, ctrl->comm);

  if (mype == 0) {
    printf("possible moves: %d\n", gpmoves);
    printf("from   ");
    for (i=0; i<nparts; i++) {
      printf("%9d ", gfrom[i]);
    }
    printf("\n");
    printf("to     ");
    for (i=0; i<nparts; i++) {
      printf("%9d ", gto[i]);
    }
    printf("\n\n");
  }

  /***************************************/
  for (i=0; i<nvtxs; i++) {
    if (rinfo[i].ndegrees > 0) {
      bestgain = rinfo[i].degrees[0].ewgt-rinfo[i].id;
      for (j=0; j<rinfo[i].ndegrees; j++)
        bestgain = amax(bestgain, rinfo[i].degrees[j].ewgt-rinfo[i].id);

      if (bestgain / HIST_FACTOR >= 10) {
        gain[19]++;
        continue;
      }

      if (bestgain / HIST_FACTOR < -10) {
        gain[0]++;
        continue;
      }

      gain[(bestgain/HIST_FACTOR)+10]++;
    }
  }

  MPI_Reduce((void *)gain, (void *)ggain, 20, IDX_DATATYPE, MPI_SUM, 0, ctrl->comm);
  if (mype == 0) {
    printf("gain histogram (buckets of %d)\n", HIST_FACTOR);
    for (i=0; i<20; i++) {
      if (i == 10 || i == 11)
        printf("    ");
      printf("%d ", ggain[i]);
    }
    printf("\n\n");
  }




  /***************************************/
  if (mype == 0) printf("******************************\n");
  if (mype == 0) printf("******************************\n");

  GKfree((void **)&gconnect, (void **)&connect, (void **)&gborder, (void **)&border, (void **)&gfrom, (void **)&from, (void **)&gto, (void **)&to, LTERM);
  return;
}
