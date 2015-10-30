/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh.c
 *
 * This file contains routines for constructing the dual graph of a mesh.
 * Assumes that each processor has at least one mesh element.
 *
 * Started 10/19/94
 * George
 *
 * $Id: mesh.c,v 1.1 2003-06-24 05:33:51 benkirk Exp $
 *
 */

#include <parmetis.h>
#define	MAXLINE	8192

/*************************************************************************
* This function converts a mesh into a dual graph
**************************************************************************/
void ParMETIS_V3_Mesh2Dual(idxtype *elmdist, idxtype *elements, int *etype, int *mgcnum,
  int *numflag, idxtype **xadj, idxtype **adjncy, MPI_Comm *comm)
{
  int i, j, jj, k, kk, kkk, m;
  int npes, mype, pe, count, mask, pass;
  int nelms, lnns, my_nns, node;
  int firstelm, lastelm, firstnode, lnode, element, nrecv, nsend;
  int *scounts, *rcounts, *sdispl, *rdispl;
  idxtype *nodedist, *nmap;
  idxtype *gnptr, *gnind, *nptr, *nind, *myxadj, *myadjncy = NULL;
  idxtype *sbuffer, *rbuffer, *htable;
  KeyValueType *nodelist, *recvbuffer;
  idxtype ind[200], wgt[200];
  int esize, esizes[5] = {-1, 3, 4, 8, 4};
  int maxnode, gmaxnode, minnode, gminnode;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);
  esize = esizes[*etype];

  if (*numflag == 1) 
    ChangeNumberingMesh(elmdist, elements, NULL, NULL, NULL, npes, mype, esize, 1);

  nelms = elmdist[mype+1]-elmdist[mype];
  mask = (1<<11)-1;

  /*****************************/
  /* Determine number of nodes */
  /*****************************/
  minnode = elements[idxamin(nelms*esize, elements)];
  MPI_Allreduce((void *)&minnode, (void *)&gminnode, 1, MPI_INT, MPI_MIN, *comm);
  for (i=0; i<nelms*esize; i++)
    elements[i] -= gminnode;

  maxnode = elements[idxamax(nelms*esize, elements)];
  MPI_Allreduce((void *)&maxnode, (void *)&gmaxnode, 1, MPI_INT, MPI_MAX, *comm);

  /**************************/
  /* Check for input errors */
  /**************************/
  ASSERTS(nelms > 0);

  nodedist = idxsmalloc(npes+1, 0, "nodedist");
  /* construct node distribution array */
  nodedist[0] = 0;
  for (i=0,j=gmaxnode+1; i<npes; i++) {
    k = j/(npes-i);
    nodedist[i+1] = nodedist[i]+k;
    j -= k;
  }
  my_nns = nodedist[mype+1]-nodedist[mype];
  firstnode = nodedist[mype];

  nodelist = (KeyValueType *)GKmalloc(nelms*esize*sizeof(KeyValueType), "nodelist");
  htable = idxsmalloc(amax(my_nns, mask+1), -1, "htable");
  scounts = imalloc(4*npes+2, "scounts");
  rcounts = scounts+npes;
  sdispl  = scounts+2*npes;
  rdispl  = scounts+3*npes+1;

  /*********************************************/
  /* first find a local numbering of the nodes */
  /*********************************************/
  for (i=0; i<nelms; i++)
    for (j=0; j<esize; j++) {
      nodelist[i*esize+j].key = elements[i*esize+j];
      nodelist[i*esize+j].val = i*esize+j;
    }

  ikeysort(nelms*esize, nodelist);

  count = 1;
  for (i=1; i<nelms*esize; i++)
    if (nodelist[i].key > nodelist[i-1].key)
      count++;

  lnns = count;
  nmap = idxmalloc(lnns, "nmap");

  /* renumber the nodes of the elements array */
  count = 1;
  nmap[0] = nodelist[0].key;
  elements[nodelist[0].val] = 0;
  for (i=1; i<nelms*esize; i++) {
    if (nodelist[i].key > nodelist[i-1].key) {
      nmap[count] = nodelist[i].key;
      count++;
    }
    elements[nodelist[i].val] = count-1;
  }
  MPI_Barrier(*comm);

  /**********************************************************/
  /* perform comms necessary to construct node-element list */
  /**********************************************************/
  pe = 0;
  iset(npes, 0, scounts);
  for (i=0; i<nelms*esize; i++) {
    while (nodelist[i].key >= nodedist[pe+1])
      pe++;

    scounts[pe] += 2;
  }
  ASSERTS(pe < npes);

  MPI_Alltoall((void *)scounts, 1, MPI_INT, (void *)rcounts, 1, MPI_INT, *comm);

  for (i=0; i<npes; i++) {
    sdispl[i] = scounts[i];
    rdispl[i] = rcounts[i];
  }

  MAKECSR(i, npes, sdispl);
  MAKECSR(i, npes, rdispl);
  ASSERTS(sdispl[npes] == nelms*esize*2);

  nrecv = rdispl[npes]/2;
  recvbuffer = (KeyValueType *)GKmalloc(amax(1, nrecv)*sizeof(KeyValueType), "recvbuffer");

  MPI_Alltoallv((void *)nodelist, scounts, sdispl, IDX_DATATYPE,
  (void *)recvbuffer, rcounts, rdispl, IDX_DATATYPE, *comm);

  /**************************************/
  /* construct global node-element list */
  /**************************************/
  gnptr = idxsmalloc(my_nns+1, 0, "gnptr");

  for (i=0; i<npes; i++) {
    for (j=rdispl[i]/2; j<rdispl[i+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      ASSERTS(lnode >= 0 && lnode < my_nns)

      gnptr[lnode]++;
    }
  }
  MAKECSR(i, my_nns, gnptr);

  gnind = idxmalloc(amax(1, gnptr[my_nns]), "gnind");
  for (pe=0; pe<npes; pe++) {
    firstelm = elmdist[pe];
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      element = recvbuffer[j].val/esize+firstelm;
      gnind[gnptr[lnode]++] = element;
    }
  }

  for (i=my_nns; i>0; i--)
    gnptr[i] = gnptr[i-1];
  gnptr[0] = 0;

  /*********************************************************/
  /* send the node-element info to the relevant processors */
  /*********************************************************/
  iset(npes, 0, scounts);

  /* use a hash table to ensure that each node is sent to a proc only once */
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      if (htable[lnode] == -1) {
        scounts[pe] += gnptr[lnode+1]-gnptr[lnode];
        htable[lnode] = 1;
      }
    }

    /* now reset the hash table */
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      htable[lnode] = -1;
    }
  }


  MPI_Alltoall((void *)scounts, 1, MPI_INT, (void *)rcounts, 1, MPI_INT, *comm);

  for (i=0; i<npes; i++) {
    sdispl[i] = scounts[i];
  }
  MAKECSR(i, npes, sdispl);

  /* create the send buffer */
  nsend = sdispl[npes];
  sbuffer = (idxtype *)realloc(nodelist, sizeof(idxtype)*amax(1, nsend));

  count = 0;
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      if (htable[lnode] == -1) {
        for (k=gnptr[lnode]; k<gnptr[lnode+1]; k++) {
          if (k == gnptr[lnode])
            sbuffer[count++] = -1*(gnind[k]+1);
          else
            sbuffer[count++] = gnind[k];
        }
        htable[lnode] = 1;
      }
    }
    ASSERTS(count == sdispl[pe+1]);

    /* now reset the hash table */
    for (j=rdispl[pe]/2; j<rdispl[pe+1]/2; j++) {
      lnode = recvbuffer[j].key-firstnode;
      htable[lnode] = -1;
    }
  }

  for (i=0; i<npes; i++) {
    rdispl[i] = rcounts[i];
  }
  MAKECSR(i, npes, rdispl);

  nrecv = rdispl[npes];
  rbuffer = (idxtype *)realloc(recvbuffer, sizeof(idxtype)*amax(1, nrecv));

  MPI_Alltoallv((void *)sbuffer, scounts, sdispl, IDX_DATATYPE,
  (void *)rbuffer, rcounts, rdispl, IDX_DATATYPE, *comm);

  k = -1;
  nptr = idxsmalloc(lnns+1, 0, "nptr");
  nind = rbuffer;
  for (pe=0; pe<npes; pe++) {
    for (j=rdispl[pe]; j<rdispl[pe+1]; j++) {
      if (nind[j] < 0) {
        k++;
        nind[j] = (-1*nind[j])-1;
      }
      nptr[k]++;
    }
  }
  MAKECSR(i, lnns, nptr);

  ASSERTS(k+1 == lnns);
  ASSERTS(nptr[lnns] == nrecv)

  myxadj = *xadj = idxsmalloc(nelms+1, 0, "xadj");
  idxset(mask+1, -1, htable);

  firstelm = elmdist[mype];
  lastelm = elmdist[mype+1]-1;

  /* Two passes -- in first pass, simply find out the memory requirements */
  for (pass=0; pass<2; pass++) {
    for (i=0; i<nelms; i++) {
      count = 0;
      for (j=0; j<esize; j++) {
        node = elements[esize*i+j];
        for (kkk=nptr[node+1]-1; kkk>=nptr[node]; kkk--) {
          kk = nind[kkk];

          if (kk == i+firstelm)
            continue;
          k = kk&mask;
          m = htable[k];

          if (m == -1) {
            ind[count] = kk;
            wgt[count] = 1;
            htable[k] = count++;
          }
          else {
            if (ind[m] == kk) {
              wgt[m]++;
            }
            else {
              for (jj=0; jj<count; jj++) {
                if (ind[jj] == kk) {
                  wgt[jj]++;
                  break;
                }
              }
              if (jj == count) {
                ind[count] = kk;
                wgt[count++] = 1;
              }
            }
          }
        }
      }
      for (j=0; j<count; j++) {
        if (wgt[j] >= *mgcnum) {
          if (pass == 0) {
            myxadj[i]++;
          }
          else {
            myadjncy[myxadj[i]++] = ind[j];
          }
        }
        htable[ind[j]&mask] = -1;
      }
    }

    if (pass == 0) {
      MAKECSR(i, nelms, myxadj);
      myadjncy = *adjncy = idxmalloc(myxadj[nelms], "adjncy");
    }
    else {
      for (i=nelms; i>0; i--)
        myxadj[i] = myxadj[i-1];
      myxadj[0] = 0;
    }
  }

  /*****************************************/
  /* correctly renumber the elements array */
  /*****************************************/
  for (i=0; i<nelms*esize; i++)
    elements[i] = nmap[elements[i]];

  if (*numflag == 1) 
    ChangeNumberingMesh(elmdist, elements, myxadj, myadjncy, NULL, npes, mype, esize, 0);

  /* do not free nodelist, recvbuffer, rbuffer */
  GKfree((void **)&scounts, (void **)&nodedist, (void **)&nmap, (void **)&sbuffer, (void **)&htable, (void **)&nptr, (void **)&nind, (void **)&gnptr, (void **)&gnind, LTERM);

  return;
}


