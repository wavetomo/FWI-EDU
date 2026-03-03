/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-05-07 14:24:54
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-11-03 10:49:49
 * @FilePath: /FWILab2dCLangV1.2.3/src/checkpoint.c
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#include "checkpoint.h"

void acoustic2dCheckPoint_init(Acoustic2dCheckPoint *chkpt)
{ 
  /*
  acoustic2dCheckPoint_init: initialize all members of an Acoustic2dCheckPoint structure to NULL.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure

  Output:
      chkpt   - structure whose internal pointers are all initialized to NULL
  */
  chkpt->status = NULL;
  chkpt->p = NULL;
  chkpt->storage_p = NULL;
}

void alloc_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, int NZ, int NX, int nt)
{ 
  /*
  alloc_acoustic2dCheckPoint: allocate memory for storing 2D acoustic wavefield data.

  Input:
      NX      - number of samples along the x-direction of the domain 
                (NX = nx + 2 * pmlThick)
      NZ      - number of samples along the z-direction of the domain 
                (NZ = nz + 2 * pmlThick)
      nt      - number of time samples

  Output:
      chkpt   - an Acoustic2dCheckPoint structure containing pointers 
                to allocated wavefield data arrays
  */
  int it, ichk, icur;
  chkpt->nt = nt;
  chkpt->dCheckPoint = MAX(roundf(sqrtf(2.0f * (nt - 1))), 3);
  /* status */
  if (chkpt->status != NULL)
  {
    free1int(chkpt->status);
  }
  chkpt->status = alloc1int(nt);
  memset(chkpt->status, 0, nt * ISIZE);
  /* p */
  if (chkpt->p != NULL)
  {
    free(chkpt->p);
  }
  chkpt->p = (float ***)malloc(nt * sizeof(float **));
  for (it = 0; it < nt; it++)
  {
    chkpt->p[it] = NULL;
  }
  chkpt->nCheckPoint = 0;
  for (it = 0; it < nt; it++)
  { /* second-order equation need 2 wavefields to reconstruct! */
    if (it % chkpt->dCheckPoint == 0 || (it + 1) % chkpt->dCheckPoint == 0)
    {
      chkpt->status[it] = 1;
      chkpt->nCheckPoint++;
    }
  }
  /* storage */
  chkpt->ntPoint = chkpt->nCheckPoint + chkpt->dCheckPoint - 2;
  chkpt->memCost = 1.0f * (chkpt->ntPoint * NZ * NX * FSIZE) / powf(1024.0f, 3);
  if (chkpt->storage_p != NULL)
  {
    free3float(chkpt->storage_p);
  }
  chkpt->storage_p = alloc3float(NZ, NX, chkpt->ntPoint);
  /***** *****/
  ichk = 0;
  for (it = 0; it < nt; it++)
  {
    if (chkpt->status[it] == 1)
    {
      chkpt->p[it] = chkpt->storage_p[ichk];
      ichk++;
      icur = chkpt->nCheckPoint;
    }
    else
    {
      chkpt->p[it] = chkpt->storage_p[icur];
      icur++;
    }
  }
  memset(chkpt->status, 0, nt * ISIZE);
}

void free_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt)
{
  /*
  free_acoustic2dCheckPoint: free the memory allocated for 2D acoustic wavefield data.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                pointers to allocated wavefield data arrays

  Output:
      chkpt   - structure whose internal pointers are reset to NULL
  */
  if (chkpt->status != NULL)
  {
    free1int(chkpt->status);
    chkpt->status = NULL;
  }
  if (chkpt->p != NULL)
  {
    free(chkpt->p);
    chkpt->p = NULL;
  }
  if (chkpt->storage_p != NULL)
  {
    free3float(chkpt->storage_p);
    chkpt->storage_p = NULL;
  }
}

void save_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, float **p0, float **p1, int NP, int it)
{ 
  /*
  save_acoustic2dCheckPoint: save the checkpoint data of a 2D acoustic wavefield.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                pointers to allocated wavefield data arrays
      p0      - 2D array of the acoustic wavefield at the previous time step
      p1      - 2D array of the acoustic wavefield at the current time step
      NP      - number of grid points in the 2D acoustic wavefield (NX x NZ)
      it      - current time index of the simulation

  Output:
      chkpt   - updated checkpoint structure with stored wavefield data
  */
  if (it >= 0 && it < chkpt->nt)
  {
    if (it % chkpt->dCheckPoint == 0)
    {
      if (it - 1 >= 0)
      {
        memcpy(chkpt->p[it - 1][0], p0[0], NP * FSIZE);
        chkpt->status[it - 1] = 1;
      }
      memcpy(chkpt->p[it][0], p1[0], NP * FSIZE);
      chkpt->status[it] = 1;
    }
  }
}

void get_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, float **p0, float **p1, int NP, int it)
{ 
  /*
  get_acoustic2dCheckPoint: retrieve the checkpoint data of a 2D acoustic wavefield.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                stored wavefield data arrays
      NP      - number of grid points in the 2D acoustic wavefield (NX x NZ)
      it      - current time index of the simulation

  Output:
      p0      - 2D array to store the acoustic wavefield at the previous time step
      p1      - 2D array to store the acoustic wavefield at the current time step
  */
  if (it >= 0 && it < chkpt->nt)
  {
    if (it % chkpt->dCheckPoint == 0)
    {
      if (it - 1 >= 0)
      {
        memcpy(p0[0], chkpt->p[it - 1][0], NP * FSIZE);
      }
      else
      {
        memset(p0[0], 0, NP * FSIZE);
      }
      memcpy(p1[0], chkpt->p[it][0], NP * FSIZE);
    }
  }
}

void save_acoustic2dWavefield(Acoustic2dCheckPoint *chkpt, float **p, int NP, int it)
{ 
  /*
  save_acoustic2dWavefield: save the 2D acoustic wavefield at a given time step.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                pointers to allocated wavefield data arrays
      p       - 2D array of the acoustic wavefield at the current time step
      NP      - number of grid points in the 2D acoustic wavefield
      it      - current time index of the simulation

  Output:
      chkpt   - updated checkpoint structure with stored wavefield data
  */
  if (it >= 0 && it < chkpt->nt)
  {
    memcpy(chkpt->p[it][0], p[0], NP * FSIZE);
    chkpt->status[it] = 1;
  }
}

void get_acoustic2dWavefield(Acoustic2dCheckPoint *chkpt, float **p, int NP, int it)
{ 
  /*
    get_acoustic2dWavefield: retrieve the 2D acoustic wavefield from a checkpoint.

    Input:
        chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                  pointers to stored wavefield data
        p       - 2D array to store the retrieved acoustic wavefield
        NP      - number of grid points in the 2D acoustic wavefield
        it      - time index of the stored wavefield to retrieve

    Output:
        p       - updated 2D array containing the retrieved wavefield data
  */
  if (it >= 0 && it < chkpt->nt)
  {
    memcpy(p[0], chkpt->p[it][0], NP * FSIZE);
  }
}

int get_acoustic2dCheckPointStatus(Acoustic2dCheckPoint *chkpt, int it)
{ 
  /*
  get_acoustic2dCheckPointStatus: check the checkpoint status of a 2D acoustic wavefield.

  Input:
      chkpt   - pointer to an Acoustic2dCheckPoint structure containing
                pointers to stored wavefield data and status flags
      it      - time index of the 2D acoustic wavefield to check

  Output:
      return  - if the checkpoint at time index 'it' exists, returns -1;
                otherwise, returns the index of the earliest time step
                that needs to be reconstructed
  */
  int itStart;
  if (chkpt->status[it] == 1)
  {
    itStart = -1;
  }
  else
  {
    itStart = it;
    while (itStart >= 0 && chkpt->status[itStart] == 0)
    {
      itStart--;
    }
  }
  return itStart;
}

float calculate_memory_cost_checkpoint(int nt, int NX, int NZ)
{
  /*
  calculate_memory_cost_checkpoint: estimate the memory cost of a checkpointing strategy 
  used in 2D acoustic wavefield simulations.

  Input:
      nt  - total number of time steps
      NX  - number of grid points along the x-direction
      NZ  - number of grid points along the z-direction

  Output:
      return - estimated memory cost (in GB) required for storing checkpointed
                wavefields during forward modeling
  */
  int it;
  int dCheckPoint, nCheckPoint, ntPoint;
  float memCost;

  nCheckPoint = 0;
  dCheckPoint = MAX(roundf(sqrtf(2.0f * (nt - 1))), 3);
  for (it = 0; it < nt; it++)
  { /* second-order equation need 2 wavefields to reconstruct! */
    if (it % dCheckPoint == 0 || (it + 1) % dCheckPoint == 0)
    {
      nCheckPoint++;
    }
  }

  ntPoint = nCheckPoint + dCheckPoint - 2;
  memCost = 1.0f * (ntPoint * NZ * NX * FSIZE) / powf(1024.0f, 3);
  return memCost;
}