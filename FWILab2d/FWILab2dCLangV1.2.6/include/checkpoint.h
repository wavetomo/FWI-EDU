/*
 * @description:
 * @Author: PingMin Cheung
 * @Date: 2024-05-07 14:24:54
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-11-03 10:49:21
 * @FilePath: /FWILab2dCLangV1.2.3/include/checkpoint.h
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#pragma once

#include <su.h>

/********** Types **********/
typedef struct
{
  int nt;
  int dCheckPoint;
  int nCheckPoint;
  int ntPoint;
  int *status;
  float memCost;
  float ***p;         /* length: nt */
  float ***storage_p; /* length: ntPoint */
} Acoustic2dCheckPoint;


void acoustic2dCheckPoint_init(Acoustic2dCheckPoint *chkpt);

/*
acoustic2dCheckPoint_init: initialize all members of an Acoustic2dCheckPoint structure to NULL.

Input:
   chkpt   - pointer to an Acoustic2dCheckPoint structure

Output:
   chkpt   - structure whose internal pointers are all initialized to NULL
*/

void alloc_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, int NZ, int NX, int nt);

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

void free_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt);

/*
free_acoustic2dCheckPoint: free the memory allocated for 2D acoustic wavefield data.

Input:
   chkpt   - pointer to an Acoustic2dCheckPoint structure containing
               pointers to allocated wavefield data arrays

Output:
   chkpt   - structure whose internal pointers are reset to NULL
*/

void save_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, float **p0, float **p1, int NP, int it);

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

void get_acoustic2dCheckPoint(Acoustic2dCheckPoint *chkpt, float **p0, float **p1, int NP, int it);

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

void get_acoustic2dWavefield(Acoustic2dCheckPoint *chkpt, float **p, int NP, int it);

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

void save_acoustic2dWavefield(Acoustic2dCheckPoint *chkpt, float **p, int NP, int it);

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

int get_acoustic2dCheckPointStatus(Acoustic2dCheckPoint *chkpt, int it);

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

float calculate_memory_cost_checkpoint(int nt, int NX, int NZ);

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