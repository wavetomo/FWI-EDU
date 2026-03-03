/*
 * @Author: QuYang Chen
 * @Date: 2025-09-17 04:02:00
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-11-21 15:05:13
 * @FilePath: /FWILab2dCLangV1.2.4/include/observation.h
 * @Description: 
 * 
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved. 
 */
#pragma once

#include <su.h>

#include "cJSON.h"
#include "common.h"

typedef struct {
    int   shotID;       // Shot ID 
    float sx;           // Source x-coordinate
    float sz;           // Source z-coordinate
    int   isx;          // Source grid index x
    int   isz;          // Source grid index z
    float fpeak;        // Dominant frequency of source wavelet (Hz)

    float *src;         // Wavelet (nt)

    float *gx;          // Receiver x-coordinates (array or scalar)
    float *gz;          // Receiver z-coordinates
    int   *igx;         // Receiver grid indices x
    int   *igz;         // Receiver grid indices z

    int sampleRate;     // Recording sample rate
    int    ntr;         // Number of receiver traces
    int    ns;          // Number of time samples
    float dtr;          // Receiver interval
    float dt;           // Time sampling interval

    int    ix_min;      // Minimum x-index of recording window
    int    ix_max;      // Maximum x-index of recording window

    float **traces_vx;  // Recorded vx component (nt x ntr)
    float **traces_vz;  // Recorded vz component (nt x ntr)
    float **traces_p;   // Recorded p component (nt x ntr)

    float **precond;    // the pseudo-Hessian preconditioning of this shot
    float **grad;       // the gradient of this shot
} Record2D;


Record2D* alloc_record2d(int Nshot);

/* 
alloc_record2d: allocate an array of Nshot Record2D structures

Inputs:
 Nshot - number of shots

Outputs:
 return - pointer to allocated Record2D array

Notes:
 - All internal pointers are initialized to NULL
 - Caller must later allocate src, gx, gz, traces_vx, etc.
 - Caller must free the returned array via free_record2d() 
*/

void free_record2d(Record2D *record2d, int Nshot);

/* 
free_record2d： free Record2D array and all allocated subarrays

Inputs:
 records - pointer to Record2D array 
 Nshot   - number of shots in the array

Outputs:
 None
 
Notes:
 - Must be called when Record2D is no longer needed 
*/

int generate_shot_positions(
    const float *x, int nx,
    const float *z, int nz,
    float fsx, float ds, float sdel,
    int Nshot_in,
    float **sx, float **sz);

/* 
generate_shot_positions: generate shot coordinates and determine number of shots.

Inputs:
  x        - array of x-coordinates in model domain
  nx       - length of x array
  z        - array of z-coordinates in model domain
  nz       - length of z array
  fsx      - x-coordinate of the first shot
  ds       - interval between adjacent shots in x-direction
  sdel     - fixed z-coordinate for all shots
  Nshot_in - desired number of shots

Outputs:
  sx       - pointer to array pointer of x-positions of shots (malloc inside)
  sz       - pointer to array pointer of z-positions of shots (malloc inside)
  Nshot    - actual number of shots used (may be reduced)

Notes:
  - Caller must free sx and sz after use.
  - If fsx/sdel out of bounds, the function will call err and exit.
*/

Record2D* copy_record2d_header(const Record2D* record2d_src, int Nshot);

/* 
copy_record2d_header: copy only the header information (not wavefield data and traces_vx, traces_vz, and traces_p)

Inputs:
 record2d_src - pointer to source Record2D array
 Nshot        - number of shots in the array
 
Outputs:
 return       - pointer to newly allocated Record2D array with copied headers 
*/

void init_observe2d(Record2D *record2d, int Nshot, float *x, int nx, float *z, int nz,
    float dt, int nt, float *sx, float *sz, float *src,
    float offsmin, float offsmax, float dtr, float gdepth, int sampleRate);

/* 
init_observe2d: Initialize 2D seismic acquisition geometry and observation records

Inputs:
  record2d   - pointer to Record2D array (size Nshot)
  Nshot      - number of shots
  x, nx      - x-coordinates of model grid and size
  z, nz      - z-coordinates of model grid and size
  dt, nt     - time sampling interval and number of samples
  sx, sz     - shot coordinates arrays (size Nshot)
  src        - source wavelet array (size nt)
  offsmin    - minimum receiver offset
  offsmax    - maximum receiver offset
  dtr        - receiver spacing in x-direction
  gdepth       - receiver depth (z-coordinate)
  sampleRate - temporal downsampling factor

Outputs:
  record2d   - filled with shot and receiver information 
*/

void create_grid_observe2d(Record2D *rec, int Nshot, float dx, float dz, int nx);

/* 
create_grid_observe2d: convert real-world coordinates to grid indices

Input:
   rec   - array of Record2D
  Nshot - number of shots
  dx    - spatial grid spacing in x
  dz    - spatial grid spacing in z
  nx    - total number of x-grid points
  
Output:
  rec   - array of Record2D
*/

void init_observe2d_from_JSON(Record2D *record2d, cJSON *records, int Nshot, float *x, int nx, float *z, int nz,
    float dt, int nt, float *src, int sampleRate);
/*
init_observe2d_from_JSON: Initialize 2D seismic acquisition geometry and observation records from a JSON object

Inputs:
    record2d   - pointer to Record2D array (size Nshot)
    records    - cJSON object containing acquisition geometry information
    Nshot      - number of shots
    x, nx      - x-coordinates of model grid and size
    z, nz      - z-coordinates of model grid and size
    dt, nt     - time sampling interval and number of samples
    src        - source wavelet array (size nt)
    sampleRate - temporal downsampling factor

Outputs:
    record2d   - filled with shot and receiver geometry parsed from the JSON object
*/