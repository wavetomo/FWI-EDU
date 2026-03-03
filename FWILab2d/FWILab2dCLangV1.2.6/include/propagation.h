/*
 * @Author: QuYang Chen
 * @Date: 2025-09-20 23:45:33
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2026-01-15 15:18:23
 * @FilePath: /FWILab2dCLangV1.2.5/include/propagation.h
 * @Description: 
 * 
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved. 
 */
#pragma once

#include <su.h>

#include "common.h"
#include "observation.h"
#include "optimization.h"
#include "dataio.h"
#include "checkpoint.h"

#define MST 5 // wavefield and gradient saving interval in time 

void acoustic_iso2d_fd_it_6order(float **restrict p2, float **restrict p1, float **restrict p0,
  float **restrict invrho,  float **restrict invrhox, float **restrict invrhoz, float **restrict K, 
  float **restrict coef1, float **restrict coef2, float **restrict coef3,
  float invdx, float invdz, float invdxdx, float invdzdz, 
  int ixStart, int ixEnd, int izStart, int izEnd
);

/*  
acoustic_iso2d_fd_it_6order: 2d isotropic second-order acoustic wave equation forward modeling 
by  6th order precision finite difference scheme

Input:
  p0        - acoustic pressure wavefiled at previous moment
  p1        - acoustic pressure wavefiled at current moment
  p2        - acoustic pressure wavefiled at next moment
  invrhox   - first order derivative of inverse of density along x direction
  invrhoz   - first order derivative of inverse of density along z direction
  invrho    - inverse of density
  K         - K = rho * vp^2
  coef1     - coefficent of finite difference scheme
  coef2     - coefficent of finite difference scheme
  coef3     - coefficent of finite difference scheme
  invdx     - inverse of grid size along x direction
  invdz     - inverse of grid size along z direction
  invdxdx   - square of inverse of grid size along x direction
  invdzdz   - square of inverse of grid size along z direction
  ixStart   - starting grid index in x-direction for computation
  ixEnd:    - ending grid index in x-direction for computation
  izStart   - starting grid index in z-direction for computation
  izEnd     - ending grid index in z-direction for computation

Output: 
  p2        - acoustic pressure wavefiled at next moment
*/

void acoustic_iso2d_freesurface(float **p2, int Nx, int Nz, int pmlThick);

/*
acoustic_iso2d_freesurface: apply free-surface boundary condition to 2D acoustic wavefield

  Input:
  p2       - 2D acoustic pressure wavefield at next moment (float**)
  Nx       - number of grid samples along x (int)
  Nz       - number of grid samples along z (int)
  pmlThick - thickness of PML boundary in grid points (int)

Output:
  p2       - updated wavefield with free-surface condition applied (modified in-place)
*/
void mealCoefs2d(float **D, float **Beta, float **vp, float **vs, float fpeak, float eal_alpha,
  float dz, float dx, int NZ, int NX, int nlayer);

/* 
meal2d_bc: Construct 2D MEAL (Modified Efficient Absorbing Layer) boundary condition

Input:
  D         - damping function values (float**, size [NZ][NX])
  Beta      - scaling factor (float**, size [NZ][NX])
  vp        - P-wave velocity model (float**, size [NZ][NX])
  vs        - S-wave velocity model (float**, size [NZ][NX])
  fpeak     - dominant frequency of the source (float)
  eal_alpha - order of Gaussian weight function (float)
  dz        - grid spacing in z-direction (1st dimension, float)
  dx        - grid spacing in x-direction (2nd dimension, float)
  NZ        - number of grid samples along z (int)
  NX        - number of grid samples along x (int)
  nlayer    - number of absorbing boundary layers (int)

Output:
  D, Beta   - updated boundary condition arrays (modified in-place)
*/

void iso_acoustic2d_propagation_engine(Record2D *record2d,
    Acoustic2dCheckPoint *chkpt, /* allocated if mode==1, else NULL */
    float **vpm, float **rhom, int nx, int nz,
    float dx, float dz, float dt, int nt,
    int pmlThick, int freesurface, int precondition,
    int mode, int iter, int verbose,
    float **grad /* allocated if mode==1, else NULL */);

/*
iso_acoustic2d_propagation_engine: 2D isotropic acoustic wave propagation engine

Input:
   record2d     - array of shot gather structures containing sources, receivers, and traces (Record2D*)
   chkpt        - pointer to a constructs which contains poiters of
                        wavefiled data (allocated if mode==1, else NULL) (Acoustic2dCheckPoint*)
   vpm          - 2D P-wave velocity model (nz x nx) (float**)
   rhom         - 2D density model (nz x nx) (float**)
   nx, nz       - number of grid points in x and z directions (int)
   dx, dz       - grid spacing in x and z directions [m] (float)
   dt           - time step [s] (float)
   nt           - number of time steps (int)
   pmlThick     - thickness of PML boundary in grid points (int)
   freesurface  - free surface flag (1 = yes, 0 = no) (int)
   precondition - apply pseudo-Hessian preconditioning (1 = yes, 0 = no) (int)
   mode         - operation mode (0: forward modeling, 1: forward+wavefield+precondition, 2: adjoint/gradient) (int)
   iter         - current iteration index (int)
   verbose      - verbosity level 
                      0 = no display
                      1 = display main information without wavefields
                      2 = display all information with wavefields
   grad         - 2D gradient array (allocated if mode==1, else NULL) (float**)

Output:
   chkpt        - updated checkpoint wavefield array if mode==1
   record2d     - updated shot gathers with synthetic traces
   grad         - computed gradient if mode==1 or 2
*/

void iso_acoustic2d_propagation(
    float **vp, float **rho, Record2D *record2d_obs,
    int Nshot, int nz, int nx,
    float dx, float dz, float dt, int nt,
    int pmlThick, int freesurface,
    float freq, float filter,
    int precondition, int problem, int iter,
    int verbose, 
    Record2D *record2d_pre, float **grad    
);

/* 
iso_acoustic2d_propagation: 2D isotropic acoustic wave propagation and gradient computation

Input:
    vp           - 2D P-wave velocity model (nz x nx) (float**)
    rho          - 2D density model (nz x nx) (float**)
    record2d_obs - array of observed 2D shot records (Record2D*)
    Nshot        - number of shots (int)
    nz, nx       - number of grid points in z and x directions (int)
    dx, dz       - grid spacing in x and z directions [m] (float)
    dt           - time step [s] (float)
    nt           - number of time steps (int)
    pmlThick     - thickness of PML absorbing boundary [grid points] (int)
    freesurface  - free surface flag (1 = yes, 0 = no) (int)
    freq         - reference frequency for residual filtering [Hz] (float)
    filter       - filter parameter for residual computation (float)
    precondition - apply pseudo-Hessian preconditioning (1 = yes, 0 = no) (int)
    problem      - operation mode flag:
                      0 = forward modeling for computing record only
                      1 = full waveform inversion (FWI) gradient computation
                      2 = forward + adjoint residual computation 
    iter         - current iteration index (int)
    verbose      - verbosity level 
                      0 = no display
                      1 = display main information without wavefields
                      2 = display all information with wavefields

Output:
    record2d_pre - predicted 2D shot records (Record2D*) if problem==0,1 or residuals if problem==2
    grad         - computed model gradient (nz x nx) (float**)

Description:
    This routine performs 2D isotropic acoustic wave propagation for forward modeling,
    residual computation, and gradient formation. Depending on the problem type, it:
      - simulates synthetic shot gathers,
      - computes adjoint sources (residuals),
      - and accumulates gradients for FWI.

Notes:
    - The environment variable "shotInterval" controls output interval for intermediate results.
    - Gradient preconditioning and low-pass filtering are applied when enabled.
    - Requires subroutines:
          iso_acoustic2d_propagation_engine()
          compute_adjoint_source2d()
          apply_gradient_precondition()
          lowpass2d()
*/