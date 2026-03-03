/*
 * @Author: QuYang Chen
 * @Date: 2025-09-17 04:20:48
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-10-31 15:44:53
 * @FilePath: /FWILab2dCLangV1.2.3/include/optimization.h
 * @Description: 
 * 
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved. 
 */
#pragma once

#include <su.h>

#include "observation.h"

void compute_adjoint_source2d(Record2D *record2d_pre, Record2D *record2d_obs,
    float **vp, float freq, int filter, Record2D *record2d_adjsource);

/*
compute_adjoint_source2d: Compute adjoint source from observed and predicted data

Inputs:
  record2d_pre   - structure with predicted data and acquisition parameters:
                  .traces_p : predicted pressure traces [nt x ntr]
                  .igx      : receiver x-indices (1-based) [1 x ntr]
                  .igz      : receiver z-indices (1-based) [1 x ntr]
                  .dt       : time sampling interval
                  .dtr      : spatial sampling interval along receiver line
  record2d_obs   - structure with observed pressure traces:
                  .traces_p : observed pressure traces [nt x ntr]
  vp             - 2D velocity model [nz x nx]
  freq           - central frequency (Hz) of the bandpass filter
  filter         - flag (0 , 1 or 2),  0: lowpass; 1: bandpass; 2: allpass

Output:
  record2d_adjsource.traces_p - adjoint source after differencing and optional filtering [nt x ntr]
*/

void apply_gradient_precondition(float *grad, const float *precond, int n, float alpha);

/* 
apply_gradient_precondition: Apply preconditioning to the gradient vector

Input:
  grad     - pointer to gradient array (float*)
  precond  - pointer to preconditioning array (float*), same size as grad
  n        - number of elements in grad and precond (int)
  alpha    - stabilization factor (float, optional, default: 0.005)

Output:
  grad     - gradient array is updated in-place with preconditioned values

Note:
  Each element of grad is divided by (precond[i] + alpha * RMS(precond)),
  where RMS(precond) is the root-mean-square of the precond array.
  This improves numerical stability during optimization.
*/

void apply_gradient_shot_precondition(float *grad, float *precond, int n);

/* 
apply_gradient_shot_precondition: Apply preconditioning to a single shot gradient

Input/Output:
  grad     - pointer to gradient array (float*), size = n
             updated in-place with preconditioned values
  precond  - pointer to preconditioning array (float*), size = n
             updated in-place for next iteration
  n        - number of elements in grad and precond (int)

Note:
  Applies stabilization-based preconditioning to the input gradient
  and updates the preconditioning array. The gradient is normalized
  by a factor related to the average of the preconditioning values,
  improving convergence behavior.
*/

float calculate_objective_value(Record2D *record2d, int Nshot);

/* 
calculate_objective_value: Compute the objective function value for FWI (L2 norm)

Input:
  record2d - array of shot gather structures (Record2D*), each containing:
                 residual_p - 2D residual matrix (observed - predicted)
  Nshot    - number of shots in record2d array (int)

Output:
  fun      - scalar objective function value (least-squares misfit) (float)
 */

float compute_model_perturbation(float fraction,float *slowness, float *gradient,float *out,int len);

/* 
compute_model_perturbation: Calculate slowness model perturbation based on the gradient

Input:
  fraction  - fraction of the average slowness to scale the perturbation (float, e.g., 0.01)
  slowness  - pointer to current slowness model array (float*)
  gradient  - pointer to gradient array of the current model (float*)
  len       - length of the model arrays (int)

Output:
  out       - pointer to output perturbation model (float*)
              [NB] can be the same as slowness or gradient

Return:
  scaling factor applied to the gradient (float)
*/

float calculate_step_length(Record2D *record2d_pre0, Record2D *record2d_pre1, int Nshot);

/* 
calculate_step_length: optimal step length based on misfit reduction

Input: 
  record2d_pre0  - predicted data before update
  record2d_pre1  - predicted data after update
  Nshot          - number of shots
  
Output: 
  steplen        - estimated optimal step length 
*/