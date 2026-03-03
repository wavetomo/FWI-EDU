/*
 * @Author: QuYang Chen
 * @Date: 2025-09-17 03:53:16
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-11-01 16:40:50
 * @FilePath: /FWILab2dCLangV1.2.3/include/common.h
 * @Description: 
 * 
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved. 
 */
#pragma once

#include <su.h>

#include <fftw3.h>
#include <omp.h>

void transpose2d(float **array1, int n1, int n2, float **array2);

/*
transpose2d: transpose a 2D float array 

Input:
   array1  - input 2D float array  (n1 x n2)
   n1      - number of array1 at second dimension 
   n2      - number of array1 at first dimension  
   
Output:
   array2  - output transposed 2D float array (n2 x n1)
*/

float find_min_float(const float *array, int nx);

/*
find_min_float: find the minimum value in a 1D float array using OpenMP

Input:
   arr   - input 1D float array
   n     - number of elements in the array 

Output:
   return - minimum value (float) found in the array
*/

void vadd(float *array1, float *array2, float *out, int nx);

/*
vadd: perform element-wise addition of two 1D float arrays

Input:
   array1  - first input array (length = nx)
   array2  - second input array (length = nx)
   nx      - number of elements in each array

Output:
   out     - output array (length = nx)
               stores the result of array1 + array2
*/

void nvmul(float *in, float scale, float *out, int nx);

/*
nvmul: multiply each element of a 1D array by a scalar

Input:
   in      - input 1D array (length = nx)
   scale   - scalar multiplier
   nx      - number of elements in the array

Output:
   out     - output 1D array (length = nx)
               each element is computed as out[ix] = scale * in[ix]
*/
   
void v2slow(float *v, float *slow, long len);

/*
v2slow: convert velocity into slowness

Input:
   v: array storing velocity
   slow: array storing slowness of velocity
   len: number of elements of arrays

Output:
   slow: array storing slowness of velocity

Note: if v=0, slow=0
*/


void slow2v(float *slow, float *v, long len);

/*
slow2v：convert slowness into velocity

Input:
   slow: array storing slowness of velocity
   v: array storing velocity
   len: number of elements of arrays

Output:
   v: array storing velocity

Note: if slow=0, v=0
*/

void pad_model2d(float **mod, int nz, int nx,
    float **modp, int nlayer);
    
/* 
pad_model2d: pad a 2D model with boundary layers

Input:
   mod     - original 2D model (nz x nx)
   nz, nx  - original model size
   nlayer  - number of padding layers

 Output:
   modp    - padded model ((nz+2*nlayer) x (nx+2*nlayer))
             memory must be allocated by caller

Note: modp must be allocated with size [Nz][Nx] beforehand.
*/

void cap_model(const float *mod0, float *mod1, int n, float cap);

/* 
cap_model: limit the relative update of model parameters

Input:
   mod0 - original model 
   mod1 - updated model 
   n    - length of array
   cap  - scalar, 0 < cap < 1  

Output:
   mod1 - updated model
*/

void clip_model(float *mod, int n, float mod_min, float mod_max);

/* 
clip_model: clip 2D model values within [mod_min, mod_max]

Input:
   mod     - input model parameters 
   n       - length of array
   mod_min - minimum allowed value
   mod_max - maximum allowed value

Output:
   mod     - output model parameters
*/

void v2slow(float *v, float *slow, long len);
/*
v2slow: convert velocity into slowness

Input:
   v: array storing velocity
   slow: array storing slowness of velocity
   len: number of elements of arrays

Output:
   slow: array storing slowness of velocity

Note: if v=0, slow=0
*/

void slow2v(float *slow, float *v, long len);
/*
  slow2v：convert slowness into velocity

  Input:
     slow: array storing slowness of velocity
     v: array storing velocity
     len: number of elements of arrays

  Output:
     v: array storing velocity

 Note: if slow=0, v=0
*/

void vp2density(const float *vp, float *rho, int n);

/* 
vp2density: convert P-wave velocity (vp) to density (rho) using Gardner relation

Input:
   vp    - pointer to P-wave velocity array (float*), must be positive
   n     - number of elements in vp array (int)

Output:
   rho   - pointer to allocated density array (float*), same size as vp

Note: Uses empirical relationship: rho = 310 * vp^0.25
*/

void calculate_gradient2d(float **f, int nx, int nz, float dx, float dz,
    float **dfdx, float **dfdz);

/*
calculate_gradient2d: compute gradients of 2D array (like Matlab gradient)
 
Input:
   f[nx][nz]   - input 2D array
   nx, nz      - dimensions
   dx, dz      - grid spacing

Output:
   dfdx[nx][nz] - gradient along x (rows)
   dfdz[nx][nz] - gradient along z (cols)
*/

void generate_frequency(float *f, float dt, int nt);

/*
generate_frequency: fill angular frequency vector for FFT

Input:
   dt - sampling interval in time domain (seconds, float)
   nt - number of time samples (int)

Output:
   f  - preallocated float array of size nt (caller allocates)
*/

float calculate_peak_frequency(const float *s, int nt, float dt);

/*
calculate_peak_frequency: estimate the peak (dominant) frequency of a time-domain signal.
Input:
   s   - input signal array (float, length nt)
   nt  - number of samples
   dt  - sampling interval in seconds
   
Output:
   return - peak frequency (Hz)
*/

void interpolate_traces(float **tr, int ntr, int nt, int r, float **trs);

/*
interpolate_traces: Interpolate traces in time direction by factor r using zero-stuffing
 and frequency-domain low-pass filtering (single-precision FFTW).

Input:
   tr   - input traces, column-major: tr[ntr][nt]
          (each tr[itr] is a pointer to nt float samples)
   ntr  - number of traces (int)
   nt   - number of time samples in input traces (int)
   r    - integer interpolation factor (>1)
   trs  - output traces, column-major: trs[ntr][ns]
          (memory must be allocated by caller, each trs[itr] points to ns floats)

Output:
   trs  - filled with interpolated traces, ns = (nt - 1)*r + 1
*/

void apply_traces_lowpass(float **tr, float **trs, int nt, int ntr,
    float fpass, float fcut, float dt);

/* 
apply_traces_lowpass: Applies a low-pass cosine-taper filter to 2D trace data

Input:
   tr     - input traces [nt x ntr], column-major (trace = tr[:, itr])
   nt     - number of time samples per trace (int)
   ntr    - number of traces (int)
   fpass  - passband frequency (Hz) (float)
   fcut   - stopband frequency (Hz), must be > fpass (float)
   dt     - sampling interval (s) (float)

Output:
   trs    - filtered traces [nt x ntr], same size as input (float **)
*/

void apply_traces_bandpass(float **tr, float **trs, int nt, int ntr, 
   float fpass, float fcut, float dt);

/*
apply_traces_bandpass: Applies a band-pass cosine-taper filter to 2D trace data

Input:
   tr     - input traces [nt x ntr], column-major (trace = tr[:, itr])
   nt     - number of time samples per trace (int)
   ntr    - number of traces (int)
   fpass  - passband frequency (Hz) (float)
   fcut   - stopband frequency (Hz), must be > fpass (float)
   dt     - sampling interval (s) (float)

Output:
   trs    - filtered traces [nt x ntr], same size as input (float **)

Note: 
      The transition zone is smoothed using a cosine-squared window.
      1.0 |   fpass
            |  / \
            | /   \
            |/     \___________
            0        fcut      (HZ)
*/

void taper_traces(float **tr, int nt, int ntr,
    int taperLenTime, int taperLenTrace);

/*
taper_traces: Apply cosine taper to the edges of seismic traces

Input/Output:
   tr              - seismic traces [nt x ntr], column-major (float **tr, tr[time][trace])

Input:
   nt              - number of time samples (int)
   ntr             - number of traces (int)
   taperLenTime    - taper length in time direction (samples) (int)
   taperLenTrace   - taper length in trace direction (traces) (int)

Output:
   tr              - tapered seismic traces 
*/

void set_waterlayer_value(float **model, float value, int waterdepth, int nx, int nz);

/*
set_waterlayer_value: Overwrite the top part of the model with a constant water layer

Input/Output:
   model          - 2D velocity/density model [nx x nz], column-major (float **model, model[x][z])

Input:
   value          - constant value to assign to the water layer (float)
   waterdepth     - number of grid cells in depth direction to overwrite (int)
   nx             - number of grid points in x-direction (int)
   nz             - number of grid points in z-direction (int)

Output:
   model          - updated model with top 'waterdepth' layers set to 'value'
*/

void smooth2D(float **in,float **out,int nx,int nz,int N,int iszdirect);

/*
smooth2D: apply cosine-squared kernel smoothing along one dimension of a 2D array

Input:
   in              - input 2D data [nz x nx], column-major (float **in, in[z][x])
   nz              - number of samples along the 1st dimension
   nx              - number of samples along the 2nd dimension
   N               - half window length of smoothing
   iszdirect       - smoothing direction flag
                      =1    smooth along the 1st dimension (z-direction)
                      else  smooth along the 2nd dimension (x-direction)

Input/Output:
   out             - output 2D data [nz x nx], must be allocated outside

Note:
   The cosine-squared smoothing acts like a low-pass filter in the frequency domain.
   The filter has amplitude 1 at 0 Hz, decreases to 0.5 at 1/(2*N*dt),
   and reaches its first notch at 1/(N*dt).
*/

void lowpass2d(
    float **s, int nz, int nx,
    float dz, float dx,
    float kpass, float kcut,
    float **s_filter
);

/*
lowpass2d: apply 2D low-pass filter in the wavenumber domain

Input:
   s               - input 2D scalar field [nz x nx], column-major (float **s, s[z][x])
   nz              - number of samples along the 1st dimension (z-direction)
   nx              - number of samples along the 2nd dimension (x-direction)
   dz              - grid spacing along the 1st dimension (z-direction)
   dx              - grid spacing along the 2nd dimension (x-direction)
   kpass           - passband wavenumber [1/m], components below are fully preserved
   kcut            - cutoff wavenumber [1/m], components above are fully suppressed

Output:
   s_filter        - output 2D scalar field [nz x nx], must be allocated outside
*/

void bandpass2d(
   float **s, int nz, int nx,
   float dz, float dx,
   float kpass, float kcut,
   float **s_filter);

/* 
bandpass2d: apply 2D band-pass filter in the wavenumber domain

Input:
      s               - input 2D scalar field [nz x nx], column-major (float **s, s[z][x])
      nz              - number of samples along the 1st dimension (z-direction)
      nx              - number of samples along the 2nd dimension (x-direction)
      dz              - grid spacing along the 1st dimension (z-direction)
      dx              - grid spacing along the 2nd dimension (x-direction)
      kpass           - passband wavenumber [1/m], components below are fully preserved
      kcut            - cutoff wavenumber [1/m], components above are fully suppressed

Output:
      s_filter        - output 2D scalar field [nz x nx], must be allocated outside
*/