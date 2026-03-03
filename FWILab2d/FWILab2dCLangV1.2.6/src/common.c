/*
 * @Author: QuYang Chen
 * @Date: 2025-09-17 19:05:53
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2026-01-22 20:16:44
 * @FilePath: /FWILab2dCLangV1.2.5/src/common.c
 * @Description:
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#include "common.h"

void transpose2d(float **array1, int n1, int n2, float **array2)
{
    /*
    transpose2d: transpose a 2D float array 

    Input:
        array1  - input 2D float array  (n1 x n2)
        n1      - number of array1 at second dimension 
        n2      - number of array1 at first dimension  

    Output:
        array2  - output transposed 2D float array (n2 x n1)
    */
    int i1, i2;

    for (i2 = 0; i2 < n2; i2++)
    {
        for (i1 = 0; i1 < n1; i1++)
        {
            array2[i1][i2] = array1[i2][i1];
        }
    }
}

float find_min_float(const float *array, int nx)
{
    /*
    find_min_float: find the minimum value in a 1D float array using OpenMP

    Input:
        arr   - input 1D float array
        n     - number of elements in the array 

    Output:
        return - minimum value (float) found in the array
    */
    float min_val = FLT_MAX;

#pragma omp parallel for reduction(min:min_val)
    for (int ix = 0; ix < nx; ix++)
    {
        if (array[ix] < min_val)
            min_val = array[ix];
    }

    return min_val;
}

void vadd(float *array1, float *array2, float *out, int nx)
{
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
    int ix;
#pragma omp parallel for private(ix) \
    firstprivate(nx)
    for (ix = 0; ix < nx; ix++)
    {
        /*out[ix]=cadd(array1[ix],array2[ix]);*/
        out[ix] = array1[ix] + array2[ix];
    }
}

void nvmul(float *in, float scale, float *out, int nx)
{
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

    int ix;
#pragma omp parallel for private(ix) \
    firstprivate(nx, scale)
    for (ix = 0; ix < nx; ix++)
    {
        /*out[ix]=crmul(in[ix],scale);*/
        out[ix] = scale * in[ix];
    }
}

void pad_model2d(float **mod, int nz, int nx,
              float **modp, int nlayer)
{
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
    int ix, iz, Nz, Nx;

    Nz = nz + 2 * nlayer;
    Nx = nx + 2 * nlayer;

    // Insert original model into center
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < nz; iz++)
        {
            modp[ix + nlayer][iz + nlayer] = mod[ix][iz];
        }
    }

    // Top padding (replicate first row of model)
    for (ix = 0; ix < Nx; ix++)
    {
        for (iz = 0; iz < nlayer; iz++)
        {
            modp[ix][iz] = modp[ix][nlayer];
        }
    }

    // Bottom padding (replicate last row of model)
    for (ix = 0; ix < Nx; ix++)
    {
        for (iz = nz + nlayer; iz < Nz; iz++)
        {
            modp[ix][iz] = modp[ix][nz + nlayer - 1];
        }
    }

    // Left padding (replicate first column)
    for (ix = 0; ix < nlayer; ix++)
    {
        for (iz = 0; iz < Nz; iz++)
        {
            modp[ix][iz] = modp[nlayer][iz];
        }
    }

    // Right padding (replicate last column)
    for (ix = nx + nlayer; ix < Nx; ix++)
    {
        for (iz = 0; iz < Nz; iz++)
        {
            modp[ix][iz] = modp[nx + nlayer - 1][iz];
        }
    }
}


void cap_model(const float *mod0, float *mod1, int n, float cap)
{
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

    int  i;
    float delta, delta_abs, delta_cap;

    if (cap <= 0.0f || cap >= 1.0f)
    {
        err("<error>: cap must be in range (0, 1)\n");
    }

#pragma omp parallel for private(i, delta, delta_abs, delta_cap) \
    firstprivate(n, cap)
    for (i = 0; i < n; i++)
    {
        delta = mod1[i] - mod0[i];
        delta_abs = fabsf(delta);
        delta_cap = cap * fabsf(mod0[i]);

        if (delta_abs > delta_cap)
        {
            // limit update direction and scale
            mod1[i] = mod0[i] + (delta > 0 ? 1.0f : -1.0f) * delta_cap;
        }
    }
}

void clip_model(float *mod, int n, float mod_min, float mod_max)
{
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
    
    int i;
#pragma ivdep
    for (i = 0; i < n; i++)
    {
        if (mod[i] > mod_max)
        {
            mod[i] = mod_max;
        }
        else if (mod[i] < mod_min)
        {
            mod[i] = mod_min;
        }
    }
}

void v2slow(float *v, float *slow, long len)
{ 
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
    #pragma ivdep
    for (long i = 0; i < len; i++)
    {
        if (v[i] != 0.0f)
            slow[i] = 1.0f / v[i];
        else
            slow[i] = 0.0f;
    }
}

void slow2v(float *slow, float *v, long len)
{ 
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
    #pragma ivdep
    for (long i = 0; i < len; i++)
    {
        if (slow[i] != 0.0f)
            v[i] = 1.0f / slow[i];
        else
            v[i] = 0.0f;
    }
}

void vp2density(const float *vp, float *rho, int n)
{
    /* 
    vp2density: convert P-wave velocity (vp) to density (rho) using Gardner relation

    Input:
        vp    - pointer to P-wave velocity array (float*), must be positive
        n     - number of elements in vp array (int)

    Output:
        rho   - pointer to allocated density array (float*), same size as vp

    Notes:
        - Uses empirical relationship: rho = 310 * vp^0.25 
    */
    
    int i;

    if (!vp || !rho)
    {
        err("<error>: function vp2density requires an input vp!");
    }

#pragma ivdep
    for (i = 0; i < n; i++)
    {
        if (vp[i] <= 0.0f)
        {
            err("<error>: all input P-wave velocities must be positive!");
        }
        rho[i] = 310.0f * powf(vp[i], 0.25f);
    }
}


void generate_frequency(float *f, float dt, int nt)
{
    /* 
    generate_frequency: fill angular frequency vector for FFT

    Input:
        dt - sampling interval in time domain (seconds, float)
        nt - number of time samples (int)

    Output:
        f  - preallocated float array of size nt (caller allocates) 
    */
    
    int it;
    float df = 2.0f * (float)PI / (nt * dt);
    int nf = nt / 2 + 1; // index of Nyquist frequency

    f[0] = 0.0f; // zero frequency

    for (it = 1; it < nf; it++)
    {
        f[it] = it * df;
        if (nt - it >= nf)
        {
            f[nt - it] = -f[it];
        }
    }
}


float calculate_peak_frequency(const float *s, int nt, float dt)
{
    /* 
    calculate_peak_frequency: estimate the peak (dominant) frequency of a time-domain signal.

    Input:
        s   - input signal array (float, length nt)
        nt  - number of samples
        dt  - sampling interval in seconds
        
    Output:
        return - peak frequency (Hz) 
    */
    
    float re, im, amp;
    int i;

    float df = 1.0f / (nt * dt); // frequency resolution
    float fpeak = 0.0f;          // peak frequency
    float max_amp = -1.0f;       // max amplitude

    fftwf_complex *out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (nt / 2 + 1)); // FFT output
    float *in = (float *)fftwf_malloc(sizeof(float) * nt);                                    // FFT input

    for (i = 0; i < nt; i++)
    {
        in[i] = s[i]; // copy input
    }

    fftwf_plan p = fftwf_plan_dft_r2c_1d(nt, in, out, FFTW_ESTIMATE); // create FFT plan
    fftwf_execute(p);                                                 // execute FFT

    for (i = 0; i <= nt / 2; i++)
    {
        re = out[i][0];
        im = out[i][1];
        amp = sqrtf(re * re + im * im); // magnitude spectrum

        if (amp > max_amp)
        {
            max_amp = amp;
            fpeak = i * df; // Hz
        }
    }

    fftwf_destroy_plan(p);
    fftwf_free(in);
    fftwf_free(out);

    return fpeak;
}


void calculate_gradient2d(float **f, int nx, int nz, float dx, float dz,
                float **dfdx, float **dfdz)
{
    
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

    int ix, iz;
    // X-direction derivative
#pragma omp parallel for private(ix, iz) \
    firstprivate(nx, nz, dx)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
            if (ix == 0)
            {
                dfdx[ix][iz] = (f[ix + 1][iz] - f[ix][iz]) / dx; // forward
            }
            else if (ix == nx - 1)
            {
                dfdx[ix][iz] = (f[ix][iz] - f[ix - 1][iz]) / dx; // backward
            }
            else
            {
                dfdx[ix][iz] = (f[ix + 1][iz] - f[ix - 1][iz]) / (2.0f * dx); // central
            }
        }
    }

    // Z-direction derivative
#pragma omp parallel for private(ix, iz) \
    firstprivate(nx, nz, dz)
    for (ix = 0; ix < nx; ix++)
    {
#pragma ivdep
        for (iz = 0; iz < nz; iz++)
        {
            if (iz == 0)
            {
                dfdz[ix][iz] = (f[ix][iz + 1] - f[ix][iz]) / dz;
            }
            else if (iz == nz - 1)
            {
                dfdz[ix][iz] = (f[ix][iz] - f[ix][iz - 1]) / dz;
            }
            else
            {
                dfdz[ix][iz] = (f[ix][iz + 1] - f[ix][iz - 1]) / (2.0f * dz);
            }
        }
    }
}


void interpolate_traces(float **tr, int ntr, int nt, int r, float **trs)
{
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
    int itr, i, k, it, ns, idx;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    
    if (ntr <= 0 || nt <= 0 || tr == NULL || trs == NULL)
    {
        err("<error>: invalid input (ntr=%d, nt=%d) in function interpolate_traces!", ntr, nt);
        return;
    }

    if (r <= 1)
    {
        // No interpolation: just copy (assume caller allocated ns == nt)
        memcpy(trs[0], tr[0], FSIZE * ntr * nt);
        return;
    }

    ns = (nt - 1) * r + 1; // output samples per trace

    // allocate complex buffer (single-precision)
    fftwf_complex *buf = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * ns);

    // create FFTW plans (in-place complex-to-complex)
    fftwf_plan p_forw = fftwf_plan_dft_1d(ns, buf, buf, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan p_back = fftwf_plan_dft_1d(ns, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);

    // build simple rectangular frequency window (length ns)
    float *win = ealloc1float(ns);
    memset(win, 0, ns * FSIZE);

    int n_keep = ns / (2 * r);
    for (k = 0; k < n_keep; ++k)
        win[k] = (float)r;
    for (k = ns - n_keep + 1; k < ns; ++k)
        win[k] = (float)r;

    // loop over traces (each trace is a column tr[itr][it])
    for (itr = 0; itr < ntr; ++itr)
    {
        // zero the complex buffer
        for (i = 0; i < ns; ++i)
        {
            buf[i][0] = 0.0f; // real
            buf[i][1] = 0.0f; // imag
        }

        // zero-stuffing: put original samples at indices it*r
        for (it = 0; it < nt; ++it)
        {
            idx = it * r;
            buf[idx][0] = tr[itr][it]; // real part
            // imag already zero
        }

        // forward FFT
        fftwf_execute(p_forw);

        // apply frequency window (multiply complex by scalar win[k])
        for (k = 0; k < ns; ++k)
        {
            buf[k][0] *= win[k];
            buf[k][1] *= win[k];
        }

        // inverse FFT
        fftwf_execute(p_back);

        // normalize and copy real part to output (FFTW's backward transform is unnormalized)
        for (it = 0; it < ns; ++it)
        {
            trs[itr][it] = buf[it][0] / (float)ns;
        }
    }

    // cleanup
    free1float(win);
    fftwf_destroy_plan(p_forw);
    fftwf_destroy_plan(p_back);
    fftwf_free(buf);
    fftwf_cleanup_threads();
}


void apply_traces_lowpass(float **tr, float **trs, int nt, int ntr,
                   float fpass, float fcut, float dt)
{
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
    
    int itr, it;
    int nfreq = nt; // full FFT length
    float df, taperLen, f;
    float *win, *in, *out;
    fftwf_complex *spec;
    fftwf_plan fplan, iplan;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    
    win = alloc1float(nfreq);

    // Allocate FFTW arrays
    spec = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (nt / 2 + 1));
    in = (float *)fftwf_malloc(sizeof(float) * nt);
    out = (float *)fftwf_malloc(sizeof(float) * nt);

    df = 1.0f / (nt * dt); // frequency resolution [Hz]
    if (fcut <= fpass)
    {
        err("<error>: fcut must be greater than fpass!");
    }

    // Build frequency window
    taperLen = fcut - fpass;
    #pragma omp parallel for private(it, f) \
    firstprivate(nfreq, df, fpass, fcut, taperLen)
    for (it = 0; it < nfreq; it++)
    {
        f = it * df;
        if (f < fpass)
        {
            win[it] = 1.0f;
        }
        else if (f >= fpass && f <= fcut)
        {
            win[it] = 0.5f * (1.0f + cosf((float)PI * (f - fpass) / taperLen));
        }
        else
        {
            win[it] = 0.0f;
        }
    }

    // Create FFT plan
    fplan = fftwf_plan_dft_r2c_1d(nt, in, spec, FFTW_ESTIMATE);
    iplan = fftwf_plan_dft_c2r_1d(nt, spec, out, FFTW_ESTIMATE);

    for (itr = 0; itr < ntr; itr++)
    {

        // Copy input trace
        for (it = 0; it < nt; it++)
        {
            in[it] = tr[itr][it];
        }

        // Forward FFT
        fftwf_execute(fplan);

        // Apply frequency window (only for positive freqs)
        for (it = 0; it < nt / 2 + 1; it++)
        {
            spec[it][0] *= win[it]; // real part
            spec[it][1] *= win[it]; // imag part
        }

        // Inverse FFT
        fftwf_execute(iplan);

        // Normalize by nt (FFTW does not normalize)
        for (it = 0; it < nt; it++)
        {
            trs[itr][it] = out[it] / nt;
        }
    }

    // Free
    fftwf_destroy_plan(fplan);
    fftwf_destroy_plan(iplan);
    fftwf_free(spec);
    fftwf_free(in);
    fftwf_free(out);
    free1float(win);
    fftwf_cleanup_threads();
}

void apply_traces_bandpass(float **tr, float **trs, int nt, int ntr, float fpass, float fcut, float dt)
{
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

    int itr, it;
    int nfreq = nt;  // full FFT length
    float df = 1.0f / (nt * dt); // frequency resolution
    float f, lowTaperLen, highTaperLen;
    float *win, *in, *out;
    fftwf_complex *spec;
    fftwf_plan fplan, iplan;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());

    win = alloc1float(nfreq);
    spec = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * (nt / 2 + 1));
    in   = (float *)fftwf_malloc(sizeof(float) * nt);
    out  = (float *)fftwf_malloc(sizeof(float) * nt);

    if (fcut <= fpass)
    {
        err("<error>: fcut must be greater than fpass!");
    }
    // ---- Build band-pass window ----
    lowTaperLen = fpass;
    highTaperLen = fcut - fpass;
    #pragma omp parallel for private(it, f) firstprivate(df, fpass, fcut, nfreq, lowTaperLen, highTaperLen)
    for (it = 0; it < nfreq; it++)
    {
        f = it * df;
        if (f < fpass)
        {
            // low-end cosine taper up
            win[it] = 0.5f * (1.0f + cosf(PI * (fpass - f) / lowTaperLen));
        }
        else if (f >= fpass && f <= fcut)
        {
            // high-end cosine taper down
            win[it] = 0.5f * (1.0f + cosf(PI * (f - fpass) / highTaperLen));
        }
        else
        {
            win[it] = 0.0f;
        }
    }

    // ---- FFT plans ----
    fplan = fftwf_plan_dft_r2c_1d(nt, in, spec, FFTW_ESTIMATE);
    iplan = fftwf_plan_dft_c2r_1d(nt, spec, out, FFTW_ESTIMATE);

    // ---- Process each trace ----
    for (itr = 0; itr < ntr; itr++)
    {
        // copy trace
        for (it = 0; it < nt; it++)
            in[it] = tr[itr][it];

        // forward FFT
        fftwf_execute(fplan);

        // apply bandpass window (positive freqs)
        for (it = 0; it < nt / 2 + 1; it++)
        {
            spec[it][0] *= win[it];
            spec[it][1] *= win[it];
        }

        // inverse FFT
        fftwf_execute(iplan);

        // normalize
        for (it = 0; it < nt; it++)
            trs[itr][it] = out[it] / nt;
    }

    // ---- Cleanup ----
    fftwf_destroy_plan(fplan);
    fftwf_destroy_plan(iplan);
    fftwf_free(spec);
    fftwf_free(in);
    fftwf_free(out);
    free1float(win);
    fftwf_cleanup_threads();
}

void taper_traces(float **tr, int nt, int ntr,
                 int taperLenTime, int taperLenTrace)
{
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
    
    int it, itr;
    float x;
    float *taperCoef;
    // --- Time-domain taper (top and bottom of each trace) ---
    if (taperLenTime > 0)
    {
        if (taperLenTime > nt / 3)
            taperLenTime = nt / 3;

        // Cosine taper coefficients
        taperCoef = alloc1float(taperLenTime + 1);
#pragma omp parallel for private(it, x) firstprivate(taperLenTime)
        for (it = 0; it <= taperLenTime; it++)
        {
            x = (float)PI * it / taperLenTime;
            taperCoef[it] = 0.5f * (1.0f + cosf(x));
        }

        // Top
#pragma omp parallel for private(itr, it) firstprivate(taperLenTime)
        for (itr = 0; itr < ntr; itr++)
        {
#pragma ivdep
            for (it = 0; it <= taperLenTime; it++)
            {
                tr[itr][it] *= taperCoef[taperLenTime - it];
            }
        }

        // bottom
#pragma omp parallel for private(itr, it) firstprivate(taperLenTime, ntr, nt)
        for (itr = 0; itr < ntr; itr++)
        {
#pragma ivdep
            for (it = 0; it <= taperLenTime; it++)
            {
                tr[itr][nt - taperLenTime - 1 + it] *= taperCoef[it];
            }
        }
        free1float(taperCoef);
    }

    // --- Trace-domain taper (left and right of each sample) ---
    if (taperLenTrace > 0)
    {
        if (taperLenTrace > ntr / 3)
            taperLenTrace = ntr / 3;

        taperCoef = alloc1float(taperLenTrace + 1);
#pragma omp parallel for private(itr, x) firstprivate(taperLenTrace)
        for (itr = 0; itr <= taperLenTrace; itr++)
        {
            x = (float)PI * itr / taperLenTrace;
            taperCoef[itr] = 0.5f * (1.0f + cosf(x));
        }

        // Left
#pragma omp parallel for private(itr, it) firstprivate(taperLenTrace, nt)
        for (itr = 0; itr <= taperLenTrace; itr++)
        {
#pragma ivdep
            for (it = 0; it < nt; it++)
            {
                tr[itr][it] *= taperCoef[taperLenTrace - itr];
            }
        }

        // Right
#pragma omp parallel for private(itr, it) firstprivate(taperLenTrace, ntr, nt)
        for (itr = 0; itr <= taperLenTrace; itr++)
        {
#pragma ivdep
            for (it = 0; it < nt; it++)
            {
                tr[ntr - taperLenTrace - 1 + itr][it] *= taperCoef[itr];
            }
        }

        free1float(taperCoef);
    }
}


void set_waterlayer_value(float **model, float value, int waterdepth, int nx, int nz)
{
    /* 
    set_waterlayer_value: overwrite the top part of the model with a constant water layer

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
    
    int ix, iz;
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < waterdepth; iz++)
        {
            model[ix][iz] = value;
        }
    }
}


void smooth2D(float **in, float **out, int nx, int nz, int N, int iszdirect)
{
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
    
    int ix, iz, icur, icur2;
    int weightStartIndex, weightEndIndex, dataStartIndex, dataEndIndex;
    double accumulator1, accumulator2;
    float *weight;

    weight = alloc1float(2 * N + 1);
#pragma ivdep
    for (icur = 0; icur < 2 * N + 1; icur++)
    {
        weight[icur] = cos(PI * (icur - N) / (2.0 * N));
        weight[icur] *= weight[icur];
    }

    if (iszdirect == 1)
    {
        for (iz = 0; iz < nz; iz++)
        {
            if (iz >= N)
            {
                weightStartIndex = 0;
                dataStartIndex = iz - N;
            }
            else
            {
                weightStartIndex = N - iz;
                dataStartIndex = 0;
            }
            if (nz - 1 - iz < N)
            {
                weightEndIndex = 2 * N - (N - (nz - 1 - iz));
                dataEndIndex = nz - 1;
            }
            else
            {
                weightEndIndex = 2 * N;
                dataEndIndex = iz + N;
            }
            for (ix = 0; ix < nx; ix++)
            {
                accumulator1 = 0.0;
                accumulator2 = 0.0;
#pragma ivdep
                for (icur = weightStartIndex, icur2 = dataStartIndex; icur <= weightEndIndex; icur++, icur2++)
                {
                    accumulator1 += weight[icur] * in[ix][icur2];
                    accumulator2 += weight[icur];
                }
                out[ix][iz] = accumulator1 / (accumulator2 + FLT_MIN);
            }
        }
    }
    else
    {
        for (ix = 0; ix < nx; ix++)
        {
            if (ix >= N)
            {
                weightStartIndex = 0;
                dataStartIndex = ix - N;
            }
            else
            {
                weightStartIndex = N - ix;
                dataStartIndex = 0;
            }
            if (nx - 1 - ix < N)
            {
                weightEndIndex = 2 * N - (N - (nx - 1 - ix));
                dataEndIndex = nx - 1;
            }
            else
            {
                weightEndIndex = 2 * N;
                dataEndIndex = ix + N;
            }
            for (iz = 0; iz < nz; iz++)
            {
                accumulator1 = 0.0;
                accumulator2 = 0.0;
#pragma ivdep
                for (icur = weightStartIndex, icur2 = dataStartIndex; icur <= weightEndIndex; icur++, icur2++)
                {
                    accumulator1 += weight[icur] * in[icur2][iz];
                    accumulator2 += weight[icur];
                }
                out[ix][iz] = accumulator1 / (accumulator2 + FLT_MIN);
            }
        }
    }
    free1float(weight);
}

void lowpass2d(
    float **s, int nz, int nx,
    float dz, float dx,
    float kpass, float kcut,
    float **s_filter)
{
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

    int iz, ix, taperLenX, taperLenZ;
    int padLen, Nz, Nx, N;
    float win, taperLen, KX, KZ, k;
    float **s_pad, *kx, *kz;

    fftwf_complex *fft_out;
    fftwf_plan fplan, iplan;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    
    // Validate cutoff settings
    if (kcut <= kpass)
    {
        err("<error>: kcut must be greater than kpass!");
    }

    // taper data
    taperLenX = roundf(4 / (kpass * dx));
    taperLenZ = roundf(4 / (kpass * dz));

    taper_traces(s, nz, nx, taperLenZ, taperLenX);

    // padding zeros to data
    padLen = 10 * MAX(taperLenX, taperLenZ);
    Nz = nz + 2 * padLen;
    Nx = nx + 2 * padLen;
    s_pad = alloc2float(Nz, Nx);
    pad_model2d(s, nz, nx, s_pad, padLen);

    // Allocate FFT arrays
    fft_out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Nx * (Nz / 2 + 1));

    // FFT plan
    fplan = fftwf_plan_dft_r2c_2d(Nx, Nz, s_pad[0], fft_out, FFTW_ESTIMATE);

    // Forward FFT
    fftwf_execute(fplan);

    // Frequency vectors
    kx = alloc1float(Nx);
    kz = alloc1float(Nz);
    generate_frequency(kx, dx, Nx);
    generate_frequency(kz, dz, Nz);

    // // Apply filter
    taperLen = kcut - kpass;

#pragma omp parallel for private(ix, iz, KX, KZ, k, win) \
    firstprivate(Nx, Nz, kpass, kcut, taperLen)
    for (ix = 0; ix < Nx; ix++)
    {
        for (iz = 0; iz <= Nz / 2; iz++)
        {
            KX = kx[ix] / (2.0f * (float)PI);
            KZ = kz[iz] / (2.0f * (float)PI);
            k = sqrtf(KX * KX + KZ * KZ);

            if (k < kpass)
            {
                win = 1.0f;
            }
            else if (k >= kpass && k <= kcut)
            {
                win = 0.5f * (1.0f + cosf((float)PI * (k - kpass) / taperLen));
            }
            else
            {
                win = 0.0f;
            }

            fft_out[ix * (Nz / 2 + 1) + iz][0] *= win;
            fft_out[ix * (Nz / 2 + 1) + iz][1] *= win;
        }
    }

    // Inverse FFT
    iplan = fftwf_plan_dft_c2r_2d(Nx, Nz, fft_out, s_pad[0], FFTW_ESTIMATE);
    fftwf_execute(iplan);

    // Normalize and crop to original size
    N = Nz * Nx;
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < nz; iz++)
        {
            s_filter[ix][iz] = s_pad[ix + padLen][iz + padLen] / N;
        }
    }

    /* Free memory */
    free2float(s_pad);
    fftwf_destroy_plan(fplan);
    fftwf_destroy_plan(iplan);
    fftwf_free(fft_out);
    free1float(kx);
    free1float(kz);
    fftwf_cleanup_threads();
}

void bandpass2d(
    float **s, int nz, int nx,
    float dz, float dx,
    float kpass, float kcut,
    float **s_filter)
{
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

    int iz, ix, taperLenX, taperLenZ;
    int padLen, Nz, Nx, N;
    float win, lowTaperLen, highTaperLen, KX, KZ, k;
    float **s_pad, *kx, *kz;

    fftwf_complex *fft_out;
    fftwf_plan fplan, iplan;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(omp_get_max_threads());
    
    // Validate cutoff settings
    if (kcut <= kpass)
    {
        err("<error>: kcut must be greater than kpass!");
    }

    // taper data
    taperLenX = roundf(4 / (kpass * dx));
    taperLenZ = roundf(4 / (kpass * dz));

    taper_traces(s, nz, nx, taperLenZ, taperLenX);

    // padding zeros to data
    padLen = 10 * MAX(taperLenX, taperLenZ);
    Nz = nz + 2 * padLen;
    Nx = nx + 2 * padLen;
    s_pad = alloc2float(Nz, Nx);
    pad_model2d(s, nz, nx, s_pad, padLen);

    // Allocate FFT arrays
    fft_out = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * Nx * (Nz / 2 + 1));

    // FFT plan
    fplan = fftwf_plan_dft_r2c_2d(Nx, Nz, s_pad[0], fft_out, FFTW_ESTIMATE);

    // Forward FFT
    fftwf_execute(fplan);

    // Frequency vectors
    kx = alloc1float(Nx);
    kz = alloc1float(Nz);
    generate_frequency(kx, dx, Nx);
    generate_frequency(kz, dz, Nz);

    // // Apply filter
    lowTaperLen = kpass;
    highTaperLen = kcut - kpass;

#pragma omp parallel for private(ix, iz, KX, KZ, k, win) \
    firstprivate(Nx, Nz, kpass, kcut, lowTaperLen, highTaperLen)
    for (ix = 0; ix < Nx; ix++)
    {
        for (iz = 0; iz <= Nz / 2; iz++)
        {
            KX = kx[ix] / (2.0f * (float)PI);
            KZ = kz[iz] / (2.0f * (float)PI);
            k = sqrtf(KX * KX + KZ * KZ);

            if (k < kpass)
            {
                win = 0.5f * (1.0f + cosf((float)PI * (kpass - k) / lowTaperLen));
            }
            else if (k >= kpass && k <= kcut)
            {
                win = 0.5f * (1.0f + cosf((float)PI * (k - kpass) / highTaperLen));
            }
            else
            {
                win = 0.0f;
            }

            fft_out[ix * (Nz / 2 + 1) + iz][0] *= win;
            fft_out[ix * (Nz / 2 + 1) + iz][1] *= win;
        }
    }

    // Inverse FFT
    iplan = fftwf_plan_dft_c2r_2d(Nx, Nz, fft_out, s_pad[0], FFTW_ESTIMATE);
    fftwf_execute(iplan);

    // Normalize and crop to original size
    N = Nz * Nx;
    for (ix = 0; ix < nx; ix++)
    {
        for (iz = 0; iz < nz; iz++)
        {
            s_filter[ix][iz] = s_pad[ix + padLen][iz + padLen] / N;
        }
    }

    /* Free memory */
    free2float(s_pad);
    fftwf_destroy_plan(fplan);
    fftwf_destroy_plan(iplan);
    fftwf_free(fft_out);
    free1float(kx);
    free1float(kz);
    fftwf_cleanup_threads();
}