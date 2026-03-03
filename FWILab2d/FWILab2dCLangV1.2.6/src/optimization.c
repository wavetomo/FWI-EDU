/*
 * @Author: QuYang Chen
 * @Date: 2025-09-25 17:22:42
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2026-01-15 16:29:10
 * @FilePath: /FWILab2dCLangV1.2.5/src/optimization.c
 * @Description:
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#include "optimization.h"

void apply_gradient_precondition(float *grad, const float *precond, int n, float alpha)
{
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
    
    int i;

    // Default stabilization factor
    if (alpha <= 0.0f)
    {
        alpha = 0.005f;
    }

    // Compute RMS of precond
    double pre_eng = 0.0;
#pragma omp parallel for reduction(+ : pre_eng)
    for (i = 0; i < n; i++)
    {
        pre_eng += (double)precond[i] * (double)precond[i];
    }
    float precond_rms = (float)sqrt(pre_eng / n);

// Apply preconditioning
#pragma omp parallel for private(i)
    for (i = 0; i < n; i++)
    {
        grad[i] = grad[i] / (precond[i] + alpha * precond_rms);
    }
}

void apply_gradient_shot_precondition(float *grad, float *precond, int n)
{
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
    
    int i;
    float avg, factor;

    // Compute average value of precond
    double sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
    for (i = 0; i < n; i++)
    {
        sum += (double)precond[i];
    }
    avg = (float)(sum / n);

// Apply preconditioning to gradient and update precond
#pragma omp parallel for private(i)
    for (i = 0; i < n; i++)
    {
        factor = 1.0f + precond[i] / (7.0f * avg);
        grad[i] = grad[i] / factor;
        precond[i] = precond[i] / factor;
    }
}

float calculate_objective_value(Record2D *record2d, int Nshot)
{
    /* 
    calculate_objective_value: Compute the objective function value for FWI (L2 norm)

    Input:
        record2d - array of shot gather structures (Record2D*), each containing:
                        residual_p - 2D residual matrix (observed - predicted)
        Nshot    - number of shots in record2d array (int)

    Output:
        fun      - scalar objective function value (least-squares misfit) (float)
    */
    
    int ishot, itrace, it;
    double fun;

    fun = 0.0;
#pragma omp parallel for reduction(+ : fun)
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        for (itrace = 0; itrace < record2d[ishot].ntr; itrace++)
        {
            for (it = 0; it < record2d[ishot].ns; it++)
            {
                fun += SQR(record2d[ishot].traces_p[itrace][it]);
            }
        }
    }

    fun *= 0.5; // Least-squares cost function
    return (float)fun;
}


float calculate_step_length(Record2D *record2d_pre0, Record2D *record2d_pre1, int Nshot)
{
    /* 
    calculate_step_length: optimal step length based on misfit reduction

    Input:
        record2d_pre0  - predicted data before update
        record2d_pre1  - predicted data after update
        Nshot          - number of shots

    Output:
        steplen        - estimated optimal step length 
    */
    
    int ishot, itr, it;
    double numerator = 0.0;
    double denominator = 0.0;
    float delta;

#pragma omp parallel for collapse(3) reduction(+ : numerator, denominator) private(ishot, itr, it, delta)
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        for (itr = 0; itr < record2d_pre0[ishot].ntr; itr++)
        {
            for (it = 0; it < record2d_pre0[ishot].ns; it++)
            {
                delta = record2d_pre0[ishot].traces_p[itr][it] - record2d_pre1[ishot].traces_p[itr][it];
                numerator += record2d_pre0[ishot].traces_p[itr][it] * delta;
                denominator += delta * delta;
            }
        }
    }
    if (denominator == 0.0)
    {
        warn("<warning>: denominator is zero in step length calculation, returning 0.\n");
        return 0.0f;
    }
    else
    {
        return (float)(numerator / denominator);
    }
}


float compute_model_perturbation(float fraction, float *slowness, float *gradient, float *out, int len)
{
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
    
    int icur, icur2, pos, peak;
    int histogram[1001];
    float Max, Min, interval, accumulator, scalar;
    float cap = 0.99; /*cap the perturbation to including 99% gradient possibility*/
    double average = 0;

    Max = 0.0;
    Min = 0.0;
#pragma ivdep
    for (icur = 0; icur < 1001; icur++) /*zero histogram*/
        histogram[icur] = 0;

    //--------------------
    // fprintf(stdout,"line:%d in file: %s\n",__LINE__,__FILE__);
    // fflush(stdout);
    //--------------------

#pragma ivdep
    for (icur = 0; icur < len; icur++)
    { /*find the maximum and minimum in the gradient*/
        if (gradient[icur] > Max)
            Max = gradient[icur];
        if (gradient[icur] < Min)
            Min = gradient[icur];
    }
    interval = (Max - Min) / 1000.0;

    //--------------------
    // fprintf(stdout,"Max=%6.5e,Min=%6.5e,interval=%6.5e,line:%d in file: %s\n",Max,Min,interval,__LINE__,__FILE__);
    // fflush(stdout);
    //--------------------

    for (icur = 0; icur < len; icur++)
    {
        pos = (int)((gradient[icur] - Min) / interval + 0.5);
        histogram[pos] = histogram[pos] + 1;
    }
    //--------------------
    // fprintf(stdout,"line:%d in file: %s\n",__LINE__,__FILE__);
    // fflush(stdout);
    //--------------------
    peak = 0;
    for (icur = 0; icur < 1001; icur++) /*zero histogram*/
    {
        if (histogram[icur] > peak)
        {
            peak = histogram[icur];
            pos = icur;
        }
    }
    icur = pos;
    icur2 = pos;

    //--------------------
    // fprintf(stdout,"line:%d in file: %s\n",__LINE__,__FILE__);
    // fflush(stdout);
    //--------------------
    /********************/
    // printf("pos=%d;\n",pos);

    accumulator = (float)(histogram[icur]) / (float)len;
    while (accumulator < cap)
    {
        if (icur < 1000)
        {
            icur++;
            accumulator += (float)(histogram[icur]) / (float)len;
        }
        if (icur2 > 0)
        {
            icur2--;
            accumulator += (float)(histogram[icur2]) / (float)len;
        }
    }
    /********************/
    // printf("icur=%d;icur2=%d;accumulator=%6.5e\n",icur,icur2,accumulator);

    //--------------------
    // fprintf(stdout,"line:%d in file: %s\n",__LINE__,__FILE__);
    // fflush(stdout);
    //--------------------

    if (ABS(Min + icur2 * interval) > ABS(Min + icur * interval))
        scalar = ABS(Min + icur2 * interval);
    else
        scalar = ABS(Min + icur * interval);
    for (icur = 0; icur < len; icur++)
        average += slowness[icur];
    scalar = (average / len * fraction) / scalar;
    for (icur = 0; icur < len; icur++)
    {
        out[icur] = -scalar * gradient[icur]; /*the minus sign means the perturbation is opposite to the gradient */
    }
    /********************
    printf("scalar=%6.5e;max=%6.5e;min=%6.5e;interval=%6.5e\n",scalar,Max,Min,interval);
    FILE* temp_fp=fopen("histogram.bin","w");
    fwrite(histogram,sizeof(int),1001,temp_fp);
    fclose(temp_fp);
    fflush(stdout);
    ***************************/
    return scalar;
}

void compute_adjoint_source2d(Record2D *record2d_pre, Record2D *record2d_obs,
                     float **vp, float freq, int filter, Record2D *record2d_adjsource)
{
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
   
    int nt, ntr;
    int itr, it;
    int ix, iz;
    int taperLenTime, taperLenTrace;
    double vp_sum;
    float vp_avg;
    float kpass;
    float **trace_transpose;

    nt = record2d_adjsource->ns;
    ntr = record2d_adjsource->ntr;

    // Step 1: Compute data residual
#pragma omp parallel for private(itr, it) \
    firstprivate(ntr, nt)
    for (itr = 0; itr < ntr; itr++)
    {
#pragma ivdep
        for (it = 0; it < nt; it++)
        {
            record2d_adjsource->traces_p[itr][it] =  record2d_pre->traces_p[itr][it] - record2d_obs->traces_p[itr][it];
        }
    }

    // Step 2: Estimate taper parameters based on wave propagation
    vp_sum = 0.0;
#pragma omp parallel for reduction(+ : vp_sum)
    for (itr = 0; itr < ntr; itr++)
    {
        ix = record2d_adjsource->igx[itr]; 
        iz = record2d_adjsource->igz[itr];
        vp_sum += vp[ix][iz];
    }
    vp_avg = (float)(vp_sum / ntr);

    if (freq == 0.0f)
    {
        freq = record2d_adjsource->fpeak;
    }
    kpass = freq / vp_avg;

    taperLenTime = (int)roundf((2.0f / freq) / record2d_adjsource->dt) + 1;
    taperLenTrace = (int)roundf((2.0f / kpass) / record2d_adjsource->dtr) + 1;

    /* Step 3: filtering */
    if (filter == 0) // lowpass
    {
        
        trace_transpose = alloc2float(ntr, nt);

        // Time-domain taper before filter
        taper_traces(record2d_adjsource->traces_p, nt, ntr, taperLenTime, 0);

        // Time-domain lowpass filtering
        apply_traces_lowpass(record2d_adjsource->traces_p, record2d_adjsource->traces_p, nt, ntr,
                      freq, 2.0f * freq, record2d_adjsource->dt);

        // Space-domain taper
        taper_traces(record2d_adjsource->traces_p, nt, ntr, 0, taperLenTrace);

        // Space-domain lowpass filtering (transpose in/out)
        transpose2d(record2d_adjsource->traces_p, nt, ntr, trace_transpose);

        apply_traces_lowpass(trace_transpose, trace_transpose, ntr, nt, kpass, 2.0f * kpass, record2d_adjsource->dtr);
        
        transpose2d(trace_transpose, ntr, nt, record2d_adjsource->traces_p);

        free2float(trace_transpose);
    }
    else if(filter == 1) // bandpass
    {
        
        trace_transpose = alloc2float(ntr, nt);

        // Time-domain taper before filter
        taper_traces(record2d_adjsource->traces_p, nt, ntr, taperLenTime, 0);

        // Time-domain bandpass filtering
        apply_traces_bandpass(record2d_adjsource->traces_p, record2d_adjsource->traces_p, nt, ntr,
                      freq, 2.0f * freq, record2d_adjsource->dt);

        // Space-domain taper
        taper_traces(record2d_adjsource->traces_p, nt, ntr, 0, taperLenTrace);

        // Space-domain bandpass filtering (transpose in/out)
        transpose2d(record2d_adjsource->traces_p, nt, ntr, trace_transpose);

        apply_traces_bandpass(trace_transpose, trace_transpose, ntr, nt, kpass, 2.0f * kpass, record2d_adjsource->dtr);

        transpose2d(trace_transpose, ntr, nt, record2d_adjsource->traces_p);

        free2float(trace_transpose);
    }

    // Step 4: Final taper in both time and space domains
    taper_traces(record2d_adjsource->traces_p, nt, ntr, taperLenTime, taperLenTrace);
}