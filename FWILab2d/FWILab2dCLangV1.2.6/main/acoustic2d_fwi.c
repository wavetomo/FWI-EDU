/*
 * @Author: QuYang Chen
 * @Date: 2024-11-25 16:01:32
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2026-01-23 15:23:38
 * @FilePath: /FWILab2dCLangV1.2.5/main/acoustic2d_fwi.c
 * @Description:
 *
 * Copyright (c) 2024 by WaveTomo, All Rights Reserved.
 */

#include <par.h>
#include <su.h>
#include <segy.h>

#include <omp.h>

#include "dataio.h"
#include "common.h"
#include "observation.h"
#include "propagation.h"
#include "optimization.h"

char *sdoc[] = {
    "acoustic2d_fwi - 2D acoustic full waveform inversion with frequency multi-scale strategy",
    "using checkpoint method for wavefield reconstruction",
    "",
    "Usage:",
    "    acoustic2d_fwi [parameters]",
    "",
    "Required Parameters:",
    "    nx=             Number of grid points in the horizontal (x) direction",
    "    nz=             Number of grid points in the vertical (z) direction",
    "    nt=             Number of time steps",
    "    dx=             Grid spacing in x direction (unit: m)",
    "    dz=             Grid spacing in z direction (unit: m)",
    "    dt=             Time step (unit: s)",
    "    Nshot=          Number of shots [required]",
    "    srcfile=        Binary file containing source wavelet (type: float)",
    "    recordfile=     Binary file containing observed shot gathers (type: float)",
    "    vpfile=         Binary file containing P-wave velocity model (type: float)",
    "",
    "Optional Parameters:",
    "    rhofile=        Binary file containing density model (type: float)",
    "                    [default: disabled if Gardner=1]",
    "    Gardner=        Apply Gardner's density relation (1=yes, 0=no)",
    "                    [default: 1]",
    "",
    "    freesurface=    Apply free-surface boundary condition (1=yes, 0=no)",
    "                    [default: 0]",
    "    pmlThick=       Thickness of the MEAL boundary (unit: cell)",
    "                    [default: 30]",
    "",
    "    vpmax=          Maximum allowed velocity value [default: 10000.0]",
    "    vpmin=          Minimum allowed velocity value [default: 0.001]",
    "    rhomax=         Maximum allowed density value [default: 10000.0]",
    "    rhomin=         Minimum allowed density value [default: 0.001]",
    "",
    "Water Layer Parameters:",
    "    waterdepth=     Depth of water layer (unit: m) [default: 0]",
    "    vp_water=       P-wave velocity in water (unit: m/s) [default: 1500.0]",
    "    rho_water=      Density in water (unit: kg/m^3) [default: 1000.0]",
    "",
    "Frequency Multi-scale Parameters:",
    "    iteration=      Number of iterations for each frequency band (type: int array)",
    "                    e.g., iteration=10,8,6",
    "    freq=           Frequency (Hz) for each frequency band, same length as iteration array",
    "                    e.g., freq=5,10,15",
    "    filter=         Apply filter for each frequency band (0=lowpass, 1=bandpass, 2=allpass), same length as iteration array",
    "                    e.g., filter=0,0,2",
    "    rx=             Smoothing radius in x direction for each frequency band (unit: cell)",
    "    rz=             Smoothing radius in z direction for each frequency band (unit: cell)",
    "",
    "Note: The lengths of iteration, freq, filter, rx, and rz must be identical.",
    "",
    "Other Parameters:",
    "    precondition=   Apply gradient spatial preconditioning (1=yes, 0=no) [default: 1]",
    "    verbose=        verbosity level",
    "                    0 = dispaly minimal information",
    "                    1 = display more information without wavefields, good for QC",
    "                    2 = display all information with wavefields, good for code development",
    "",
    "Notes:",
    "    - Gardner=1 applies Gardner's equation to compute density automatically.",
    "    - All grid-based binary files must be stored in C float (32-bit) column-major format (Seismic Unix style).",
    "    - Frequency-domain staging helps stabilize inversion by starting from low to high frequencies.",
    "    - Using the MEAL absorbing boundary, following the reference:",
    "      ZHANG PingMin and YAO Gang. 2024. Numerical simulation of seismic wavefield based on modified EAL boundary condition. ",
    "      Chinese Journal of Geophysics (in Chinese), 67(1): 261-276, doi: 10.6038/cjg2023Q0746.",
    "",
    NULL};

int main(int argc, char *argv[])
{
    int verbose, freesurface;
    int nx, nz, nt, pmlThick;
    float dx, dz, dt;

    int ix, iz, ishot, i, k, iter;

    // source
    int Nshot;
    float *src;

    // model
    int Gardner, waterdepth;
    float vp_water, rho_water;
    float **vp, **rho, **slow, **grad, **perturb;

    // inversion parameters
    float vpmin, vpmax, rhomin, rhomax, scalar, steplen;
    float *freq, *fun1, *fun2, *timeCost;
    int *iteration, *filter, *rx, *rz;
    int precondition, paralength, iterNum;

    // filename
    char filename[64];
    char *srcfile, *recordfile, *vpfile, *rhofile;

    // time recording
    double start, end, totalTime, seconds;
    int hours, minutes;
    
    // checkpoint memory cost (unit: GB)
    float memCost;

    Record2D *record2d_obs, *record2d_pre0, *record2d_pre1;

    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose))
        verbose = 0;

    // Get model grid parameter
    if (!getparint("nx", &nx))
        err("<error>: must specify nx!");
    if (!getparint("nz", &nz))
        err("<error>: must specify nz!");
    if (!getparint("nt", &nt))
        err("<error>: must specify nt!");
    if (!getparfloat("dx", &dx))
        err("<error>: must specify dx!");
    if (!getparfloat("dz", &dz))
        err("<error>: must specify dz!");
    if (!getparfloat("dt", &dt))
        err("<error>: must specify dt!");
    
    // Get boundary parameter
    if (!getparint("freesurface", &freesurface))
        freesurface = 0;
    if (!getparint("pmlThick", &pmlThick))
        pmlThick = 30;

    // Get iteration values
    if ((paralength = countparval("iteration")) != 0)
    {
        iteration = ealloc1int(paralength);
        getparint("iteration", iteration);
    }
    else
        err("<error> length of iteration is 0!");

    iterNum = 0;
    for (i = 0; i < paralength; i++)
        iterNum += iteration[i];

    // Get frequency values
    if ((paralength = countparval("freq")) != countparval("iteration"))
        err("<error> length of freq != length of iteration");
    else
    {
        freq = ealloc1float(paralength);
        getparfloat("freq", freq);
    }

    // Get band values
    if ((paralength = countparval("filter")) != countparval("iteration"))
        err("<error> length of iteration != length of filter");
    else
    {
        filter = ealloc1int(paralength);
        getparint("filter", filter);
    }

    // Get smooth parameter values
    if ((paralength = countparval("rx")) != countparval("iteration"))
        err("<error> length of iteration != length of rx");
    else
    {
        rx = ealloc1int(paralength);
        getparint("rx", rx);
    }

    if ((paralength = countparval("rz")) != countparval("iteration"))
        err("<error> length of iteration != length of rz");
    else
    {
        rz = ealloc1int(paralength);
        getparint("rz", rz);
    }

    if (!getparint("Gardner", &Gardner))
        Gardner = 1;
    if (!getparfloat("vpmax", &vpmax))
        vpmax = 10000.0f;
    if (!getparfloat("vpmin", &vpmin))
        vpmin = 0.001;
    if (!getparfloat("rhomax", &rhomax))
        rhomax = 10000.0f;
    if (!getparfloat("rhomin", &rhomin))
        rhomin = 0.001;
    if (!getparstring("srcfile", &srcfile))
    {
        err("<error>: didn't specify source signature file!");
    }
    if (!getparstring("recordfile", &recordfile))
    {
        err("<error>: didn't specify recordfile!");
    }
    if (!getparstring("vpfile", &vpfile))
    {
        err("<error>: didn't specify data file of p-velocity model!");
    }
    if (Gardner != 0)
    {
        rhofile = NULL;
        fprintf(stderr,"<warning>: Gardner density is applied!");
    }
    else
    {
        if (!getparstring("rhofile", &rhofile))
        {
            err("<error>: didn't specify data file of density model!");
        }
    }

    // Get water layer infomation
    if (!getparint("waterdepth", &waterdepth))
    {
        fprintf(stderr,"<warning>: didn't specify water depth! Set it as 0m.");
        waterdepth = 0;
    }
    if (!getparfloat("vp_water", &vp_water))
    {
        fprintf(stderr,"<warning>: didn't specify water velocity! Set it as 1500 m/s.");
        vp_water = 1500.0f;
    }
    if (!getparfloat("rho_water", &rho_water))
    {
        fprintf(stderr,"<warning>: didn't specify water density! Set it as 1000 kg/m^3.");
        rho_water = 1000.0f;
    }

    // Get shot number
    if (!getparint("Nshot", &Nshot))
    {
        fprintf(stderr,"<error>: didn't specify the number of shots!");
    }

    if (!getparint("precondition", &precondition))
    {
        fprintf(stderr,"<warning>: didn't specify precondition! Set it as 1.");
        precondition = 1;
    }

    /* Print all input parameters to stdout */
    fprintf(stdout, "\n==================== Input Parameters ====================\n");

    /* Basic model setup */
    fprintf(stdout, "Model grid: nx = %d, nz = %d, nt = %d\n", nx, nz, nt);
    fprintf(stdout, "Sampling:   dx = %.3f m, dz = %.3f m, dt = %.6f s\n", dx, dz, dt);
    fprintf(stdout, "Boundary:   freesurface = %d, pmlThick = %d\n", freesurface, pmlThick);

    /* Inversion control parameters */
    fprintf(stdout, "\n------------ Frequency Multi-scale Parameters ------------\n");
    fprintf(stdout, "Number of frequency bands = %d\n", paralength);
    fprintf(stdout, "Group\tFrequency\tIteration\tfilter\tsmooth radius\n");
    for(int i=0; i<paralength; i++)
    {
        fprintf(stdout, "%3d\t%6.2f\t\t%5d\t\t%3d\t(%d,%d)\n",i+1,freq[i],iteration[i],filter[i],rx[i], rz[i]);
    }

    // fprintf(stdout, "Iteration counts per group: ");
    // for (i = 0; i < paralength; i++)
    //     fprintf(stdout, "%d ", iteration[i]);
    // fprintf(stdout, "\nTotal iterations = %d\n", iterNum);

    // fprintf(stdout, "Filtering frequency (Hz):     ");
    // for (i = 0; i < paralength; i++)
    //     fprintf(stdout, "%.2f ", freq[i]);
    // fprintf(stdout, "\n");

    // fprintf(stdout, "Filter type per stage:      ");
    // for (i = 0; i < paralength; i++)
    //     fprintf(stdout, "%d ", filter[i]);
    // fprintf(stdout, "\n");

    // fprintf(stdout, "Smoothing radius (rx, rz):  ");
    // for (i = 0; i < paralength; i++)
    //     fprintf(stdout, "(%d,%d) ", rx[i], rz[i]);
    // fprintf(stdout, "\n");

    /* Physical property limits */
    fprintf(stdout, "\n----------- Model Parameter Bounds -----------\n");
    fprintf(stdout, "Velocity range: vpmin = %.3f, vpmax = %.3f\n", vpmin, vpmax);
    fprintf(stdout, "Density range:  rhomin = %.3f, rhomax = %.3f\n", rhomin, rhomax);
    fprintf(stdout, "Gardner relation: %s\n", (Gardner != 0) ? "enabled" : "disabled");

    /* File paths */
    fprintf(stdout, "\n----------------- Input Files ----------------\n");
    fprintf(stdout, "Velocity model file: %s\n", vpfile);
    if (Gardner == 0)
        fprintf(stdout, "Density model file:  %s\n", rhofile);
    fprintf(stdout, "Source wavelet file: %s\n", srcfile);
    fprintf(stdout, "Record file:         %s\n", recordfile);

    /* Water layer */
    fprintf(stdout, "\n----------------- Water Layer ----------------\n");
    fprintf(stdout, "Depth = %d m\n", waterdepth);
    fprintf(stdout, "P-wave velocity = %.2f m/s\n", vp_water);
    fprintf(stdout, "Density = %.2f kg/m^3\n", rho_water);

    /* Source-receiver setup */
    fprintf(stdout, "\n------------ Acquisition Geometry ------------\n");
    fprintf(stdout, "Number of shots = %d\n", Nshot);

    /* Preconditioning */
    fprintf(stdout, "\n------------- Numerical Settings -------------\n");
    fprintf(stdout, "Gradient precondition = %d\n", precondition);

    /* Verbose mode */
    fprintf(stdout, "\n----------------- Info Output -----------------\n");
    fprintf(stdout, "Verbose mode = %d\n", verbose);

    fprintf(stdout, "==========================================================\n\n");
    fflush(stdout);

    // allocate memory
    src = alloc1float(nt);
    vp = alloc2float(nz, nx);
    rho = alloc2float(nz, nx);
    slow = alloc2float(nz, nx);
    grad = alloc2float(nz, nx);
    perturb = alloc2float(nz, nx);
    fun1 = alloc1float(iterNum);
    fun2 = alloc1float(iterNum);
    timeCost = alloc1float(iterNum);

    memset(timeCost, 0, FSIZE * iterNum);

    // input file
    read_bin1d(srcfile, nt, src);
    read_bin2d(vpfile, nz, nx, vp);
    if (rhofile == NULL || Gardner != 0)
        vp2density(vp[0], rho[0], nz * nx);
    else
        read_bin2d(rhofile, nz, nx, rho);

    // fixed the set_waterlayer_value
    set_waterlayer_value(vp, vp_water, waterdepth, nx, nz);
    set_waterlayer_value(rho, rho_water, waterdepth, nx, nz);

    // calculate slowness
    v2slow(vp[0], slow[0], nz * nx);
    
    // read shot gather and source from file
    record2d_obs = input_record2d(recordfile, srcfile, nt, dt, Nshot);

    if (verbose > 0)
    {
        // save vp file
        sprintf(filename, "vp_init_%d_%d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
        write_bin(filename, vp[0], nz * nx);
        // save density file
        sprintf(filename, "rho_init_%d_%d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
        write_bin(filename, rho[0], nz * nx);
        // save source file
        sprintf(filename, "wavelet_%d_%.2fms.bin", nt, dt * 1000);
        write_bin(filename, src, nt);
    }

    // define griding observation
    create_grid_observe2d(record2d_obs, Nshot, dx, dz, nx);

    // copy header of observed data to predicted data structure
    record2d_pre0 = copy_record2d_header(record2d_obs, Nshot);
    record2d_pre1 = copy_record2d_header(record2d_obs, Nshot);
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        record2d_pre0[ishot].traces_p = alloc2float(record2d_obs[ishot].ns, record2d_obs[ishot].ntr);
        record2d_pre1[ishot].traces_p = alloc2float(record2d_obs[ishot].ns, record2d_obs[ishot].ntr);
    }
    
    // calculate max memory cost of checkpoint method
    memCost = calculate_memory_cost_checkpoint(nt, nx + 2 * pmlThick, nz + 2 * pmlThick);
    fprintf(stdout, "Estimated checkpoint memory cost: %.3f GB\n", memCost);
    fprintf(stdout, "Max threads = %d\n", omp_get_max_threads());
    // inversion loop
    fprintf(stdout, "Inversion starts...\n");
    fflush(stdout);
    
    iter = 0;
    for (k = 0; k < paralength; k++)
    {
        for (i = 0; i < iteration[k]; i++)
        {
            start = omp_get_wtime();

            // GRADIENT
            iso_acoustic2d_propagation(
                vp, rho, record2d_obs,
                Nshot, nz, nx, dx, dz, dt, nt,
                pmlThick, freesurface, freq[k], filter[k], precondition, 1,
                iter + 1, verbose,
                record2d_pre0, grad);

            smooth2D(grad, perturb, nx, nz, rz[k], 1);
            smooth2D(perturb, grad, nx, nz, rx[k], 0);

            // STEP-LENGTH
            fun1[iter] = calculate_objective_value(record2d_pre0, Nshot);
            scalar = compute_model_perturbation(0.01, slow[0], grad[0], perturb[0], nz * nx);

            set_waterlayer_value(grad, 0.0f, waterdepth, nx, nz);
            
            // Perturb Model 
            vadd(slow[0], perturb[0], vp[0], nx * nz);

            slow2v(vp[0], vp[0], nz * nx);
            set_waterlayer_value(vp, vp_water, waterdepth, nx, nz);
            clip_model(vp[0], nz * nx, vpmin, vpmax);

            if (Gardner == 1)
            {
                vp2density(vp[0], rho[0], nz * nx);
                set_waterlayer_value(rho, rho_water, waterdepth, nx, nz);
                clip_model(rho[0], nz * nx, rhomin, rhomax);
            }

            iso_acoustic2d_propagation(
                vp, rho, record2d_obs,
                Nshot, nz, nx, dx, dz, dt, nt,
                pmlThick, freesurface, freq[k], filter[k], precondition, 2,
                iter + 1, verbose,
                record2d_pre1, NULL);
                
            fun2[iter] = calculate_objective_value(record2d_pre1, Nshot);
            steplen = calculate_step_length(record2d_pre0, record2d_pre1, Nshot);

            // UPDATE Model
            nvmul(perturb[0], steplen, perturb[0], nx * nz);
            vadd(slow[0], perturb[0], slow[0], nx * nz);
            slow2v(slow[0], vp[0], nz * nx);
            set_waterlayer_value(vp, vp_water, waterdepth, nx, nz);
            clip_model(vp[0], nz * nx, vpmin, vpmax);

            if (Gardner == 1)
            {
                vp2density(vp[0], rho[0], nz * nx);
                set_waterlayer_value(rho, rho_water, waterdepth, nx, nz);
                clip_model(rho[0], nz * nx, rhomin, rhomax);
            }
            
            v2slow(vp[0], slow[0], nz * nx);

            // output
            sprintf(filename, "grad_iter%03d_%05d_%05d_%.2fm_%.2fm.bin", iter + 1, nz, nx, dz, dx);
            write_bin(filename, grad[0], nz * nx);
            sprintf(filename, "vp_iter%03d_%05d_%05d_%.2fm_%.2fm.bin",   iter + 1, nz, nx, dz, dx);
            write_bin(filename, vp[0], nz * nx);
            sprintf(filename, "rho_iter%03d_%05d_%05d_%.2fm_%.2fm.bin",  iter + 1, nz, nx, dz, dx);
            write_bin(filename, rho[0], nz * nx);

            end = omp_get_wtime();
            timeCost[iter] = end - start;
            fprintf(stdout, "Iteration %03d: freq = %4.1fHz, steplength = %e, fun1 = %e, fun2 = %e, takes %f seconds.\n",
                    iter + 1, freq[k], steplen, fun1[iter], fun2[iter], timeCost[iter]);
            fflush(stdout);
            iter++;
        }
    }

    // calculate final predicted data
    iso_acoustic2d_propagation(
        vp, rho, record2d_obs,
        Nshot, nz, nx, dx, dz, dt, nt,
        pmlThick, freesurface, 0, 0, 0, 0,
        0, verbose,
        record2d_pre0, NULL);

    // final output
    write_acoustic2d_traces("shots_p_predicted", record2d_pre0, Nshot);
    sprintf(filename, "objective_function_%03d.bin", iterNum);
    write_bin(filename, fun1, iterNum);
    sprintf(filename, "vp_inversion_%05d_%05d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
    write_bin(filename, vp[0], nz * nx);
    sprintf(filename, "rho_inversion_%05d_%05d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
    write_bin(filename, rho[0], nz * nx);

    // compute hours, minutes, seconds
    totalTime = sasum(iterNum, timeCost, 0);
    hours = (int)(totalTime / 3600);
    minutes = (int)((totalTime - hours * 3600) / 60);
    seconds = totalTime - hours * 3600 - minutes * 60;

    fprintf(stdout, "Inversion Takes: %d h %d min %.3f sec\n", hours, minutes, seconds);
    fflush(stdout);

    free1int(iteration);
    free1float(freq);
    free1int(filter);
    free1int(rx);
    free1int(rz);
    free1float(src);
    free2float(vp);
    free2float(rho);
    free2float(slow);
    free2float(grad);
    free2float(perturb);
    free1float(fun1);
    free1float(fun2);
    free1float(timeCost);
    free_record2d(record2d_obs, Nshot);
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        free2float(record2d_pre0[ishot].traces_p);
        free2float(record2d_pre1[ishot].traces_p);
    }

    return (CWP_Exit());
}