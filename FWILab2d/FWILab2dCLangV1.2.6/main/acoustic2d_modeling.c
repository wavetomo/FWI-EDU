/*
 * @Author: QuYang Chen
 * @Date: 2024-11-25 16:01:32
 * @LastEditors: Quyang Chen
 * @LastEditTime: 2025-12-23 16:31:54
 * @FilePath: /FWILab2dCLangV1.2.5/FWILab2dCLangV1.2.5/main/acoustic2d_modeling.c
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
#include "cJSON.h"

char *sdoc[] =  {
    "acoustic2d_modeling - 2D acoustic wave propagation modeling with optional free surface and MEAL boundaries.",
    "",
    "Usage:",
    "    acoustic2d_modeling [parameters]",
    "",
    "Required Parameters:",
    "    nx=             Number of grid points in the horizontal (x) direction",
    "    nz=             Number of grid points in the vertical (z) direction",
    "    nt=             Number of time steps",
    "    dx=             Grid spacing in x direction (unit: m)",
    "    dz=             Grid spacing in z direction (unit: m)",
    "    dt=             Time step (unit: s)",
    "    srcfile=        Binary file containing source wavelet (data type: float)",
    "    vpfile=         Binary file containing P-wave velocity model (data type: float)",
    "",
    "Optional Parameters:",
    "    rhofile=        Binary file containing density model (data type: float)",
    "                    [default: disabled if Gardner=1]",
    "    Gardner=        Apply Gardner's density relation (1=yes, 0=no)",
    "                    [default: 1]",
    "",
    "    freesurface=    Apply free-surface boundary condition (1=yes, 0=no)",
    "                    [default: 0]",
    "    pmlThick=       Thickness of the MEAL boundary (unit: cell)",
    "                    [default: 30]",
    "",
    "    waterdepth=     Depth of water layer (unit: m)",
    "                    [default: 0]",
    "    vp_water=       P-wave velocity in water (unit: m/s)",
    "                    [default: 1500.0]",
    "    rho_water=      Density in water (unit: kg/m^3)",
    "                    [default: 1000.0]",
    "",
    "Survey Geometry:",
    "    Geometry can be read from a JSON file or directly from parameters.",
    "",
    "    jsonfile=       JSON file containing acquisition geometry",
    "                    [read from parameters if not provided]",
    "",
    "    If jsonfile is NOT specified, the following parameters will be used:",
    "        Nshot=      Number of shots [default: 99999]",
    "        fsx=        First shot's x-position (unit: cell) [default: 1]",
    "        ds=         Shot interval (unit: cell) [default: 1]",
    "        sdepth=     Shot depth (unit: cell) [default: 2]",
    "        offsmin=    Minimum receiver offset (unit: cell) [default: 0]",
    "        offsmax=    Maximum receiver offset (unit: cell) [default: 99999]",
    "        dtr=        Receiver interval (unit: cell) [default: 1]",
    "        gdepth=     Receiver depth (unit: cell) [default: 2]",
    "",
    "Other Parameters:",
    "    sampleRate=     Output sampling rate of traces [default: 1]",
    "                    e.g., =2: record sampling rate is 2*dt",
    "    verbose=        verbosity level",
    "                    0 = dispaly minimal information",
    "                    1 = display more information without wavefields, good for QC",
    "                    2 = display all information with wavefields, good for code development",
    "",
    "Notes:",
    "    - When Gardner=1, the density model will be automatically derived from velocity using Gardner's equation.",
    "    - If the JSON file is provided, acquisition geometry will be completely read from it.",
    "    - All grid-based input binary files should be stored in C float (32-bit) format with column-major order (Seismic Unix style).",
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
    float *x, *z;

    int ix, iz, it, ishot;

    // source
    int Nshot;
    float fsx, ds, sdepth;
    float *src, *sx, *sz;

    // receiver
    int sampleRate;
    float offsmin, offsmax, dtr, gdepth;

    // model
    int Gardner, waterdepth;
    float vp_water, rho_water;
    float **vp, **rho;

    // filename
    char filename[64];
    char *srcfile, *vpfile, *rhofile;
    char *jsonfile, *jsonstr;
    int isReadJson;
    cJSON *jsondata, *records;

    // openmp
    double start, end, totalTime, seconds;
    int hours, minutes;

    Record2D *record2d, *record2d_pre;

    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose))
        verbose = 0;

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
    if (!getparint("freesurface", &freesurface))
        freesurface = 0;
    if (!getparint("pmlThick", &pmlThick))
        pmlThick = 30;

    if (!getparint("Gardner", &Gardner))
        Gardner = 1;
    if (!getparstring("srcfile", &srcfile))
    {
        err("<error>: didn't specify source signature file!");
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

    //
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

    // read survey geometry from JSON file or parameters
    // isReadJson = 0: read from parameters
    // isReadJson = 1: read from JSON file
    if (!getparstring("jsonfile", &jsonfile) || jsonfile[0] == '\0')
    {
        fprintf(stderr,"<warning>: didn't specify acquisition geometry JSON file. It will be read from the parameters.");
        isReadJson = 0;
        // source
        if (!getparint("Nshot", &Nshot))
        {
            fprintf(stderr,"<warning>: didn't specify the number of shots! Set it as 99999.");
            Nshot = 99999;
        }
        if (!getparfloat("fsx", &fsx))
        {
            fprintf(stderr,"<warning>: didn't specify first shot's x-position! Set it at No.1 cell.");
            fsx = 1.0f;
        }
        if (!getparfloat("ds", &ds))
        {
            fprintf(stderr,"<warning>: didn't specify the interval of shot! Set it as 1 cell.");
            ds = 1.0f;
        }
        if (!getparfloat("sdepth", &sdepth))
        {
            fprintf(stderr,"<warning>: didn't specify depth of shots! Set it at No.2 cell.");
            sdepth = 2.0f;
        }

        // receiver
        if (!getparfloat("offsmin", &offsmin))
        {
            fprintf(stderr,"<warning>: didn't specify the minimum offset of receiver! Set it as 0 cell.");
            offsmin = 0.0f;
        }
        if (!getparfloat("offsmax", &offsmax))
        {
            fprintf(stderr,"<warning>: didn't specify the minimum offset of receiver! Set it as 99999 cell.");
            offsmax = 99999.0f;
        }
        if (!getparfloat("dtr", &dtr))
        {
            fprintf(stderr,"<warning>: didn't specify interval of trace! Set it as 1 cell.");
            dtr = 1.0f;
        }
        if (!getparfloat("gdepth", &gdepth))
        {
            fprintf(stderr,"<warning>: didn't specify depth of receiver! Set it at No.2 cell.");
            gdepth = 2.0f;
        }
        
        // tramsform coordinate 
        fsx *=  dx;
        ds *= dx;
        sdepth *= dz;
        offsmin *= dx;
        offsmax *= dx;
        dtr *=  dx;
        gdepth *= dz;
    }
    else
    {
        isReadJson = 1;
        if (!getparint("Nshot", &Nshot))
        {
            err("<error>: didn't specify the number of shots!");
        }
    }

    if (!getparint("sampleRate", &sampleRate))
    {
        fprintf(stderr,"<warning>: didn't specify the time sample rate of trace! Set it as 1.");
        sampleRate = 1;
    }

    /* print input parameters to stdout */
    fprintf(stdout, "\n==================== Input Parameters ====================\n");
    fprintf(stdout, "Model grid: nx = %d, nz = %d, nt = %d\n", nx, nz, nt);
    fprintf(stdout, "Sampling:   dx = %.3f m, dz = %.3f m, dt = %.6f s\n", dx, dz, dt);
    fprintf(stdout, "Boundary:   freesurface = %d, pmlThick = %d\n", freesurface, pmlThick);

    fprintf(stdout, "Density:    Gardner = %d\n", Gardner);
    if (Gardner != 0)
        fprintf(stdout, "             -> Using Gardner relation for density.\n");
    else
        fprintf(stdout, "             -> Density file: %s\n", rhofile);

    fprintf(stdout, "Velocity:   vpfile = %s\n", vpfile);
    fprintf(stdout, "Source:     srcfile = %s\n", srcfile);

    fprintf(stdout, "Water layer:\n");
    fprintf(stdout, "             depth = %d m\n", waterdepth);
    fprintf(stdout, "             vp_water = %.2f m/s\n", vp_water);
    fprintf(stdout, "             rho_water = %.2f kg/m^3\n", rho_water);

    if (isReadJson == 1)
    {
        fprintf(stdout, "Acquisition geometry: read from JSON file = %s\n", jsonfile);
        fprintf(stdout, "Number of shots = %d\n", Nshot);
    }
    else
    {
        fprintf(stdout, "Acquisition geometry: read from parameters\n");
        fprintf(stdout, "Number of shots = %d\n", Nshot);
        fprintf(stdout, "Shot geometry: fsx = %.3f m, ds = %.3f m, sdepth = %.3f m\n", fsx, ds, sdepth);
        fprintf(stdout, "Receiver geometry: offsmin = %.3f m, offsmax = %.3f m, dtr = %.3f m, gdepth = %.3f m\n",
                offsmin, offsmax, dtr, gdepth);
    }

    fprintf(stdout, "Trace sample rate = %d\n", sampleRate);
    fprintf(stdout, "Verbose mode = %d\n", verbose);
    fprintf(stdout, "==========================================================\n\n");
    fflush(stdout);

    // allocate memory
    x = alloc1float(nx);
    z = alloc1float(nz);
    src = alloc1float(nt);
    vp = alloc2float(nz, nx);
    rho = alloc2float(nz, nx);

    // coordinate define
    for (ix = 0; ix < nx; ix++)
        x[ix] = (float)ix * dx;
    for (iz = 0; iz < nz; iz++)
        z[iz] = (float)iz * dz;

    // input file
    read_bin1d(srcfile, nt, src);
    read_bin2d(vpfile, nz, nx, vp);
    if (rhofile == NULL || Gardner != 0)
        vp2density(vp[0], rho[0], nz * nx);
    else
        read_bin2d(rhofile, nz, nx, rho);

    // fixed the waterlayer
    set_waterlayer_value(vp, vp_water, waterdepth, nx ,nz);
    set_waterlayer_value(rho, rho_water, waterdepth, nx ,nz);

    if (verbose > 0)
    {
        // save vp file
        sprintf(filename, "vp_%d_%d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
        write_bin(filename, vp[0], nz * nx);
        // save density file
        sprintf(filename, "rho_%d_%d_%.2fm_%.2fm.bin", nz, nx, dz, dx);
        write_bin(filename, rho[0], nz * nx);
        sprintf(filename, "wavelet_%d_%.2fms.bin", nt, dt * 1000);
        write_bin(filename, src, nt);
    }

    if (isReadJson == 1)
    {
        // read from json file
        jsonstr = read_string_file(jsonfile);
        jsondata = cJSON_Parse(jsonstr);
        if (!jsondata)
        {
            free(jsonstr);
            err("<error>: Parse JSON from file (%s) Error!", filename);
        }
        records = cJSON_GetObjectItem(jsondata, "record");
        record2d = alloc_record2d(Nshot);
        init_observe2d_from_JSON(record2d, records, Nshot, x, nx, z, nz,
                                dt, nt, src, sampleRate);
    }
    else
    {
        // read from input parameters
        Nshot = generate_shot_positions(x, nx, z, nz, fsx, ds, sdepth, Nshot, &sx, &sz);
        record2d = alloc_record2d(Nshot);
        init_observe2d(record2d, Nshot, x, nx, z, nz,
                      dt, nt, sx, sz, src, offsmin, offsmax, dtr, gdepth, sampleRate);
    }

    // define griding observation
    create_grid_observe2d(record2d, Nshot, dx, dz, nx);

    // copy header of observed data to predicted data structure
    record2d_pre = copy_record2d_header(record2d, Nshot);
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        record2d_pre[ishot].traces_p = alloc2float(record2d[ishot].ns, record2d[ishot].ntr);
    }

    // start modeling
    fprintf(stdout, "Max threads = %d\n", omp_get_max_threads());
    fprintf(stdout, "Modeling starts...\n");
    fflush(stdout);
    
    start = omp_get_wtime();

    iso_acoustic2d_propagation(
        vp, rho, record2d,
        Nshot, nz, nx, dx, dz, dt, nt,
        pmlThick, freesurface, 0, 0, 0, 0,
        0, verbose,
        record2d_pre, NULL);

    // gather output
    write_acoustic2d_traces("shots_p", record2d_pre, Nshot);

    // compute running time
    end = omp_get_wtime();
    totalTime = end - start;

    // compute hours, minutes, seconds
    hours = (int)(totalTime / 3600);
    minutes = (int)((totalTime - hours * 3600) / 60);
    seconds = totalTime - hours * 3600 - minutes * 60;

    fprintf(stdout, "Modeling takes: %d h %d min %.3f sec\n", hours, minutes, seconds);
    fflush(stdout);

    // free memory
    free1float(x);
    free1float(z);
    free1float(src);
    free2float(vp);
    free2float(rho);
    free_record2d(record2d, Nshot);
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        free2float(record2d_pre[ishot].traces_p);
    }

    if (isReadJson == 1)
    {
        free(jsonstr);
        cJSON_Delete(jsondata);
    }
    else
    {
        free1float(sx);
        free1float(sz);
    }
    return (CWP_Exit());
}
