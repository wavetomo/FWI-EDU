/*
 * @Author: QuYang Chen
 * @Date: 2025-09-22 23:22:49
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-11-03 14:36:56
 * @FilePath: /FWILab2dCLangV1.2.3/src/observation.c
 * @Description:
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#include "observation.h"


Record2D *alloc_record2d(int Nshot)
{
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
    
    if (Nshot <= 0)
    {
        err("alloc_record2d: Nshot must be positive!\n");
    }

    Record2D *record2d = (Record2D *)calloc(Nshot, sizeof(Record2D));
    if (!record2d)
    {
        err("alloc_record2d: Failed to allocate Record2D array!\n");
    }

    for (int ishot = 0; ishot < Nshot; ishot++)
    {
        record2d[ishot].shotID = ishot;
        record2d[ishot].src = NULL;
        record2d[ishot].gx = NULL;
        record2d[ishot].gz = NULL;
        record2d[ishot].igx = NULL;
        record2d[ishot].igz = NULL;
        record2d[ishot].traces_vx = NULL;
        record2d[ishot].traces_vz = NULL;
        record2d[ishot].traces_p = NULL;
    }

    return record2d;
}


void free_record2d(Record2D *record2d, int Nshot)
{
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
    
    if (!record2d)
        return;

    for (int ishot = 0; ishot < Nshot; ishot++)
    {

        free1float(record2d[ishot].src);
        free1float(record2d[ishot].gx);
        free1float(record2d[ishot].gz);
        free1int(record2d[ishot].igx);
        free1int(record2d[ishot].igz);
        
        if (record2d[ishot].traces_p)
        {
            free2float(record2d[ishot].traces_p);
        }
        // if (record2d[ishot].traces_vx) {
        //     free2float(record2d[ishot].traces_vx);
        // }
        // if (record2d[ishot].traces_vz) {
        //     free2float(record2d[ishot].traces_vz);
        // }
    }

    free(record2d);
}

int generate_shot_positions(
    const float *x, int nx,
    const float *z, int nz,
    float fsx, float ds, float sdepth,
    int Nshot_in,
    float **sx, float **sz)
{
    /* 
    generate_shot_positions: generate shot coordinates and determine number of shots.

    Inputs:
        x        - array of x-coordinates in model domain
        nx       - length of x array
        z        - array of z-coordinates in model domain
        nz       - length of z array
        fsx      - x-coordinate of the first shot
        ds       - interval between adjacent shots in x-direction
        sdepth   - fixed z-coordinate for all shots
        Nshot_in - desired number of shots

    Outputs:
        sx       - pointer to array pointer of x-positions of shots (malloc inside)
        sz       - pointer to array pointer of z-positions of shots (malloc inside)

    retrun:
        Nshot    - actual number of shots used (may be reduced)

    Notes:
        - Caller must free sx and sz after use.
        - If fsx/sdepth out of bounds, the function will call err and exit.
    */
    
    int ishot, maxShots, Nshot;

    if (nx < 1 || nz < 1)
    {
        err("Model domain arrays x or z are empty!\n");
    }

    /* Check if fsx is within x bounds */
    if (fsx < x[0] || fsx > x[nx - 1])
    {
        err("x-coordinate of first shot is out of range!\n");
    }

    /* Check if sdepth is within z bounds */
    if (sdepth < z[0] || sdepth > z[nz - 1])
    {
        err("z-coordinate of shots is out of range!\n");
    }

    /* Compute maximum possible number of shots from fsx to x[nx-1] */
    maxShots = (int)((x[nx - 1] - fsx) / ds) + 1;

    /* Use the smaller of user-requested Nshot_in and available maxShots */
    Nshot = (Nshot_in < maxShots) ? Nshot_in : maxShots;

    *sx = alloc1float(Nshot);
    *sz = alloc1float(Nshot);

    /* Fill shot positions */
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        (*sx)[ishot] = fsx + ishot * ds;
        (*sz)[ishot] = sdepth;
    }

    return Nshot;
}


Record2D *copy_record2d_header(const Record2D *record2d_src, int Nshot)
{
    /* 
    copy_record2d_header: copy only the header information (not wavefield data and traces_vx, traces_vz, and traces_p)

    Inputs:
        record2d_src - pointer to source Record2D array
        Nshot        - number of shots in the array

    Outputs:
        return       - pointer to newly allocated Record2D array with copied headers 
    */
    
    // Allocate new array of Record2D
    Record2D *record2d_dest = alloc_record2d(Nshot);

    for (int ishot = 0; ishot < Nshot; ishot++)
    {
        // Copy scalar values
        record2d_dest[ishot].shotID = record2d_src[ishot].shotID;
        record2d_dest[ishot].sx = record2d_src[ishot].sx;
        record2d_dest[ishot].sz = record2d_src[ishot].sz;
        record2d_dest[ishot].isx = record2d_src[ishot].isx;
        record2d_dest[ishot].isz = record2d_src[ishot].isz;
        record2d_dest[ishot].fpeak = record2d_src[ishot].fpeak;

        // Shallow copy for pointer src
        record2d_dest[ishot].src = record2d_src[ishot].src;

        // Copy receiver geometry pointers
        record2d_dest[ishot].gx = record2d_src[ishot].gx;
        record2d_dest[ishot].gz = record2d_src[ishot].gz;
        record2d_dest[ishot].igx = record2d_src[ishot].igx;
        record2d_dest[ishot].igz = record2d_src[ishot].igz;

        // Sampling info
        record2d_dest[ishot].sampleRate = record2d_src[ishot].sampleRate;
        record2d_dest[ishot].ntr = record2d_src[ishot].ntr;
        record2d_dest[ishot].ns = record2d_src[ishot].ns;
        record2d_dest[ishot].dt = record2d_src[ishot].dt;
        record2d_dest[ishot].dtr = record2d_src[ishot].dtr;

        // Recording window
        record2d_dest[ishot].ix_min = record2d_src[ishot].ix_min;
        record2d_dest[ishot].ix_max = record2d_src[ishot].ix_max;

        // traces_vx/vz/p
        record2d_dest[ishot].traces_vx = NULL;
        record2d_dest[ishot].traces_vz = NULL;
        record2d_dest[ishot].traces_p = NULL;
    }

    return record2d_dest;
}


void init_observe2d(Record2D *record2d, int Nshot, float *x, int nx, float *z, int nz,
                   float dt, int nt, float *sx, float *sz, float *src,
                   float offsmin, float offsmax, float dtr, float gdepth, int sampleRate)
{
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
    
    int ishot, itrace, it, ioffs, noffs;
    float offset, gx_abs;
    float *gx_offset;

    // Generate receiver offsets relative to shot
    noffs = (offsmin == 0.0f) ? 2 * ((int)(offsmax - offsmin) / dtr + 1) - 1 : 2 * ((int)(offsmax - offsmin) / dtr + 1);
    gx_offset = alloc1float(noffs);
    ioffs = 0;
    for (offset = -offsmax; offset <= -offsmin; offset += dtr)
    {
        gx_offset[ioffs] = offset;
        ioffs++;
    }
    if (offsmin == 0.0f)
        ioffs--;
    for (offset = offsmin; offset <= offsmax; offset += dtr)
    {
        gx_offset[ioffs] = offset;
        ioffs++;
    }

    // Loop over shots
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        record2d[ishot].shotID = ishot + 1;
        record2d[ishot].sx = sx[ishot];
        record2d[ishot].sz = sz[ishot];

        // Compute absolute receiver positions
        record2d[ishot].ntr = 0;
        for (ioffs = 0; ioffs < noffs; ioffs++)
        {
            gx_abs = record2d[ishot].sx + gx_offset[ioffs];
            if (gx_abs >= x[0] && gx_abs <= x[nx - 1])
                record2d[ishot].ntr++;
        }

        if (record2d[ishot].ntr == 0)
            err("the number of trace in shot %d is 0", ishot, record2d[ishot].ntr);

        record2d[ishot].gx = alloc1float(record2d[ishot].ntr);
        record2d[ishot].gz = alloc1float(record2d[ishot].ntr);

        itrace = 0;
        for (int ioffs = 0; ioffs < noffs; ioffs++)
        {
            gx_abs = record2d[ishot].sx + gx_offset[ioffs];
            if (gx_abs >= x[0] && gx_abs <= x[nx - 1])
            {
                record2d[ishot].gx[itrace] = gx_abs;
                record2d[ishot].gz[itrace] = gdepth;
                itrace++;
            }
        }

        // Sampling info
        record2d[ishot].sampleRate = sampleRate;
        record2d[ishot].ns = (nt - 1) / sampleRate + 1;
        record2d[ishot].dt = dt * sampleRate;

        record2d[ishot].dtr = 0;
        for (itrace = record2d[ishot].ntr - 1; itrace > 0; itrace--)
        {
            record2d[ishot].dtr += (record2d[ishot].gx[itrace] - record2d[ishot].gx[itrace - 1]);
        }
        record2d[ishot].dtr /= record2d[ishot].ntr - 1;

        // trace alloc
        record2d[ishot].traces_p = alloc2float(record2d[ishot].ns, record2d[ishot].ntr);

        // Source info
        record2d[ishot].src = alloc1float(nt);
        for (it = 0; it < nt; it++)
            record2d[ishot].src[it] = src[it];

        record2d[ishot].fpeak = calculate_peak_frequency(src, nt, dt);
    }

    free1float(gx_offset);
}


void create_grid_observe2d(Record2D *record2d, int Nshot, float dx, float dz, int nx)
{
    /* 
    create_grid_observe2d: convert real-world coordinates to grid indices

    Input:
        record2d    - array of Record2D
        Nshot       - number of shots
        dx          - spatial grid spacing in x
        dz          - spatial grid spacing in z
        nx          - total number of x-grid points

    Output:
        record2d    - array of Record2D 
    */
    
    int isx, isz, igx_min, igx_max, ix_min, ix_max;
    int ntr, ishot, itrace;

    for (ishot = 0; ishot < Nshot; ishot++)
    {

        // Validate receiver arrays
        if (!record2d[ishot].gx || !record2d[ishot].gz)
        {
            err("gx and gz must be allocated for shot %d\n", ishot);
        }

        ntr = record2d[ishot].ntr;
        if (ntr <= 0)
        {
            err("ntr must be positive for shot %d\n", ishot);
        }

        // Allocate receiver grid indices if not allocated
        if (!record2d[ishot].igx)
            record2d[ishot].igx = alloc1int(ntr);
        if (!record2d[ishot].igz)
            record2d[ishot].igz = alloc1int(ntr);

        // Convert source coordinates to global grid indices
        isx = (int)floorf(record2d[ishot].sx / dx);
        isz = (int)floorf(record2d[ishot].sz / dz);

        // Convert receiver coordinates to global grid indices
        igx_min = nx - 1; // initialize for min search
        igx_max = 0;      // initialize for max search

        for (itrace = 0; itrace < ntr; itrace++)
        {
            record2d[ishot].igx[itrace] = (int)floorf(record2d[ishot].gx[itrace] / dx);
            record2d[ishot].igz[itrace] = (int)floorf(record2d[ishot].gz[itrace] / dz);

            if (record2d[ishot].igx[itrace] < igx_min)
                igx_min = record2d[ishot].igx[itrace];
            if (record2d[ishot].igx[itrace] > igx_max)
                igx_max = record2d[ishot].igx[itrace];
        }

        // Define local simulation window with 30-point padding
        ix_min = igx_min - 30;
        if (ix_min < 0)
            ix_min = 0;

        ix_max = igx_max + 30;
        if (ix_max > nx - 1)
            ix_max = nx - 1;

        // Adjust source and receiver indices to local subgrid
        record2d[ishot].isx = isx - ix_min;
        record2d[ishot].isz = isz;

        for (itrace = 0; itrace < ntr; itrace++)
        {
            record2d[ishot].igx[itrace] -= ix_min; // local indices
            // igz stays global
        }

        // Store subgrid boundaries
        record2d[ishot].ix_min = ix_min;
        record2d[ishot].ix_max = ix_max;
    }
}

void init_observe2d_from_JSON(Record2D *record2d, cJSON *records, int Nshot, float *x, int nx, float *z, int nz,
                             float dt, int nt, float *src, int sampleRate)
{
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
    int ishot, itrace, it;
    int Nshot_actually, shotID, ntr, gx_len, gz_len;
    float sx, sz;
    cJSON *shot, *gx_array, *gz_array;

    Nshot_actually = cJSON_GetArraySize(records);
    if (Nshot_actually != Nshot)
    {
        err("<error>: The number of shots in JSON (%d) is not equal to Nshot (%d)!\n", Nshot_actually, Nshot);
    }
    // Loop over shots
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        // read shot info from cJSON type
        shot = cJSON_GetArrayItem(records, ishot);
        shotID = cJSON_GetObjectItem(shot, "shotID")->valueint;
        sx = (float)cJSON_GetObjectItem(shot, "sx")->valuedouble;
        sz = (float)cJSON_GetObjectItem(shot, "sz")->valuedouble;
        ntr = cJSON_GetObjectItem(shot, "ntr")->valueint;

        if (ntr <= 0)
        {
            err("<error>: the number of trace in shot %d is not more than 0", ishot, ntr);
        }
        if (sx < x[0] || sx > x[nx - 1])
        {
            err("<error>: x-coordinate of shot %d is out of range!\n", ishot);
        }
        if (sz < z[0] || sz > z[nz - 1])
        {
            err("<error>: z-coordinate of shot %d is out of range!\n", ishot);
        }
        // save shot info to record2d
        record2d[ishot].shotID = shotID;
        record2d[ishot].sx = sx;
        record2d[ishot].sz = sz;

        record2d[ishot].ntr = ntr;

        record2d[ishot].gx = alloc1float(record2d[ishot].ntr);
        record2d[ishot].gz = alloc1float(record2d[ishot].ntr);

        // read receiver coordinates from cJSON type
        gx_array = cJSON_GetObjectItem(shot, "gx");
        gx_len = cJSON_GetArraySize(gx_array);
        if (gx_len != ntr)
        {
            err("<error>: The number of gx (%d) is not equal to ntr (%d) in shot %d at JSON!\n", gx_len, ntr, ishot);
        }
        for (itrace = 0; itrace < gx_len; itrace++)
        {
            record2d[ishot].gx[itrace] = (float)cJSON_GetArrayItem(gx_array, itrace)->valuedouble;
            if (record2d[ishot].gx[itrace] < x[0] || record2d[ishot].gx[itrace] > x[nx - 1])
                err("<error>: x-coordinate of receiver %d in shot %d is out of range!\n", itrace, ishot);
        }

        gz_array = cJSON_GetObjectItem(shot, "gz");
        gz_len = cJSON_GetArraySize(gz_array);
        if (gz_len != ntr)
        {
            err("<error>: The number of gz (%d) is not equal to ntr (%d) in shot %d at JSON!\n", gz_len, ntr, ishot);
        }
        for (itrace = 0; itrace < gx_len; itrace++)
        {
            record2d[ishot].gz[itrace] = (float)cJSON_GetArrayItem(gz_array, itrace)->valuedouble;
            if (record2d[ishot].gz[itrace] < z[0] || record2d[ishot].gz[itrace] > z[nz - 1])
                err("<error>: z-coordinate of receiver %d in shot %d is out of range!\n", itrace, ishot);
        }

        // Sampling info
        record2d[ishot].sampleRate = sampleRate;
        record2d[ishot].ns = (nt - 1) / sampleRate + 1;
        record2d[ishot].dt = dt * sampleRate;

        record2d[ishot].dtr = 0;
        for (itrace = record2d[ishot].ntr - 1; itrace > 0; itrace--)
        {
            record2d[ishot].dtr += (record2d[ishot].gx[itrace] - record2d[ishot].gx[itrace - 1]);
        }
        record2d[ishot].dtr /= record2d[ishot].ntr - 1;

        // trace alloc
        record2d[ishot].traces_p = alloc2float(record2d[ishot].ns, record2d[ishot].ntr);

        // Source info
        record2d[ishot].src = alloc1float(nt);
        for (it = 0; it < nt; it++)
            record2d[ishot].src[it] = src[it];

        record2d[ishot].fpeak = calculate_peak_frequency(src, nt, dt);
    }
}