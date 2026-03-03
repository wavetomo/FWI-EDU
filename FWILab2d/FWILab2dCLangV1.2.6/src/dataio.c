/*
 * @Author: QuYang Chen
 * @Date: 2025-09-23 00:30:46
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-12-24 16:25:34
 * @FilePath: /FWILab2dCLangV1.2.5/src/dataio.c
 * @Description:
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */

#include "dataio.h"

int get_precision_bytes(const char *precision)
{
    /* 
    get_precision_bytes: return number of bytes per element for a given precision

    Input:
        precision - string describing data type, e.g., "float32", "double", "int32", etc.

    Output:
        return    - number of bytes per element (int)

    Notes:
        - Uses SU-style err() to report unsupported precision.
        - Accepted strings (case-insensitive):
        "float32", "single" -> 4
        "float64", "double" -> 8
        "int32" -> 4
        "int16" -> 2
        "uint8" -> 1 
    */ 
    
    if (precision == NULL)
    {
        err("get_precision_bytes: precision string is NULL\n");
    }

    if (strcmp(precision, "float32") == 0 || strcmp(precision, "single") == 0)
    {
        return 4;
    }
    else if (strcmp(precision, "float64") == 0 || strcmp(precision, "double") == 0)
    {
        return 8;
    }
    else if (strcmp(precision, "int32") == 0)
    {
        return 4;
    }
    else if (strcmp(precision, "int16") == 0)
    {
        return 2;
    }
    else if (strcmp(precision, "uint8") == 0)
    {
        return 1;
    }
    else
    {
        err("get_precision_bytes: unsupported precision: %s\n", precision);
        return 0; // will not reach here because err() exits
    }
}

long long get_file_size(const char *fileName)
{
    /* 
    get_file_size: get file size in bytes

    Input:
        fileName  - path to the file (const char *fileName)

    Output:
        return    - file size in bytes (long long), or -1LL if fseek/ftell fails
    */
    
    FILE *fp;
    long long fileSize;

    fp = fopen(fileName, "rb"); // open file in binary read mode

    if (!fp)
    {
        err("Error: Could not open file: %s\n", fileName);
        return -1LL;
    }

    // move file pointer to end
    if (fseek(fp, 0, SEEK_END) != 0)
    {
        fclose(fp);
        return -1LL;
    }

    // get position (file size)
    fileSize = ftell(fp);

    // close file
    fclose(fp);

    return fileSize;
}


long long count_file_elements(const char *fileName, int elementSize)
{
    /* 
    count_file_elements: calculate number of elements in a binary file

    Input:
        fileName     - path to the file (const char*)
        elementSize  - size of each element in bytes (int, must be > 0)

    Output:
        return       - number of elements (long long), or -1LL if an error occurs
    */
    
    if (elementSize <= 0)
    {
        err("Error: element size must be a positive integer.\n");
        return -1LL;
    }

    long long fileSize = get_file_size(fileName);
    if (fileSize < 0)
    {
        return -1LL;
    }

    long long numElements = fileSize / elementSize;

    if (fileSize % elementSize != 0)
    {
        err("Warning: file size is not a multiple of element size. Result may be truncated.\n");
        return -1LL;
    }

    return numElements;
}

void read_bin1d(const char *fileName, int n, float *data)
{
    /* 
    read_bin1d: read 1D float binary data from file

    Input:
        fileName  - path to binary file (const char*)
        n         - number of float elements to read (int)
                    if n <= 0, it will be inferred from file size

    Output:
        data      - pointer to allocated float array
    */

    // infer n if not provided
    if (n <= 0)
    {
        int elemSize = get_precision_bytes("float");
        n = count_file_elements(fileName, elemSize);
    }

    // open file
    FILE *fp = fopen(fileName, "rb");
    if (!fp)
    {
        err("read_bin1d: failed to open file (%s)\n", fileName);
    }

    // read data
    size_t count = fread(data, sizeof(float), n, fp);
    if (count != (size_t)n)
    {
        fclose(fp);
        err("read_bin1d: expected %d floats, but only read %zu\n", n, count);
    }

    fclose(fp);
}

void read_bin2d(const char *fileName, int n1, int n2, float **data)
{
    /* 
    read_bin2d: read 2D float binary data from file

    Input:
        fileName  - path to binary file (const char*)
        n1        - first dimension (e.g., number of cloumns) (int)
        n2        - second dimension (int)
                    if n2 <= 0, it will be inferred from file size

    Output:
        data      - pointer to allocated 2D float array (float**)
    */

    // infer n2 if not provided
    if (n2 <= 0)
    {
        int elemSize = get_precision_bytes("float");
        int totalElements = count_file_elements(fileName, elemSize);
        if (totalElements % n1 != 0)
        {
            err("Cannot evenly divide file into [%d, ?]. File size mismatch", n1);
        }
        n2 = totalElements / n1;
    }

    // open file
    FILE *fp = fopen(fileName, "rb");
    if (!fp)
    {
        err("read_bin2d: failed to open file: %s\n", fileName);
    }

    // read data cloumn by cloumn
    size_t count = fread(data[0], FSIZE, n1 * n2, fp);
    if (count != (size_t)(n1 * n2))
    {
        fclose(fp);
        err("read_bin2d: expected to read %d elements, but only read %zu\n", n1 * n2, count);
    }

    fclose(fp);
}

void read_bin3d(const char *fileName, int n1, int n2, int n3, float ***data)
{
    /* 
    read_bin3d: read 3D float binary data from file

    Input:
        fileName  - path to binary file (const char*)
        n1        - first dimension (e.g., depth) (int)
        n2        - second dimension (int)
        n3        - third dimension (int)
                    if n3 <= 0, it will be inferred from file size

    Output:
        data      - pointer to allocated 3D float array (float***)
                    Access: data[iz][iy][ix],
                    where iz=0..n1-1, ix=0..n2-1, iy=0..n3-1
    */

    // infer n3 if not provided
    if (n3 <= 0)
    {
        int elemSize = get_precision_bytes("float");
        int totalElements = count_file_elements(fileName, elemSize);
        if (totalElements % (n1 * n2) != 0)
        {
            err("Cannot evenly divide file into [%d, %d, ?]. File size mismatch", n1, n2);
        }
        n3 = totalElements / (n1 * n2);
    }

    // open file
    FILE *fp = fopen(fileName, "rb");
    if (!fp)
    {
        err("read_bin3d: failed to open file: %s\n", fileName);
    }

    // read data: block by block (each block is n1*n2)

    size_t count = fread(data[0][0], FSIZE, n1 * n2 * n3, fp);
    if (count != (size_t)(n1 * n2 * n3))
    {
        fclose(fp);
        err("read_bin3d: expected to read %d elements, but only read %zu\n",
            n1 * n2 * n3, count);
    }

    fclose(fp);
}

void write_bin(const char *fileName, float *data, int n)
{
    /* 
    write_bin: write 1D float array to binary file

    Input:
        fileName - path to output file (const char*)
        data     - pointer to float array
        n        - number of elements

    Output:
        None (writes binary data to disk as files)
    */
    
    FILE *fp = fopen(fileName, "wb");
    if (!fp)
    {
        err("write_bin: failed to open file: %s\n", fileName);
    }
    size_t count = fwrite(data, FSIZE, n, fp);
    if (count != (size_t)n)
    {
        fclose(fp);
        err("write_bin: expected to write %d elements, but only wrote %zu\n", n, count);
    }

    fclose(fp);
}

void write_acoustic2d_traces(char *prefix, Record2D *record2d, int Nshot)
{
    /* 
    write_acoustic2d_traces: write 2D shot gathers to disk

    Input:
        prefix      - file name prefix for output files (const char*)
        record2d    - array of shot gather structures containing traces and metadata (Record2D*)
        Nshot       - number of shots in the record2d array (int)

    Output:
        None (writes traces to disk as files)
    */
    
    FILE *fp = NULL;
    segy tr;
    char p_gatherName[256];
    int ishot, itrace;
    int ntr, nt;
    int temp;
    float sx, sz, gx, gz;
    float dt, scalel, scalco;
    
    ntr = 0;
#pragma omp parallel for reduction(+ : ntr) private(ishot)
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        ntr += record2d[ishot].ntr;
    }
    nt = record2d[0].ns;
    dt = record2d[0].dt * 1000;
    sprintf(p_gatherName, "%s_%d_%d_%d_%05.2fms.su", prefix, Nshot, nt, ntr, dt);
    fp = efopen(p_gatherName, "wb");
    memset(&tr, 0, sizeof(segy));
    tr.trid = 1;
    tr.tracl = 0;
    tr.tracr = 0;
    tr.scalel = -1000;
    tr.scalco = -1000;
    temp = SGN(-tr.scalel) * log10(abs((int)-tr.scalel));
    scalel = powf(10.0f, (float)temp);
    temp = SGN(-tr.scalco) * log10(abs((int)-tr.scalco));
    scalco = powf(10.0f, (float)temp);
    /* */
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        sx = record2d[ishot].sx;
        sz = record2d[ishot].sz;
        tr.fldr = record2d[ishot].shotID;
        tr.tracf = 0;
        tr.ntr = record2d[ishot].ntr;
        tr.ns = record2d[ishot].ns;
        tr.dt = record2d[ishot].dt * 1000000; /* micro-seconds */
        tr.d1 = record2d[ishot].dt;
        tr.d2 = record2d[ishot].dtr;
        tr.sdepth = sz * scalel;
        tr.sx = sx * scalco;
        
        for (itrace = 0; itrace < tr.ntr; itrace++)
        {
            tr.tracl++;
            tr.tracr++;
            tr.tracf++;
            gx = record2d[ishot].gx[itrace];
            gz = record2d[ishot].gz[itrace];    
            tr.gelev = gz * scalel;
            tr.gx = gx * scalco;
            tr.offset = tr.gx - tr.sx;
            tr.cdp = (tr.gx + tr.sx) / 2;
            memcpy(tr.data, record2d[ishot].traces_p[itrace], tr.ns * FSIZE);
            fputtr(fp, &tr);
        }
    }
    efclose(fp);
}

float scale_header(int keyword, int scale)
{
    /* 
    scale_header: apply scaling to a seismic trace header keyword

    Input:
        keyword - original header value (int)
        scale   - scaling factor applied to the header (int)

    Output:
        return  - scaled header value (float) 
    */
    
    int temp;
    float val;
    if (scale == 0)
    {
        err("<error> : scale must not be %d!", scale);
    }
    else
    {
        temp = SGN(scale) * log10(abs((int)scale));
        val = keyword * powf(10.0f, temp);
    }
    return val;
}

int get_shots_number(const char *filename)
{
    /* 
    get_shots_number: count the number of shots in a record file

    Input:
        filename - path to the record file (const char *)

    Output:
        return   - number of shots contained in the file (int) 
    */
    
    FILE *fp = NULL;
    segy tr;
    int numShots;
    memset(&tr, 0, sizeof(segy));
    fp = fopen(filename, "rb");
    if (fp == NULL)
    {
        err("<error>: failed to open file: %s!", filename);
        numShots = -1;
    }
    else
    {
        numShots = 0;
        while (fvgettr(fp, &tr))
        {
            if (tr.ntr == 0)
                err("<error>: the ntr of shot gather is 0!");
            numShots++;
            fseek(fp, (tr.ntr - 1) * (240 + tr.ns * FSIZE), SEEK_CUR);
        }
        fclose(fp);
    }
    return numShots;
}

void read_gather2d_from_sufile(Record2D *record2d, const char *recordfile,
                      int Nshot, int nt)
{
    /* 
    read_gather2d_from_sufile: read 2D gather file into Record2D structure

    Input:
        record2d   - pointer to Record2D array to store the gather data
        recordfile - path to the gather file (const char *)
        Nshot      - number of shots to read
        nt         - number of time samples per trace

    Output:
        record2d   - filled with receiver positions and trace data 
    */

    FILE *fp;
    segy tr;
    int ishot, itrace;
    double dtr;

    if (recordfile != NULL)
    {
        fp = fopen(recordfile, "rb");
        if (fp == NULL)
        {
            err("<error>: fail to open recordfile: %s!", recordfile);
        }
        else
        {
            for (ishot = 0; ishot < Nshot; ishot++)
            {
                fvgettr(fp, &tr);
                /* */
                if (tr.ntr <= 0)
                    err("<error>: the number of traces in shot%d is %d.", ishot + 1, tr.ntr);
                // else printf("Shot %04d: the number of traces is %d.\n", ishot + 1, tr.ntr);
                /* */
                record2d[ishot].shotID = ishot + 1;
                record2d[ishot].ntr = tr.ntr;
                record2d[ishot].sx = scale_header(tr.sx, tr.scalco);
                record2d[ishot].sz = scale_header(tr.sdepth, tr.scalel);
                record2d[ishot].gx = alloc1float(record2d[ishot].ntr);
                record2d[ishot].gz = alloc1float(record2d[ishot].ntr);
                record2d[ishot].ns = tr.ns;
                record2d[ishot].dt = (float)tr.dt / 1e6f;
                record2d[ishot].sampleRate = (nt - 1) / (record2d[ishot].ns - 1);
                record2d[ishot].traces_p = alloc2float(record2d[ishot].ns, record2d[ishot].ntr);
                dtr = 0;
                for (itrace = 0; itrace < record2d[ishot].ntr; itrace++)
                {
                    if (itrace > 0)
                    {
                        fvgettr(fp, &tr);
                    }
                    /* */
                    record2d[ishot].gx[itrace] = scale_header(tr.gx, tr.scalco);
                    record2d[ishot].gz[itrace] = scale_header(tr.gelev, tr.scalel);
                    memcpy(record2d[ishot].traces_p[itrace], tr.data, tr.ns * FSIZE);
                    dtr += (double)tr.d2;
                }
                record2d[ishot].dtr = dtr / record2d[ishot].ntr;
            }
        }
        fclose(fp);
    }
    else
    {
        err("<error>: didn't specify recordfile");
    }
}

Record2D *input_record2d(const char *recordfile, const char *srcfile,
                         int nt, float dt, int Nshot)
{
    /* 
    input_record2d: load 2D shot records and attach source wavelet

    Input:
        recordfile - path to shot gather file (const char *)
        srcfile    - path to binary file containing 1D source wavelet (const char *)
        nt         - number of time samples in wavelet and gathers (int)
        dt         - temporal sampling interval [s] (float)
        Nshot      - number of shots to read (int)

    Output:
        return     - pointer to allocated array of Record2D structures, each including:
                    .sx, .sz       : source positions
                    .rx, .rz       : receiver positions
                    .src           : source wavelet (nt x 1 float array)
                    .fpeak         : estimated peak frequency [Hz] 
    */
    
    int ishot, it, Nshot_act;
    float *src;

    src = alloc1float(nt);

    Nshot_act = get_shots_number(recordfile);
    if (Nshot != Nshot_act)
        err("<error>: The number of input shots (%d) does not match the number of shots in the shot gather file (%d)!", Nshot, Nshot_act);
    Record2D *record2d = alloc_record2d(Nshot);
    // import shot information
    read_gather2d_from_sufile(record2d, recordfile, Nshot, nt);

    // attach source wavelet and compute fpeak for each shot
    read_bin1d(srcfile, nt, src);
    for (ishot = 0; ishot < Nshot; ishot++)
    {
        // src
        record2d[ishot].src = alloc1float(nt);
        for (it = 0; it < nt; it++)
            record2d[ishot].src[it] = src[it];
        record2d[ishot].fpeak = calculate_peak_frequency(src, nt, dt);
    }

    free1float(src);
    return record2d;
}

char *read_string_file(const char *filename)
{
    /* 
    read_string_file: load entire file into a null-terminated string buffer

    Input:
        filename  - path to input file (const char *)

    Output:
        return    - pointer to allocated char buffer containing file content
              (null-terminated, must be freed by caller) 
    */
    
    FILE *fp;
    long len;
    char *data;

    fp = fopen(filename, "rb");
    if (!fp)
    {
        err("<error>: failed to open JSON file (%s).", filename);
        return NULL;
    }
    fseek(fp, 0, SEEK_END);
    len = ftell(fp);
    rewind(fp);
    data = (char *)malloc(len + 1);
    fread(data, 1, len, fp);
    data[len] = '\0';
    fclose(fp);
    return data;
}