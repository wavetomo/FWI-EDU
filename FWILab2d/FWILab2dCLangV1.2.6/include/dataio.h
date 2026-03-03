/*
 * @Author: QuYang Chen
 * @Date: 2025-09-25 17:18:33
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2025-12-23 16:53:55
 * @FilePath: /FWILab2dCLangV1.2.5/include/dataio.h
 * @Description: 
 * 
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved. 
 */
#pragma once

#include <su.h>
#include <segy.h>

#include "common.h"
#include "observation.h"

int get_precision_bytes(const char *precision);

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

long long get_file_size(const char *fileName);

/* 
get_file_size: get file size in bytes

Input:
  fileName  - path to the file (const char *fileName)

Output:
  return    - file size in bytes (long long), or -1LL if fseek/ftell fails 
*/

long long count_file_elements(const char *fileName, int elementSize);

/* 
count_file_elements: calculate number of elements in a binary file

Input:
  fileName     - path to the file (const char*)
  elementSize  - size of each element in bytes (int, must be > 0)

Output:
  return       - number of elements (long long), or -1LL if an error occurs
*/

void read_bin1d(const char *fileName, int n, float *data);

/* 
read_bin1d: read 1D float binary data from file

Input:
  fileName  - path to binary file (const char*)
  n         - number of float elements to read (int)
              if n <= 0, it will be inferred from file size

Output:
  data      - pointer to allocated float array 
*/

void read_bin2d(const char *fileName, int n1, int n2, float **data);

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

void read_bin3d(const char *fileName, int n1, int n2, int n3, float ***data);

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


void write_bin(const char *fileName, float *data, int n);

/* 
write_bin: write 1D float array to binary file

Input:
  fileName - path to output file (const char*)
  data     - pointer to float array
  n        - number of elements

Output:
    None (writes binary data to disk as files)
*/

void write_acoustic2d_traces(char *prefix, Record2D *record2d, int Nshot);

/* 
write_acoustic2d_traces: write 2D shot gathers to disk

Input:
  prefix      - file name prefix for output files (const char*)
  record2d    - array of shot gather structures containing traces and metadata (Record2D*)
  Nshot       - number of shots in the record2d array (int)

Output:
  None (writes traces to disk as files)
*/

float scale_header(int keyword, int scale);

/* 
scale_header: Apply scaling to a seismic trace header keyword

Input:
  keyword - original header value (int)
  scale   - scaling factor applied to the header (int)

Output:
  return  - scaled header value (float) 
*/

int get_shots_number(const char *filename);

/* 
get_shots_number: count the number of shots in a record file

Input:
  filename - path to the record file (const char *)

Output:
  return   - number of shots contained in the file (int) 
*/

void read_gather2d_from_sufile(Record2D *record2d, const char *recordfile,
    int Nshot, int nt);

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

Record2D* input_record2d(const char *recordfile, const char *srcfile,
                         int nt, float dt, int Nshot);

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

char *read_string_file(const char *filename);

/* 
read_string_file: load entire file into a null-terminated string buffer

Input:
  filename  - path to input file (const char *)

Output:
  return    - pointer to allocated char buffer containing file content
              (null-terminated, must be freed by caller)
*/