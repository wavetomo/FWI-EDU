/*
 * @Author: QuYang Chen
 * @Date: 2025-09-25 17:22:42
 * @LastEditors: QuYang Chen
 * @LastEditTime: 2026-01-23 20:27:42
 * @FilePath: /FWILab2dCLangV1.2.5/src/propagation.c
 * @Description:
 *
 * Copyright (c) 2025 by WaveTomo, All Rights Reserved.
 */
#include "propagation.h"

void acoustic_iso2d_fd_it_6order(float **restrict p2, float **restrict p1, float **restrict p0,
                                float **restrict invrho, float **restrict invrhox, float **restrict invrhoz, float **restrict K,
                                float **restrict coef1, float **restrict coef2, float **restrict coef3,
                                float invdx, float invdz, float invdxdx, float invdzdz,
                                int ixStart, int ixEnd, int izStart, int izEnd)
{
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
// finite-difference coefficients
const float C1 = 0.77793115f;
const float C2 = -0.17388691f;
const float C3 = 0.02338713f;
const float D0 = -2.8215452f;
const float D1 = 1.57663238f;
const float D2 = -0.18347238f;
const float D3 = 0.01761260f;

int ix, iz;
float dpdz, dpdx, dpdzdz, dpdxdx;

#pragma omp parallel for private(ix, iz, dpdz, dpdx, dpdzdz, dpdxdx) \
    firstprivate(ixStart, ixEnd, izStart, izEnd, invdz, invdx, invdzdz, invdxdx, C1, C2, C3, D0, D1, D2, D3)
for (ix = ixStart; ix <= ixEnd; ++ix)
{
#pragma ivdep
    for (iz = izStart; iz <= izEnd; ++iz)
    {
    // finite difference
    dpdz = (C1 * (p1[ix][iz + 1] - p1[ix][iz - 1]) + C2 * (p1[ix][iz + 2] - p1[ix][iz - 2]) + C3 * (p1[ix][iz + 3] - p1[ix][iz - 3])) * invdz;
    dpdx = (C1 * (p1[ix + 1][iz] - p1[ix - 1][iz]) + C2 * (p1[ix + 2][iz] - p1[ix - 2][iz]) + C3 * (p1[ix + 3][iz] - p1[ix - 3][iz])) * invdx;
    dpdzdz = (D0 * p1[ix][iz] + D1 * (p1[ix][iz + 1] + p1[ix][iz - 1]) + D2 * (p1[ix][iz + 2] + p1[ix][iz - 2]) + D3 * (p1[ix][iz + 3] + p1[ix][iz - 3])) * invdzdz;
    dpdxdx = (D0 * p1[ix][iz] + D1 * (p1[ix + 1][iz] + p1[ix - 1][iz]) + D2 * (p1[ix + 2][iz] + p1[ix - 2][iz]) + D3 * (p1[ix + 3][iz] + p1[ix - 3][iz])) * invdxdx;

    p2[ix][iz] = coef1[ix][iz] * p1[ix][iz] + coef2[ix][iz] * p0[ix][iz] + coef3[ix][iz] * K[ix][iz] * (invrhox[ix][iz] * dpdx + invrhoz[ix][iz] * dpdz + invrho[ix][iz] * (dpdxdx + dpdzdz));
    }
}
}

void acoustic_iso2d_freesurface(float **p2, int Nx, int Nz, int pmlThick)
{
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

int ix, iz;
#pragma omp parallel for private(ix, iz) \
    firstprivate(Nx, Nz, pmlThick)
for (ix = 0; ix < Nx; ix++)
{
    p2[ix][pmlThick] = 0.0f;
#pragma ivdep
    for (iz = 0; iz <= pmlThick; iz++)
    p2[ix][pmlThick - iz] = -p2[ix][pmlThick + iz];
}
}

void acoustic_iso2d_rgfd_parameter(float **vp, float **rho, int nx, int nz,
                                float dx, float dz,
                                float **K, float **invrho,
                                float **invrhox, float **invrhoz)
{
    /*
    acoustic_iso2d_rgfd_parameter: compute acoustic model parameters for 2D grid

    Input:
        vp[nx][nz]    - velocity model
        rho[nx][nz]   - density model
        nx, nz        - model size
        dx, dz        - grid spacing

    Output:
        K[nx][nz]         - bulk modulus
        invrho[nx][nz]    - reciprocal of density
        invrhox[nx][nz]   - x-derivative of 1/rho (optimized stencil)
        invrhoz[nx][nz]   - z-derivative of 1/rho (optimized stencil)
    */

    // coefficients for optimized 1st-order derivative (L0 norm)
    const float C1 = 0.91067892f;
    const float C2 = -0.34187892f;
    const float C3 = 0.13833962f;
    const float C4 = -0.04880710f;
    const float C5 = 0.01302148f;
    const float C6 = -0.00199047f;

    float invdx = 1.0f / dx;
    float invdz = 1.0f / dz;

    // compute K and invrho
    for (int ix = 0; ix < nx; ix++)
    {
#pragma ivdep
        for (int iz = 0; iz < nz; iz++)
        {
            K[ix][iz] = rho[ix][iz] * vp[ix][iz] * vp[ix][iz];
            invrho[ix][iz] = 1.0f / rho[ix][iz];
        }
    }

    calculate_gradient2d(invrho, nx, nz, dx, dz, invrhox, invrhoz);

    // compute invrhox and invrhoz using optimized 6th-order central differences
    for (int ix = 6; ix < nx - 6; ix++)
    {
#pragma ivdep
        for (int iz = 6; iz < nz - 6; iz++)
        {
            invrhox[ix][iz] = (C1 * (invrho[ix + 1][iz] - invrho[ix - 1][iz]) +
                            C2 * (invrho[ix + 2][iz] - invrho[ix - 2][iz]) +
                            C3 * (invrho[ix + 3][iz] - invrho[ix - 3][iz]) +
                            C4 * (invrho[ix + 4][iz] - invrho[ix - 4][iz]) +
                            C5 * (invrho[ix + 5][iz] - invrho[ix - 5][iz]) +
                            C6 * (invrho[ix + 6][iz] - invrho[ix - 6][iz])) *
                            invdx;

            invrhoz[ix][iz] = (C1 * (invrho[ix][iz + 1] - invrho[ix][iz - 1]) +
                            C2 * (invrho[ix][iz + 2] - invrho[ix][iz - 2]) +
                            C3 * (invrho[ix][iz + 3] - invrho[ix][iz - 3]) +
                            C4 * (invrho[ix][iz + 4] - invrho[ix][iz - 4]) +
                            C5 * (invrho[ix][iz + 5] - invrho[ix][iz - 5]) +
                            C6 * (invrho[ix][iz + 6] - invrho[ix][iz - 6])) *
                            invdz;
        }
    }
}

static float absorbWeight(float alpha, int l, int nlayer)
{
    /*
    absorbWeight: compute the damping weight for an absorbing boundary layer

    Inputs:
        alpha   - absorption strength exponent (controls damping profile shape)
        l       - current layer index (0 ≤ l < nlayer)
        nlayer  - total number of absorbing layers

    Output:
        return  - damping weight (float), increases smoothly from 0 to 1 across boundary
    */
    float w, d;
    d = (0.1f * l) / (0.1f * nlayer);
    w = expf(logf(2.0f) * powf(d, alpha)) - 1.0f;
    return w;
}

void mealCoefs2d(float **D, float **Beta, float **vp, float **vs, float fpeak, float eal_alpha,
                float dz, float dx, int NZ, int NX, int nlayer)
{
    /*
    mealCoefs2d: Construct 2D MEAL (Modified Efficient Absorbing Layer) boundary condition

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

    int ix, iz, l;
    float weight;
    float pd, Rcoef, PPWfd, PPWfcs;
    float thicknessx, thicknessz;
    float d0x, d0z;
    float betax, betaz;
    float **Vs = NULL;
    float **Dx = NULL, **Dz = NULL;
    float **Betax = NULL, **Betaz = NULL;
    pd = 2;
    Rcoef = powf(10.0f, -1.0f * ((log10f(nlayer) - 1.0f) / (log10f(2.0f))) - 3.0f - 0.5f * pd);
    PPWfd = 10.0f;
    thicknessx = nlayer * dx;
    thicknessz = nlayer * dz;
    Vs = alloc2float(NZ, NX);
    Dx = alloc2float(NZ, NX);
    Dz = alloc2float(NZ, NX);
    Betax = alloc2float(NZ, NX);
    Betaz = alloc2float(NZ, NX);
    if (eal_alpha == -1.0f)
    {
        eal_alpha = 2.5f;
    }
    if (vs == NULL)
    {
#pragma omp parallel for private(ix, iz)
        for (ix = 0; ix < NX; ix++)
        {
#pragma ivdep
            for (iz = 0; iz < NZ; iz++)
            {
                Vs[ix][iz] = vp[ix][iz] / sqrtf(3.0f);
            }
        }
    }
    else
    {
#pragma omp parallel for private(ix, iz)
        for (ix = 0; ix < NX; ix++)
        {
#pragma ivdep
            for (iz = 0; iz < NZ; iz++)
            {
                if (vs[ix][iz] == 0.0f)
                {
                    Vs[ix][iz] = vp[ix][iz] / sqrtf(3.0f);
                }
                else
                {
                    Vs[ix][iz] = vs[ix][iz];
                }
            }
        }
    }
    /* initialization */
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < NX; ix++)
    {
#pragma ivdep
        for (iz = 0; iz < NZ; iz++)
        {
            D[ix][iz] = 0.0f;
            Dx[ix][iz] = 0.0f;
            Dz[ix][iz] = 0.0f;
            Beta[ix][iz] = 1.0f;
            Betax[ix][iz] = 1.0f;
            Betaz[ix][iz] = 1.0f;
        }
    }
    /* top and bottom boundary */
#pragma omp parallel for private(ix, iz, l, weight, d0z, PPWfcs, betaz) \
    firstprivate(PPWfd, Rcoef, thicknessz, dz)
    for (ix = 0; ix < NX; ix++)
    {
#pragma ivdep
        for (l = 1; l <= nlayer; l++)
        {
            // weight = powf((1.0f * l) / (1.0f * nlayer), pd);
            weight = absorbWeight(eal_alpha, l, nlayer);
            /* top */
            iz = nlayer - l;
            d0z = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
            PPWfcs = Vs[ix][iz] / (dz * fpeak);
            betaz = 2.0f * PPWfcs / PPWfd;
            Dz[ix][iz] = d0z * weight;
            Betaz[ix][iz] = 1.0f + (betaz - 1.0f) * weight;
            /* bottom */
            iz = NZ - nlayer - 1 + l;
            d0z = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessz);
            PPWfcs = Vs[ix][iz] / (dz * fpeak);
            betaz = 2.0f * PPWfcs / PPWfd;
            Dz[ix][iz] = d0z * weight;
            Betaz[ix][iz] = 1.0f + (betaz - 1.0f) * weight;
        }
    }
    /* left and right boundary */
#pragma omp parallel for private(ix, iz, l, weight, d0x, PPWfcs, betax) \
    firstprivate(PPWfd, Rcoef, thicknessx, dx)
    for (l = 1; l <= nlayer; l++)
    {
        // weight = powf((1.0f * l) / (1.0f * nlayer), pd);
        weight = absorbWeight(eal_alpha, l, nlayer);
#pragma ivdep
        for (iz = 0; iz < NZ; iz++)
        {
            /* left */
            ix = nlayer - l;
            d0x = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
            PPWfcs = Vs[ix][iz] / (dx * fpeak);
            betax = 2.0f * PPWfcs / PPWfd;
            Dx[ix][iz] = d0x * weight;
            Betax[ix][iz] = 1.0f + (betax - 1.0f) * weight;
            /* right */
            ix = NX - nlayer - 1 + l;
            d0x = -1.0f * (1.0f + pd) * vp[ix][iz] * logf(Rcoef) / (2.0f * thicknessx);
            PPWfcs = Vs[ix][iz] / (dx * fpeak);
            betax = 2.0f * PPWfcs / PPWfd;
            Dx[ix][iz] = d0x * weight;
            Betax[ix][iz] = 1.0f + (betax - 1.0f) * weight;
        }
    }
    /* absorbing coefficent */
#pragma omp parallel for private(ix, iz)
    for (ix = 0; ix < NX; ix++)
    {
#pragma ivdep
        for (iz = 0; iz < NZ; iz++)
        {
            D[ix][iz] = sqrtf(powf(Dx[ix][iz], 2.0f) + powf(Dz[ix][iz], 2.0f));
            Beta[ix][iz] =
                sqrtf(powf(Betax[ix][iz] - 1.0f, 2.0f) + powf(Betaz[ix][iz] - 1.0f, 2.0f)) + 1.0f;
        }
    }
    /* free memory */
    free2float(Vs);
    free2float(Dx);
    free2float(Dz);
    free2float(Betax);
    free2float(Betaz);
}

void iso_acoustic2d_propagation_engine(Record2D *record2d,
                                    Acoustic2dCheckPoint *chkpt, /* allocated if mode==1, else NULL */
                                    float **vpm, float **rhom, int nx, int nz,
                                    float dx, float dz, float dt, int nt,
                                    int pmlThick, int freesurface, int precondition,
                                    int mode, int iter, int verbose,
                                    float **grad /* allocated if mode==1 or 2, else NULL */)
{
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

    // finite-difference coefficients
    const int N = 3;

    int Nx, Nz, Np;

    int ix, iz, it, k, itr, ishot, igx, igz, isx, isz;
    int mst, Nt, ntr, sampleRate, ns;
    int ixStart, ixEnd, izStart, izEnd;

    int ichkpt, itStart, itEnd, itt;

    float **vp, **rho;
    float **K, **invrho, **invrhox, **invrhoz;
    float **D, **Beta;
    float **coef1, **coef2, **coef3;
    float dpdz, dpdx, dpdzdz, dpdxdx;
    float invdx, invdz, invdxdx, invdzdz, invdtdt;

    float **p0, **p1, **p2, **p0_fwd, **p1_fwd, **p2_fwd, **ptmp;

    float **adjsrc;

    double **precond_cumul, **grad_cumul;

    char filename[256];

    int shotInterval, snapInterval;
    char *pStrTemp;

    if (verbose > 1)
    {
        // read parameters of shot output interval
        if (getenv("shotInterval") != NULL)
        {
            pStrTemp = getenv("shotInterval");
            shotInterval = atoi(pStrTemp);
        }
        else
            shotInterval = 20;
        // read parameters of snap output interval
        if (getenv("snapInterval") != NULL)
        {
            pStrTemp = getenv("snapInterval");
            snapInterval = atoi(pStrTemp);
        }
        else
            snapInterval = 1000;
    }

    invdx = 1.0f / dx;
    invdz = 1.0f / dz;
    invdxdx = 1.0f / (dx * dx);
    invdzdz = 1.0f / (dz * dz);
    invdtdt = 1.0f / (dt * dt);

    Nx = nx + 2 * pmlThick;
    Nz = nz + 2 * pmlThick;
    Np = Nx * Nz;

    ishot = record2d->shotID - 1;
    isx = record2d->isx + pmlThick;
    isz = record2d->isz + pmlThick;
    ntr = record2d->ntr;
    sampleRate = record2d->sampleRate;

    // alloc memory
    vp = alloc2float(Nz, Nx);
    rho = alloc2float(Nz, Nx);
    D = alloc2float(Nz, Nx);
    Beta = alloc2float(Nz, Nx);
    coef1 = alloc2float(Nz, Nx);
    coef2 = alloc2float(Nz, Nx);
    coef3 = alloc2float(Nz, Nx);
    K = alloc2float(Nz, Nx);
    invrho = alloc2float(Nz, Nx);
    invrhox = alloc2float(Nz, Nx);
    invrhoz = alloc2float(Nz, Nx);
    p0 = alloc2float(Nz, Nx);
    p1 = alloc2float(Nz, Nx);
    p2 = alloc2float(Nz, Nx);
    if (mode == 2)
    {
        adjsrc = alloc2float(nt, ntr);
        p0_fwd = alloc2float(Nz, Nx);
        p1_fwd = alloc2float(Nz, Nx);
        p2_fwd = alloc2float(Nz, Nx);
    }

    if (mode == 1)
    {
        mst = MST;
        
        if (precondition == 1) // save preconditon in grad[0]
        {
        precond_cumul = alloc2double(nz, nx);
        memset(precond_cumul[0], 0, nx * nz * DSIZE);
        memset(grad[0], 0, nx * nz * FSIZE);
        }
            
    }
    else if (mode == 2)
    {
        mst = MST;
        grad_cumul = alloc2double(nz, nx);
        memset(grad_cumul[0], 0, nx * nz * DSIZE);
        memset(grad[0], 0, nx * nz * FSIZE);
        memset(p0_fwd[0], 0, Nx * Nz * FSIZE);
        memset(p1_fwd[0], 0, Nx * Nz * FSIZE);
        memset(p2_fwd[0], 0, Nx * Nz * FSIZE);
    }
    memset(p0[0], 0, Nx * Nz * FSIZE);
    memset(p1[0], 0, Nx * Nz * FSIZE);
    memset(p2[0], 0, Nx * Nz * FSIZE);

    pad_model2d(vpm, nz, nx, vp, pmlThick);
    pad_model2d(rhom, nz, nx, rho, pmlThick);
    
    mealCoefs2d(D, Beta, vp, NULL, record2d->fpeak, 2.5f, dz, dx, Nz, Nx, pmlThick);

    float denom;
#pragma omp parallel for private(ix, iz, denom) \
    firstprivate(Nx, Nz, dt)
    for (ix = 0; ix < Nx; ix++)
    {
#pragma ivdep
        for (iz = 0; iz < Nz; iz++)
        {
            D[ix][iz] /= Beta[ix][iz];
            vp[ix][iz] /= Beta[ix][iz];
            rho[ix][iz] *= Beta[ix][iz];

            coef1[ix][iz] = (2.0 - dt * dt * D[ix][iz] * D[ix][iz]) / (1.0 + dt * D[ix][iz]);
            coef2[ix][iz] = -1.0 * (1.0 - dt * D[ix][iz]) / (1.0 + dt * D[ix][iz]);
            coef3[ix][iz] = (dt * dt) / (1.0 + dt * D[ix][iz]);
        }
    }

    // grid parameters
    acoustic_iso2d_rgfd_parameter(vp, rho, Nx, Nz, dx, dz, K, invrho, invrhox, invrhoz);

    ixStart = N;
    ixEnd = Nx - N - 1;
    izStart = (freesurface == 0) ? N : pmlThick;
    izEnd = Nz - N - 1;

    if (mode == 0 || mode == 1)
    {
        for (it = 0; it < nt; ++it)
        {

            // save wavefield with check-pointing method
            if (mode == 1)
            {
                save_acoustic2dCheckPoint(chkpt, p0, p1, Np, it);
            }

            // finite-difference time extrapolation
            acoustic_iso2d_fd_it_6order(p2, p1, p0,
                                        invrho, invrhox, invrhoz, K,
                                        coef1, coef2, coef3,
                                        invdx, invdz, invdxdx, invdzdz,
                                        ixStart, ixEnd, izStart, izEnd);

            // load source
            p2[isx][isz] += 1.0 * record2d->src[it] * K[isx][isz] * dt * dt * invdxdx * invdzdz;

            // free-surface boundary condition
            if (freesurface == 1)
            {
                acoustic_iso2d_freesurface(p2, Nx, Nz, pmlThick);
            }

            // receiver
            if (it % sampleRate == 0)
            {
                k = it / sampleRate;
                for (itr = 0; itr < ntr; ++itr)
                {
                    igx = record2d->igx[itr] + pmlThick;
                    igz = record2d->igz[itr] + pmlThick;
                    record2d->traces_p[itr][k] = p1[igx][igz];
                }
            }

            // save precondition
            if (mode == 1 && (it % mst) == 0)
            {
                if (precondition == 1)
                {
#pragma omp parallel for private(ix, iz) \
    firstprivate(Nx, Nz, pmlThick, invdtdt)
                    for (ix = pmlThick; ix < pmlThick + nx; ix++)
                    {
#pragma ivdep
                        for (iz = pmlThick; iz < pmlThick + nz; iz++)
                        {
                            precond_cumul[ix - pmlThick][iz - pmlThick] += (double)(SQR((p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) * invdtdt));
                        }
                    }
                }
            }

            // output wavfield
            if (verbose > 1 && ishot % shotInterval == 0 && it % snapInterval == 0)
            {
                sprintf(filename, "wavefield_forward_iter%03d_ishot%04d_it%05d_%d_%d_%.2fm_%.2fm.bin", iter, ishot, it, Nz, Nx, dz, dx);
                write_bin(filename, p1[0], Nz * Nx);
            }

            // swap wavefield
            ptmp = p0;
            p0 = p1;
            p1 = p2;
            p2 = ptmp;
        } /* end time loop */

        if (mode == 1 && precondition == 1)
        {
#pragma omp parallel for private(ix, iz) \
    firstprivate(nx, nz)
            for (ix = 0; ix < nx; ix++)
            {
#pragma ivdep
                for (iz = 0; iz < nz; iz++)
                {
                precond_cumul[ix][iz] = 4.0 * precond_cumul[ix][iz]  / ((double)(SQR(rhom[ix][iz] * vpm[ix][iz]))); // precondition
                grad[ix][iz] = (float)precond_cumul[ix][iz];
                }
            }
        }

        // if (ishot == Nshot / 2) {
        //     sprintf(filename, "shot%d_%d_%d.bin", ishot, record2d->ns, ntr);
        //     write_bin(filename, record2d->traces_p[0],  record2d->ns * record2d->ntr);
        // }

    } /* end forward mode */

    /* ---------- backward/adjoint (mode==2) ---------- */
    if (mode == 2)
    {
        interpolate_traces(record2d->traces_p, ntr, record2d->ns, sampleRate, adjsrc); /* adapt signature as needed */

        /* time-reversed loop */
        for (int it = nt - 1; it >= 0; --it)
        {
            /* calculate forward wavefield with check-point method*/
            itStart = get_acoustic2dCheckPointStatus(chkpt, it);
            if (itStart != -1)
            {
                itEnd = it;
                get_acoustic2dCheckPoint(chkpt, p0_fwd, p1_fwd, Np, itStart);
                for (itt = itStart; itt < itEnd; itt++)
                {
                    // finite-difference time extrapolation
                    acoustic_iso2d_fd_it_6order(p2_fwd, p1_fwd, p0_fwd,
                                                invrho, invrhox, invrhoz, K,
                                                coef1, coef2, coef3,
                                                invdx, invdz, invdxdx, invdzdz,
                                                ixStart, ixEnd, izStart, izEnd);

                    // free-surface boundary condition
                    if (freesurface == 1)
                    {
                        acoustic_iso2d_freesurface(p2_fwd, Nx, Nz, pmlThick);
                    }

                    // load source
                    p2_fwd[isx][isz] += record2d->src[it] * K[isx][isz] * dt * dt * invdxdx * invdzdz;

                    // save forward wavefield at it+1
                    if ((itt+1) % mst == 0)
                        save_acoustic2dWavefield(chkpt, p2_fwd, Np, itt + 1);

                    // swap wavefield
                    ptmp = p0_fwd;
                    p0_fwd = p1_fwd;
                    p1_fwd = p2_fwd;
                    p2_fwd = ptmp;
                }
            }

            /* calculate backward wavefield*/
            for (itr = 0; itr < ntr; ++itr)
            {
                // load source
                igx = record2d->igx[itr] + pmlThick;
                igz = record2d->igz[itr] + pmlThick;
                p1[igx][igz] -= adjsrc[itr][it] * K[igx][igz] * dt * dt * invdxdx * invdzdz;
            }

            // finite-difference time extrapolation
            acoustic_iso2d_fd_it_6order(p2, p1, p0,
                                        invrho, invrhox, invrhoz, K,
                                        coef1, coef2, coef3,
                                        invdx, invdz, invdxdx, invdzdz,
                                        ixStart, ixEnd, izStart, izEnd);

            // free-surface boundary condition
            if (freesurface == 1)
            {
                acoustic_iso2d_freesurface(p2, Nx, Nz, pmlThick);
            }



            // calculate gradient
            if (it % mst == 0)
            {
                // get forward wavefield
                get_acoustic2dWavefield(chkpt, p1_fwd, Np, it);                
#pragma omp parallel for private(ix, iz) \
    firstprivate(nx, nz, pmlThick, invdtdt)
                for (ix = pmlThick; ix < pmlThick + nx; ++ix)
                {
#pragma ivdep
                    for (iz = pmlThick; iz < pmlThick + nz; ++iz)
                    {
                    grad_cumul[ix - pmlThick][iz - pmlThick] += (double)((p2[ix][iz] - 2.0f * p1[ix][iz] + p0[ix][iz]) * invdtdt * p1_fwd[ix][iz]);
                    }
                }
            }

            // output wavfield
            if (verbose > 1 && ishot % shotInterval == 0 && it % snapInterval == 0)
            {
                sprintf(filename, "wavefield_forward_reconstruction_iter%03d_ishot%04d_it%05d_%d_%d_%.2fm_%.2fm.bin", iter, ishot, it, Nz, Nx, dz, dx);
                write_bin(filename, p1_fwd[0], Nz * Nx);
                sprintf(filename, "wavefield_backward_iter%03d_ishot%04d_it%05d_%d_%d_%.2fm_%.2fm.bin", iter, ishot, it, Nz, Nx, dz, dx);
                write_bin(filename, p1[0], Nz * Nx);
            }

            // swap wavefield
            ptmp = p0;
            p0 = p1;
            p1 = p2;
            p2 = ptmp;
        } /* end reverse time loop */

        // finalize gradient scaling
#pragma omp parallel for private(ix, iz) \
    firstprivate(nx, nz)
        for (ix = 0; ix < nx; ix++)
        {
#pragma ivdep
            for (iz = 0; iz < nz; iz++)
            {
            grad_cumul[ix][iz] = 2.0 * grad_cumul[ix][iz] / ((double)(rhom[ix][iz] * vpm[ix][iz]));
            grad[ix][iz] =  (float)grad_cumul[ix][iz];
            }
        }   
    } /* end mode==2 */

    free2float(vp);
    free2float(rho);
    free2float(D);
    free2float(Beta);
    free2float(coef1);
    free2float(coef2);
    free2float(coef3);
    free2float(K);
    free2float(invrho);
    free2float(invrhox);
    free2float(invrhoz);
    free2float(p0);
    free2float(p1);
    free2float(p2);
    if (mode == 1)
    {
    if(precondition == 1)
    {
        free2double(precond_cumul);
    }
    }
    else if (mode == 2)
    {
        free2float(adjsrc);
        free2double(grad_cumul);
        free2float(p0_fwd);
        free2float(p1_fwd);
        free2float(p2_fwd);
    }
}

void iso_acoustic2d_propagation(
    float **vp, float **rho, Record2D *record2d_obs,
    int Nshot, int nz, int nx,
    float dx, float dz, float dt, int nt,
    int pmlThick, int freesurface,
    float freq, float filter,
    int precondition, int problem, int iter,
    int verbose,
    Record2D *record2d_pre, float **grad)
{
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
        record2d_pre - predicted 2D shot records (Record2D*) if  problem==0,1 or residuals if problem==2
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

    int ix, iz, ix_min, ix_max, nxx, itr, ishot;
    int Nx, Nz, Nt;
    float **precond_per_shot, **grad_per_shot, **precond;
    double **grad_cumul, **precond_cumul;
    float vp_min, kpass;
    char filename[256];

    Acoustic2dCheckPoint chkpt;

    int shotInterval;
    char *pStrTemp;
    // read parameters of shot output interval
    if (getenv("shotInterval") != NULL)
    {
        pStrTemp = getenv("shotInterval");
        shotInterval = atoi(pStrTemp);
    }
    else
        shotInterval = 20;

    // Case 1: forward modeling only (problem == 0 or 2)
    if (problem == 0 || problem == 2)
    {
        for (int ishot = 0; ishot < Nshot; ishot++)
        {
            ix_min = record2d_obs[ishot].ix_min;
            ix_max = record2d_obs[ishot].ix_max;
            nxx = ix_max - ix_min + 1;

            // Forward modeling
            iso_acoustic2d_propagation_engine(&(record2d_pre[ishot]), NULL,
                                            vp + ix_min, rho + ix_min, nxx, nz, dx, dz, dt, nt,
                                            pmlThick, freesurface, precondition, 0, iter, verbose, NULL);

            // If problem == 2, compute adjoint source (residual)
            if (problem == 2)
            {
                compute_adjoint_source2d(&record2d_pre[ishot], &record2d_obs[ishot], vp, freq, filter, &record2d_pre[ishot]);
            }
        }
    }
    // Case 2: full waveform inversion gradient calculation (problem == 1)
    else if (problem == 1)
    {

        if (precondition == 1)
        {
            precond = alloc2float(nz, nx);
            precond_cumul = alloc2double(nz, nx);
            memset(precond_cumul[0], 0, nx * nz * DSIZE);
        }
        else
        {
            precond = NULL;
            precond_cumul = NULL;
        }

        grad_cumul = alloc2double(nz, nx);
        memset(grad_cumul[0], 0, nx * nz * DSIZE);

        // Loop over all shots
        for (ishot = 0; ishot < Nshot; ishot++)
        {
            ix_min = record2d_obs[ishot].ix_min;
            ix_max = record2d_obs[ishot].ix_max;
            nxx = ix_max - ix_min + 1;
            Nx = nxx + 2 * pmlThick;
            Nz = nz + 2 * pmlThick;

            acoustic2dCheckPoint_init(&chkpt);
            alloc_acoustic2dCheckPoint(&chkpt, Nz, Nx, nt);
            grad_per_shot = alloc2float(nz, nxx);

            if (precondition == 1)
                precond_per_shot = alloc2float(nz, nxx);
            else
                precond_per_shot = NULL;

            // Forward propagation
            iso_acoustic2d_propagation_engine(&(record2d_pre[ishot]), &chkpt,
                                            vp + ix_min, rho + ix_min, nxx, nz, dx, dz, dt, nt,
                                            pmlThick, freesurface, precondition, 1, iter, verbose, precond_per_shot);

            if (verbose > 0 && ishot % shotInterval == 0)
            {
                sprintf(filename, "predicted_iter%03d_shot%05d_%d_%d_%.2fms_%.2fm.bin", iter, ishot + 1, record2d_pre[ishot].ns,
                        record2d_pre[ishot].ntr, record2d_pre[ishot].dt * 1e3, record2d_pre[ishot].dtr);
                write_bin(filename, record2d_pre[ishot].traces_p[0], record2d_pre[ishot].ns * record2d_pre[ishot].ntr);
            }

            // Compute adjoint source (residual = record2d_pre)
            compute_adjoint_source2d(&record2d_pre[ishot], &record2d_obs[ishot], vp, freq, filter, &record2d_pre[ishot]);

            if (verbose > 0 && ishot % shotInterval == 0)
            {
                sprintf(filename, "adjsource_iter%03d_shot%05d_%d_%d_%.2fms_%.2fm.bin", iter, ishot + 1, record2d_pre[ishot].ns,
                        record2d_pre[ishot].ntr, record2d_pre[ishot].dt * 1e3, record2d_pre[ishot].dtr);
                write_bin(filename, record2d_pre[ishot].traces_p[0], record2d_pre[ishot].ns * record2d_pre[ishot].ntr);
            }

            // Backward propagation
            iso_acoustic2d_propagation_engine(&(record2d_pre[ishot]), &chkpt,
                                            vp + ix_min, rho + ix_min, nxx, nz, dx, dz, dt, nt,
                                            pmlThick, freesurface, precondition, 2, iter, verbose, grad_per_shot);

            if (verbose > 0 && ishot % shotInterval == 0)
            {
                sprintf(filename, "gradient_raw_iter%03d_shot%05d_%d_%d_%.2fm_%.2fm.bin", iter, ishot + 1, nz,
                        nxx, dz, dx);
                write_bin(filename, grad_per_shot[0], nz * nxx);
            }

            // Gradient preconditioning (per shot)
            if (precondition == 1)
            {
                apply_gradient_shot_precondition(grad_per_shot[0], precond_per_shot[0], nxx * nz);
            }

            // Accumulate shot gradient into global gradient and pseudo-Hessian matrix
#pragma omp parallel for private(ix, iz) \
            firstprivate(nz, ix_min, ix_max)
            for (ix = ix_min; ix <= ix_max; ix++)
            {
#pragma ivdep
                for (iz = 0; iz < nz; iz++)
                {
                    grad_cumul[ix][iz] += grad_per_shot[ix - ix_min][iz];
                }
            }

            if(precondition == 1)
            {
#pragma omp parallel for private(ix, iz) \
            firstprivate(nz, ix_min, ix_max)
                for (ix = ix_min; ix <= ix_max; ix++)
                {
#pragma ivdep
                    for (iz = 0; iz < nz; iz++)
                    {
                        precond_cumul[ix][iz] += precond_per_shot[ix - ix_min][iz];
                    }
                }
            }

            free_acoustic2dCheckPoint(&chkpt);
            free2float(grad_per_shot);
            if (precondition == 1)
            {
                free2float(precond_per_shot);
            }
        }

        for (ix = 0; ix < nz * nx; ix++)
            grad[0][ix] = (float)grad_cumul[0][ix];

        if(precondition == 1)
        {
            for (ix = 0; ix < nz * nx; ix++)
                precond[0][ix] = (float)precond_cumul[0][ix];
        }

        // Output
        if (verbose > 0)
        {
            sprintf(filename, "gradient_raw_iter%03d_%d_%d_%.2fm_%.2fm.bin", iter, nz, nx, dz, dx);
            write_bin(filename, grad[0], nz * nx);
            if (precondition == 1)
            {
                sprintf(filename, "precond_iter%03d_%d_%d_%.2fm_%.2fm.bin", iter, nz, nx, dz, dx);
                write_bin(filename, precond[0], nz * nx);
            }
        }
        
        // Apply preconditioning to global gradient
        if (precondition == 1)
        {
            apply_gradient_precondition(grad[0], precond[0], nz * nx, 0.005f);
        }

        // // Apply low-pass filtering to gradient
        // if (filter == 0)
        // {
        //     vp_min = find_min_float(vp[0], nz * nx);
        //     kpass = freq / vp_min;
        //     lowpass2d(grad, nz, nx, dz, dx, kpass, 2.0f * kpass, grad);
        // }
        // // Apply band-pass filtering to gradient
        // else if (filter == 1)
        // {
        //     vp_min = find_min_float(vp[0], nz * nx);
        //     kpass = freq / vp_min;
        //     bandpass2d(grad, nz, nx, dz, dx, kpass, 2.0f * kpass, grad);
        // }
        
        if (precondition == 1)
        {
            free2float(precond);
            free2double(precond_cumul);
        }

        free2double(grad_cumul);
    }
}
