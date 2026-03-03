function [record2d, wavefields, grad] = iso_acoustic2d_propagation_engine(record2d, wavefields, vpm, rhom, ...
    dx, dz, dt, nt, pmlThick, freesurface, precondition, mode, verbose)
%% ISO_ACOUSTIC2D_PROPAGATION_ENGINE - Acoustic wave equation modeling engine
%
%   This function performs isotropic acoustic wavefield modeling and
%   computes forward or backward propagation for full waveform inversion.
%
%   [record2d, wavefields, grad] = iso_acoustic2d_propagation_engine(...)
%
%   INPUTS:
%       record2d     - structure array of shot gathers. Each element includes:
%                          .sx, .sz       : source location indices
%                          .rx, .rz       : receiver location indices
%                          .traces_p      : (nt x nrec) recorded data (only used in mode 2)
%                          .source_p      : (nt x 1) source time function
%
%       wavefields   - structure with fields used for wavefield storage:
%                          .u             : (nz x nx x nt) forward wavefield (used in mode 1/2)
%                          .lambda        : (nz x nx x nt) adjoint wavefield (used in mode 2)
%
%       vpm          - (nz x nx) P-wave velocity model [m/s]
%       rhom         - (nz x nx) density model [kg/m³]
%       dx, dz       - grid spacing in x and z directions [m]
%       dt           - temporal sampling interval [s]
%       nt           - total number of time steps
%       pmlThick     - PML thickness in number of grid points
%       freesurface  - flag to enable free surface (1 = yes, 0 = no)
%       precondition - apply pseudo-Hessian preconditioning (1 = yes, 0 = no)
%       mode         - operation mode:
%                          0 = forward modeling only (record only)
%                          1 = forward modeling (record + wavefield + pseudo-Hessian preconditioning)
%                          2 = adjoint modeling (gradient computation)
%       verbose      - display level:
%                          0 = silent
%                          1 = essential info
%                          2 = full visualization/debug info
%
%   OUTPUTS:
%       record2d     - updated with synthetic traces (.traces_p)
%       wavefields   - updated with wavefields (u, lambda depending on mode)
%       grad         - (nz x nx) gradient of misfit with respect to vpm (if mode == 2),
%                      zero matrix otherwise
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

[nz, nx] = size(vpm);

invdx = 1.0 / dx;
invdz = 1.0 / dz;
invdxdx = 1.0 / dx / dx;
invdzdz = 1.0 / dz / dz;
invdtdt = 1.0 / dt / dt;

Nx = nx + 2 * pmlThick;
Nz = nz + 2 * pmlThick;

isx = record2d.isx + pmlThick;
isz = record2d.isz + pmlThick;
ig = record2d.igz + pmlThick + (record2d.igx + pmlThick - 1) * Nz;

x = ((record2d.ix_min - 1) + (0:nx - 1)) * dx / 1000;
z = (0:nz - 1) * dz / 1000;

%% MEAL
vp = pad_model2d(vpm, pmlThick);
rho = pad_model2d(rhom, pmlThick);

[D, Beta] = mealCoefs2d(vp, [], record2d.fpeak, 2.5, dx, dz, pmlThick);
D = D ./ Beta;
vp = vp ./ Beta;
rho = rho .* Beta;

coef1 = (2.0 - dt * dt .* (D.^2)) ./ (1.0 + dt .* D);
coef2 = -1.0 * (1.0 - dt .* D) ./ (1.0 + dt .* D);
coef3 = dt * dt ./ (1.0 + dt .* D);

%% grid parameters
[K, invrho, invrhox, invrhoz] = acoustic_iso2d_rgfd_parameter(vp, rho, dx, dz);

%% alloc memory
if mode == 0 || mode == 1
    record2d.traces_p = zeros(record2d.ns, record2d.ntr);
end
if mode == 0
    wavefields = 0;
    grad = 0;
elseif mode == 1
    mst = 5;
    Nt = floor((nt - 1)/mst) + 1;
    wavefields = zeros(Nz, Nx, Nt);
    if precondition == 1
        grad = zeros(nz, nx);
    else
        grad = 0;
    end
elseif mode == 2
    mst = 5;
    grad = zeros(nz, nx);
end

%% initial wavefields
p0 = zeros(Nz, Nx);
p1 = zeros(Nz, Nx);
p2 = zeros(Nz, Nx);

%%
% kx = generate_frequency(dx, Nx)';
% kz = generate_frequency(dz, Nz);
N = 3;
C1 = 0.77793115; % optimized cofficients for 1st order derivative based on L0 norm
C2 = -0.17388691;
C3 = 0.02338713;
D0 = -2.8215452; % optimized cofficients for 2nd order derivative based on L0 norm
D1 = 1.57663238;
D2 = -0.18347238;
D3 = 0.01761260;
% C1 = 0.91067892; % omptimized cofficients for 1st order derivative based on L0 norm
% C2 = -0.34187892;
% C3 = 0.13833962;
% C4 = -0.04880710;
% C5 = 0.01302148;
% C6 = -0.00199047;
% D0 = -3.12513824; % optimized cofficients for 2nd order derivative based on L0 norm
% D1 = 1.84108651;
% D2 = -0.35706478;
% D3 = 0.10185626;
% D4 = -0.02924772;
% D5 = 0.00696837;
% D6 = -0.00102952;

ixStart = N + 1;
ixEnd = Nx - N;
if freesurface == 0
    izStart = N + 1;
else
    izStart = pmlThick + 1;
end
izEnd = Nz - N;
if mode == 0 || mode == 1
    if verbose > 1
        hFig = figure('Name', 'Forward Wavefield');
    end
    for it = 1:nt
        % finite difference
        for ix = ixStart:ixEnd
            for iz = izStart:izEnd
                dpdz = (C1 * (p1(iz+1, ix) - p1(iz-1, ix)) + C2 * (p1(iz+2, ix) - p1(iz-2, ix)) + ...
                    C3 * (p1(iz+3, ix) - p1(iz-3, ix))) * invdz;
                dpdx = (C1 * (p1(iz, ix+1) - p1(iz, ix-1)) + C2 * (p1(iz, ix+2) - p1(iz, ix-2)) + ...
                    C3 * (p1(iz, ix+3) - p1(iz, ix-3))) * invdx;
                dpdzdz = (D0 * p1(iz, ix) + D1 * (p1(iz+1, ix) + p1(iz-1, ix)) + ...
                    D2 * (p1(iz+2, ix) + p1(iz-2, ix)) + D3 * (p1(iz+3, ix) + p1(iz-3, ix))) * invdzdz;
                dpdxdx = (D0 * p1(iz, ix) + D1 * (p1(iz, ix+1) + p1(iz, ix-1)) + ...
                    D2 * (p1(iz, ix+2) + p1(iz, ix-2)) + D3 * (p1(iz, ix+3) + p1(iz, ix-3))) * invdxdx;
                % dpdx = (C1 * (p1(iz, ix+1) - p1(iz, ix-1)) + C2 * (p1(iz, ix+2) - p1(iz, ix-2)) + ...
                %     C3 * (p1(iz, ix+3) - p1(iz, ix-3)) + C4 * (p1(iz, ix+4) - p1(iz, ix-4)) + ...
                %     C5 * (p1(iz, ix+5) - p1(iz, ix-5)) + C6 * (p1(iz, ix+6) - p1(iz, ix-6))) * invdx;
                % dpdz = (C1 * (p1(iz+1, ix) - p1(iz-1, ix)) + C2 * (p1(iz+2, ix) - p1(iz-2, ix)) + ...
                %     C3 * (p1(iz+3, ix) - p1(iz-3, ix)) + C4 * (p1(iz+4, ix) - p1(iz-4, ix)) + ...
                %     C5 * (p1(iz+5, ix) - p1(iz-5, ix)) + C6 * (p1(iz+6, ix) - p1(iz-6, ix))) * invdz;
                % dpdxdx = (D0 * p1(iz, ix) + D1 * (p1(iz, ix+1) + p1(iz, ix-1)) + ...
                %     D2 * (p1(iz, ix+2) + p1(iz, ix-2)) + D3 * (p1(iz, ix+3) + p1(iz, ix-3)) + ...
                %     D4 * (p1(iz, ix+4) + p1(iz, ix-4)) + D5 * (p1(iz, ix+5) + p1(iz, ix-5)) + ...
                %     D6 * (p1(iz, ix+6) + p1(iz, ix-6))) * invdxdx;
                % dpdzdz = (D0 * p1(iz, ix) + D1 * (p1(iz+1, ix) + p1(iz-1, ix)) + ...
                %     D2 * (p1(iz+2, ix) + p1(iz-2, ix)) + D3 * (p1(iz+3, ix) + p1(iz-3, ix)) + ...
                %     D4 * (p1(iz+4, ix) + p1(iz-4, ix)) + D5 * (p1(iz+5, ix) + p1(iz-5, ix)) + ...
                %     D6 * (p1(iz+6, ix) + p1(iz-6, ix))) * invdzdz;
                p2(iz, ix) = coef1(iz, ix) * p1(iz, ix) + coef2(iz, ix) * p0(iz, ix) + coef3(iz, ix) * ...
                    K(iz, ix) * (invrhox(iz, ix) * dpdx + invrhoz(iz, ix) * dpdz + ...
                    invrho(iz, ix) * (dpdxdx + dpdzdz));
            end
        end
        % spec = fft(p1, Nz, 1);
        % dpdz = real(ifft(1i*kz.*spec, Nz, 1));
        % dpdzdz = real(ifft(-kz.*kz.*spec, Nz, 1));
        % spec = fft(p1, Nx, 2);
        % dpdx = real(ifft(1i*kx.*spec, Nx, 2));
        % dpdxdx = real(ifft(-kx.*kx.*spec, Nx, 2));
        % p2 = coef1 .* p1 + coef2 .* p0 + coef3 .* ...
        %     K .* (invrhox .* dpdx + invrhoz .* dpdz + invrho .* (dpdxdx + dpdzdz));
        % load source
        p2(isz, isx) = p2(isz, isx) + record2d.src(it) * K(isz, isx) * dt * dt * invdxdx * invdzdz;
        % free-surface boundary condition
        if freesurface == 1
            p2(pmlThick+1, :) = 0;
            for iz = 1:pmlThick
                p2(pmlThick+1-iz, :) = -p2(pmlThick+1+iz, :);
            end
        end
        % receiver
        if mod(it-1, record2d.sampleRate) == 0
            k = (it - 1) / record2d.sampleRate + 1;
            record2d.traces_p(k, :) = p1(ig);
        end
        % save wavefield
        if mode == 1 && mod(it-1, mst) == 0
            k = (it - 1) / mst + 1;
            wavefields(:, :, k) = p1;
            if precondition == 1
                dpdt2 = (p2 - 2 .* p1 + p0) * invdtdt;
                grad = grad + dpdt2(pmlThick+1:pmlThick+nz, pmlThick+1:pmlThick+nx).^2;
            end
        end
        % output wavfield
        if verbose > 1 && mod(it-1, 1000) == 0 && it > 1
            snap = p1(pmlThick+1:pmlThick+nz, pmlThick+1:pmlThick+nx);
            hFig = plot_tiled_images(1, 1, hFig, x/1000, z/1000, 1, 0.5, snap);
            set_all_axes_labels(hFig, 'Distance / km', 'Depth / km');
            set_subplot_titles('Forward Wavefield');
            colormap(gca, "gray"), drawnow;
        end
        % swap wavefield
        p0 = p1;
        p1 = p2;
    end
    if mode == 1 && precondition == 1
        grad = 4.0 .* grad ./ (rhom .* vpm).^2; % precondition
    end
elseif mode == 2
    if verbose > 1
        hFig = figure('Name', 'Backward Wavefield');
    end
    src = interpolate_traces(record2d.residual_p, record2d.sampleRate);
    for it = nt:-1:1
        % load source
        p1(ig) = p1(ig) - src(it, :)'.* K(ig) * dt * dt * invdxdx * invdzdz;
        % finite difference
        for ix = ixStart:ixEnd
            for iz = izStart:izEnd
                dpdz = (C1 * (p1(iz+1, ix) - p1(iz-1, ix)) + C2 * (p1(iz+2, ix) - p1(iz-2, ix)) + ...
                    C3 * (p1(iz+3, ix) - p1(iz-3, ix))) * invdz;
                dpdx = (C1 * (p1(iz, ix+1) - p1(iz, ix-1)) + C2 * (p1(iz, ix+2) - p1(iz, ix-2)) + ...
                    C3 * (p1(iz, ix+3) - p1(iz, ix-3))) * invdx;
                dpdzdz = (D0 * p1(iz, ix) + D1 * (p1(iz+1, ix) + p1(iz-1, ix)) + ...
                    D2 * (p1(iz+2, ix) + p1(iz-2, ix)) + D3 * (p1(iz+3, ix) + p1(iz-3, ix))) * invdzdz;
                dpdxdx = (D0 * p1(iz, ix) + D1 * (p1(iz, ix+1) + p1(iz, ix-1)) + ...
                    D2 * (p1(iz, ix+2) + p1(iz, ix-2)) + D3 * (p1(iz, ix+3) + p1(iz, ix-3))) * invdxdx;
                p2(iz, ix) = coef1(iz, ix) * p1(iz, ix) + coef2(iz, ix) * p0(iz, ix) + coef3(iz, ix) * ...
                    K(iz, ix) * (invrhox(iz, ix) * dpdx + invrhoz(iz, ix) * dpdz + ...
                    invrho(iz, ix) * (dpdxdx + dpdzdz));
            end
        end
        % free-surface boundary condition
        if freesurface == 1
            p2(pmlThick+1, :) = 0;
            for iz = 1:pmlThick
                p2(pmlThick+1-iz, :) = -p2(pmlThick+1+iz, :);
            end
        end
        % gradient
        if mod(it-1, mst) == 0
            k = (it - 1) / mst + 1;
            dpdt2 = (p2 - 2 .* p1 + p0) * invdtdt;
            grad = grad + dpdt2(pmlThick+1:pmlThick+nz, pmlThick+1:pmlThick+nx) .* ...
                wavefields(pmlThick+1:pmlThick+nz, pmlThick+1:pmlThick+nx, k);
        end
        % output wavfield
        if verbose > 1 && mod(it-1, 1000) == 0 && it < nt
            snap = p1(pmlThick+1:pmlThick+nz, pmlThick+1:pmlThick+nx);
            hFig = plot_tiled_images(1, 1, hFig, x/1000, z/1000, 1, 0.5, snap);
            set_all_axes_labels(hFig, 'Distance / km', 'Depth / km');
            set_subplot_titles('Backward Wavefield');
            colormap(gca, "gray"), drawnow;
        end
        % swap wavefield
        p0 = p1;
        p1 = p2;
    end
    grad = 2.0 .* grad ./ (rhom .* vpm); % gradient
    wavefields = 0;
end
if verbose > 1
    close(hFig);
end

end