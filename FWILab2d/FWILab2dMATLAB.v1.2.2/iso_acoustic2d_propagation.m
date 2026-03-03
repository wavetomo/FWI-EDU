function [record2d_pre, grad] = iso_acoustic2d_propagation(vp, rho, record2d_obs, ...
    dx, dz, dt, nt, pmlThick, freesurface, freq, filter, precondition, problem, isParallel, verbose)
%% ISO_ACOUSTIC2D_PROPAGATION 2D isotropic acoustic wave propagation and gradient calculation
%
%   [record2d_pre, grad] = iso_acoustic2d_propagation(vp, rho, record2d_obs, ...
%       dx, dz, dt, nt, pmlThick, freesurface, freq, filter, precondition, problem, parallel, verbose)
%
%   INPUTS:
%       vp           - (nz x nx) matrix of P-wave velocity model [m/s]
%       rho          - (nz x nx) matrix of density model [kg/m³]
%       record2d_obs - structure array containing observed seismic data.
%                      Each element corresponds to one shot gather and contains:
%                          - .sx, .sz: source location (in grid index)
%                          - .rx, .rz: receiver locations (arrays)
%                          - .traces_p: observed pressure wavefields [nt x nrec]
%
%       dx           - grid spacing in x-direction [m]
%       dz           - grid spacing in z-direction [m]
%       dt           - time sampling interval [s]
%       nt           - number of time steps
%       pmlThick     - thickness of absorbing boundary (PML) [grid points]
%       freesurface  - flag for free surface (1 = free surface, 0 = no)
%       freq         - reference frequency [Hz] for low-pass filtering (used in gradient)
%       filter       - filter parameters (used in backward propagation)
%       precondition - flag to apply pseudo-Hessian preconditioning (1 = yes, 0 = no)
%       problem      - problem type:
%                         0 = forward modeling for computing record only
%                         1 = compute FWI gradient
%                         2 = forward modeling for computing record and residual
%       parallel     - enable parallel computing using `parfor` (1 = yes, 0 = no)
%       verbose      - verbosity level:
%                         0 = no display
%                         1 = display main images without wavefields
%                         2 = display all images with wavefields
%
%   OUTPUTS:
%       record2d_pre - structure array of predicted data, same structure as record2d_obs,
%                      including .traces_p (modeled wavefields)
%       grad         - (nz x nx) FWI gradient (zero if problem ~= 1)
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

if isParallel && verbose > 1
    error('verbose must be 0 or 1, when parallel = 1!');
end
if isParallel && isempty(gcp('nocreate'))
    warning('Parallel mode requested, but no parallel pool found. Starting pool...');
    parpool;
end

[nz, nx] = size(vp);

x = (0:nx - 1) * dx / 1000;
z = (0:nz - 1) * dz / 1000;
t = (0:nt - 1) * dt;

Nshot = length(record2d_obs);
record2d_pre = copy_record2d_header(record2d_obs);
if problem == 0 || problem == 2
    if isParallel == 1
        parfor ishot = 1:Nshot
            ix_min = record2d_obs(ishot).ix_min;
            ix_max = record2d_obs(ishot).ix_max;
            vpm = vp(:, ix_min:ix_max);
            rhom = rho(:, ix_min:ix_max);
            [record2d_pre(ishot), ~, ~] = iso_acoustic2d_propagation_engine ...
                (record2d_pre(ishot), 0, vpm, rhom, dx, dz, dt, nt, ...
                pmlThick, freesurface, precondition, 0, verbose);
            if problem == 2
                record2d_pre(ishot) = compute_adjoint_source2d(record2d_pre(ishot), record2d_obs(ishot), vpm, freq, filter);
            end
        end
    else
        if verbose > 0
            hFig = figure('Name', 'Traces');
        end
        for ishot = 1:Nshot
            ix_min = record2d_obs(ishot).ix_min;
            ix_max = record2d_obs(ishot).ix_max;
            vpm = vp(:, ix_min:ix_max);
            rhom = rho(:, ix_min:ix_max);
            [record2d_pre(ishot), ~, ~] = iso_acoustic2d_propagation_engine ...
                (record2d_pre(ishot), 0, vpm, rhom, dx, dz, dt, nt, ...
                pmlThick, freesurface, precondition, 0, verbose);
            if problem == 2
                record2d_pre(ishot) = compute_adjoint_source2d(record2d_pre(ishot), record2d_obs(ishot), vpm, freq, filter);
            end
            if verbose > 0
                hFig = plot_tiled_images(1, 1, hFig, record2d_pre(ishot).gx/1000, t, 0, 0.5, record2d_pre(ishot).traces_p);
                set_all_axes_labels(hFig, 'Distance / km', 'Time / s');
                set_subplot_titles(['shotID = ', num2str(ishot), ': Observed Traces']);
                colormap(gca, "gray"), drawnow;
            end
        end
        if verbose > 0
            close(hFig);
        end
    end
    grad = 0;
elseif problem == 1
    for ishot = 1:Nshot
        nxx = record2d_pre(ishot).ix_max - record2d_pre(ishot).ix_min + 1;
        record2d_pre(ishot).grad = zeros(nz, nxx);
    end
    if precondition == 1
        for ishot = 1:Nshot
            nxx = record2d_pre(ishot).ix_max - record2d_pre(ishot).ix_min + 1;
            record2d_pre(ishot).precond = zeros(nz, nxx);
        end
    else
        for ishot = 1:Nshot
            record2d_pre(ishot).precond = 0;
        end
    end
    if isParallel == 1
        parfor ishot = 1:Nshot
            ix_min = record2d_obs(ishot).ix_min;
            ix_max = record2d_obs(ishot).ix_max;
            vpm = vp(:, ix_min:ix_max);
            rhom = rho(:, ix_min:ix_max);
            % forward
            [record2d_pre(ishot), wavefields, record2d_pre(ishot).precond] = ...
                iso_acoustic2d_propagation_engine(record2d_pre(ishot), 0, vpm, rhom, ...
                dx, dz, dt, nt, pmlThick, freesurface, precondition, 1, verbose);
            % adjoint source
            record2d_pre(ishot) = compute_adjoint_source2d(record2d_pre(ishot), record2d_obs(ishot), vpm, freq, filter);
            % backward
            [~, ~, record2d_pre(ishot).grad] = iso_acoustic2d_propagation_engine ...
                (record2d_pre(ishot), wavefields, vpm, rhom, dx, dz, dt, nt, ...
                pmlThick, freesurface, precondition, 2, verbose);
            wavefields = 0;
            if precondition == 1
                [record2d_pre(ishot).grad, record2d_pre(ishot).precond] = ...
                    apply_gradient_shot_precondition(record2d_pre(ishot).grad, record2d_pre(ishot).precond);
            end
        end
    else
        if verbose > 0
            hFig1 = figure('Name', 'Adjoint Source');
            hFig2 = figure('Name', 'Pesudo Hessaian and Gradient of Single Shot');
        end
        for ishot = 1:Nshot
            ix_min = record2d_obs(ishot).ix_min;
            ix_max = record2d_obs(ishot).ix_max;
            vpm = vp(:, ix_min:ix_max);
            rhom = rho(:, ix_min:ix_max);
            % forward
            [record2d_pre(ishot), wavefields, record2d_pre(ishot).precond] = ...
                iso_acoustic2d_propagation_engine(record2d_pre(ishot), 0, vpm, rhom, ...
                dx, dz, dt, nt, pmlThick, freesurface, precondition, 1, verbose);
            % adjoint source
            record2d_pre(ishot) = compute_adjoint_source2d(record2d_pre(ishot), record2d_obs(ishot), vpm, freq, filter);
            if verbose > 0
                hFig1 = plot_tiled_images(1, 1, hFig1, record2d_pre(ishot).gx/1000, t, 0, 0.5, record2d_pre(ishot).residual_p);
                set_all_axes_labels(hFig1, 'Distance / km', 'Time / s');
                set_subplot_titles(['shotID = ', num2str(ishot), ': Adjoint Source']);
                colormap(gca, "gray"), drawnow;
            end
            % backward
            [~, ~, record2d_pre(ishot).grad] = iso_acoustic2d_propagation_engine ...
                (record2d_pre(ishot), wavefields, vpm, rhom, dx, dz, dt, nt, ...
                pmlThick, freesurface, precondition, 2, verbose);
            clear wavefields;
            if precondition == 1
                [record2d_pre(ishot).grad, record2d_pre(ishot).precond] = ...
                    apply_gradient_shot_precondition(record2d_pre(ishot).grad, record2d_pre(ishot).precond);
            end
            if verbose > 0
                hFig2 = plot_tiled_images(2, 1, hFig2, x(ix_min:ix_max)/1000, z(1:end)/1000, 1, 0.5, ...
                    record2d_pre(ishot).precond, record2d_pre(ishot).grad);
                set_all_axes_labels(hFig2, 'Distance / km', 'Depth / km');
                set_subplot_titles('Pesudo Hessaian', 'Gradient of V_P');
                drawnow;
            end
        end
        if verbose > 0
            close(hFig1);
            close(hFig2);
        end
    end
    grad = zeros(nz, nx);
    for ishot = 1:Nshot
        ix_min = record2d_obs(ishot).ix_min;
        ix_max = record2d_obs(ishot).ix_max;
        grad(:, ix_min:ix_max) = grad(:, ix_min:ix_max) + record2d_pre(ishot).grad;
    end
    if precondition == 1
        precond = zeros(nz, nx);
        for ishot = 1:Nshot
            ix_min = record2d_obs(ishot).ix_min;
            ix_max = record2d_obs(ishot).ix_max;
            precond(:, ix_min:ix_max) = precond(:, ix_min:ix_max) + record2d_pre(ishot).precond;
        end
        grad = apply_gradient_precondition(grad, precond);
    end
    % vp_min = min(vp, [], 'all');
    % kpass = freq / vp_min;
    % if filter == 1
    %     grad = bandpass2d(grad, dx, dz, kpass, 2*kpass);
    % elseif filter == 0
    %     grad = lowpass2d(grad, dx, dz, kpass, 2*kpass);
    % end
    if verbose > 0
        hFig = figure('Name', 'Pesudo Hessaian and Gradient');
        hFig = plot_tiled_images(2, 1, hFig, x/1000, z/1000, 1, 0.5, precond, grad);
        set_all_axes_labels(hFig, 'Distance / km', 'Depth / km');
        set_subplot_titles('Pesudo Hessaian', 'Gradient of V_P');
        drawnow;
        close(hFig);
    end
end
end
