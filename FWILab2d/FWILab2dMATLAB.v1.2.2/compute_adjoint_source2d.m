function record2d_pre = compute_adjoint_source2d(record2d_pre, record2d_obs, vp, freq, filter)
%% COMPUTE_ADJOINT_SOURCE2D Compute adjoint source from observed and predicted data.
%
%   record2d_pre = compute_adjoint_source2d(record2d_pre, record2d_obs, vp, freq, filter)
%
%   This function computes the adjoint source term used in adjoint-state
%   methods by subtracting the observed and predicted pressure traces.
%   Optional time and space domain tapering and lowpass filtering are applied.
%
%   Inputs:
%     record2d_pre   - structure with predicted data and acquisition parameters:
%                      .traces_p : predicted pressure traces [nt x ntr]
%                      .igx      : receiver x-indices (1-based) [1 x ntr]
%                      .igz      : receiver z-indices (1-based) [1 x ntr]
%                      .dt       : time sampling interval
%                      .dtr      : spatial sampling interval along receiver line
%     record2d_obs   - structure with observed pressure traces:
%                      .traces_p : observed pressure traces [nt x ntr]
%     vp             - 2D velocity model [nz x nx]
%     freq           - central frequency (Hz) of the bandpass filter
%     filter         - flag (0 or 1), apply filtering if nonzero
%
%   Output:
%     record2d_pre.residual_p - adjoint source after differencing and optional filtering [nt x ntr]
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

%% Step 1: Compute data residual
record2d_pre.residual_p = record2d_pre.traces_p - record2d_obs.traces_p;

%% Step 2: Estimate taper parameters based on wave propagation
% Convert (igx, igz) to linear indices for velocity sampling
ig_lin = (record2d_pre.igx - 1) * size(vp, 1) + record2d_pre.igz;
vp_avg = mean(vp(ig_lin), 'all'); % Average velocity at receivers
if freq == 0
    freq = record2d_pre.fpeak;
end
kpass = freq / vp_avg; % Wavenumber cutoff

% Calculate taper window lengths
taperLenTime = round((2.0 / freq)/record2d_pre.dt) + 1;
taperLenTrace = round((2.0 / kpass)/record2d_pre.dtr) + 1;

%% Step 3: Apply filtering if requested
if filter == 0 % lowpass
    % Time-domain taper before filter
    record2d_pre.residual_p =  taper_traces(record2d_pre.residual_p, taperLenTime, 0);

    % Time-domain lowpass filtering
    record2d_pre.residual_p = apply_traces_lowpass(record2d_pre.residual_p, freq, 2.0*freq, record2d_pre.dt);

    % Space-domain taper
    record2d_pre.residual_p = taper_traces(record2d_pre.residual_p, 0, taperLenTrace);

    % Space-domain lowpass filtering (transpose in/out)
    record2d_pre.residual_p = apply_traces_lowpass(record2d_pre.residual_p', kpass, 2.0*kpass, record2d_pre.dtr)';
elseif filter == 1 % bandpass
    % Time-domain taper before filter
    record2d_pre.residual_p =  taper_traces(record2d_pre.residual_p, taperLenTime, 0);

    % Time-domain bandpass filtering
    record2d_pre.residual_p = apply_traces_bandpass(record2d_pre.residual_p, freq, 2.0*freq, record2d_pre.dt);

    % Space-domain taper
    record2d_pre.residual_p = taper_traces(record2d_pre.residual_p, 0, taperLenTrace);

    % Space-domain bandpass filtering (transpose in/out)
    record2d_pre.residual_p = apply_traces_bandpass(record2d_pre.residual_p', kpass, 2.0*kpass, record2d_pre.dtr)';
end
%% Step 4: Final taper in both time and space domains
record2d_pre.residual_p = taper_traces(record2d_pre.residual_p, taperLenTime, taperLenTrace);

end
