function record2d = init_observe2d(x, z, dt, nt, sx, sz, src, offsmin, offsmax, dtr, gdepth, sampleRate)
%% INIT_OBSERVE2D Initialize a structure array for 2D seismic observation geometry.
%
% Inputs:
%   x               - vector of x-coordinates in the model
%   z               - vector of z-coordinates in the model
%   dt              - time step of the simulation (s)
%   nt              - total number of time steps
%   sx, sz          - vectors of shot x and z locations
%   src             - source time function (length = nt)
%   offsmin         - minimum offset from shot to receiver (can be negative)
%   offsmax         - maximum offset from shot to receiver (can be positive)
%   dtr             - receiver spacing in x-direction
%   gdepth          - z-coordinate of receivers
%   sampleRate      - temporal decimation factor (e.g., 2 = sample every 2 steps)
%
% Output:
%   record2d        - array of structures storing observation geometry and data fields
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-11-3
%   Copyright (c) 2025 WaveTomo. All rights reserved.

Nshot = length(sx);

% Generate receiver offsets relative to source
gx_offset = unique([-offsmax:dtr:-offsmin, offsmin:dtr:offsmax]);

% Input validation
if gdepth < z(1) || gdepth > z(end)
    error('gdepth is not in the range of z-coordinates!');
end
if nt ~= length(src)
    error('Length of src must match nt!');
end

% Preallocate structure array
record2d = alloc_record2d(Nshot);

% Loop over each shot to initialize fields
for ishot = 1:Nshot
    record2d(ishot).shotID = ishot;
    record2d(ishot).sx = sx(ishot);
    record2d(ishot).sz = sz(ishot);

    % Compute absolute receiver positions for this shot
    gx_full = gx_offset + sx(ishot);

    % Clip receivers that are outside the model domain
    gx_valid_idx = gx_full >= x(1) & gx_full <= x(end);
    gx_clipped = gx_full(gx_valid_idx);

    % Assign receiver positions
    record2d(ishot).gx = gx_clipped;
    record2d(ishot).ntr = length(gx_clipped);
    record2d(ishot).gz = gdepth * ones(1, record2d(ishot).ntr);

    % Compute time and space sampling info
    record2d(ishot).sampleRate = sampleRate;
    record2d(ishot).ns = floor((nt - 1)/sampleRate) + 1;
    if record2d(ishot).ntr > 1
        record2d(ishot).dtr = mean(diff(gx_clipped));
    else
        record2d(ishot).dtr = 0; % or set to dtr if known
    end
    record2d(ishot).dt = sampleRate * dt;

    % Source info
    record2d(ishot).src = src;
    record2d(ishot).fpeak = calculate_peak_frequency(src, dt); % peak frequency of source wavelet
end
end
