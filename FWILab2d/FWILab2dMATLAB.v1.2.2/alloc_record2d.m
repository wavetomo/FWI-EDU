function record2d = alloc_record2d(Nshot)
%% ALLOC_RECORD2D Allocate and initialize an array of record2d structures.
%
%   record2d = alloc_record2d(Nshot)
%
%   This function pre-allocates a 1-by-Nshot array of structures used for
%   storing seismic shot records, including acquisition geometry and
%   wavefield traces (pressure, vx, vz).
%
%   Input:
%     Nshot    - Number of seismic shots (structures to allocate)
%
%   Output:
%     record2d - 1×Nshot structure array with initialized fields
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Define the default structure template for one shot record
recordTemplate = struct( ...
    'shotID', 0, ... % Shot ID (integer)
    'sx', 0, ... % Source x-coordinate (real units)
    'sz', 0, ... % Source z-coordinate
    'isx', 0, ... % Source grid index x
    'isz', 0, ... % Source grid index z
    'fpeak', 0, ... % Dominant frequency of source wavelet (Hz)
    'src', 0, ... % Adjoint source or wavelet (nt x ntr or vector)
    ...
    'gx', 0, ... % Receiver x-coordinates (array or scalar)
    'gz', 0, ... % Receiver z-coordinates
    'igx', 0, ... % Receiver grid indices x
    'igz', 0, ... % Receiver grid indices z
    ...
    'sampleRate', 0, ... % Recording sample rate
    'ntr', 0, ... % Number of receiver traces
    'ns', 0, ... % Number of time samples
    'dtr', 0, ... % Receiver interval (grid or physical space)
    'dt', 0, ... % Time sampling interval
    ...
    'ix_min', 0, ... % Minimum x-index of recording window
    'ix_max', 0, ... % Maximum x-index of recording window
    ...
    'traces_vx', 0, ... % Recorded vx component (nt x ntr)
    'traces_vz', 0, ... % Recorded vz component (nt x ntr)
    'traces_p', 0, ... % Recorded pressure (nt x ntr)
    ...
    'residual_vx', 0, ... % Residual vx component (nt x ntr)
    'residual_vz', 0, ... % Residual vz component (nt x ntr)
    'residual_p', 0);% Residual pressure (nt x ntr)

% Preallocate structure array using repmat
record2d = repmat(recordTemplate, 1, Nshot);

end
