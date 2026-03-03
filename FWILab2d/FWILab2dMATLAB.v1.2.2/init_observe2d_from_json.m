function [record2d, Nshot] = init_observe2d_from_json(dt, nt, src, jsonfile, sampleRate)
%% INIT_OBSERVE2D Initialize a structure array for 2D seismic observation geometry.
%
% Inputs:
%   dt          - time step of the simulation (s)
%   nt          - total number of time steps
%   src         - source time function (length = nt)
%   
%   sampleRate  - temporal decimation factor (e.g., 2 = sample every 2 steps)
%   jsonfile        - the path of acquisition json file
% Output:
%   record2d    - array of structures storing observation geometry and data fields
% 
%   Author: Chen QuYang
%   Date: 2025-07-22
%   Last Modified: 2025-11-3
%   Copyright (c) 2025 WaveTomo. All rights reserved.
jsonStr = fileread(jsonfile);
dataStruct = jsondecode(jsonStr);
recordStruct = dataStruct.record;

Nshot = length(recordStruct);

if nt ~= length(src)
    error('Length of src must match nt!');
end

% Preallocate structure array
record2d = alloc_record2d(Nshot);

% Loop over each shot to initialize fields
for ishot = 1:Nshot
    record2d(ishot).shotID = ishot;
    record2d(ishot).sx = recordStruct(ishot).sx;
    record2d(ishot).sz = recordStruct(ishot).sz;

    % Assign receiver positions
    record2d(ishot).gx = recordStruct(ishot).gx;
    record2d(ishot).ntr = length(recordStruct(ishot).gx);
    record2d(ishot).gz = recordStruct(ishot).gz;

    % Compute time and space sampling info
    record2d(ishot).sampleRate = sampleRate;
    record2d(ishot).ns = floor((nt - 1)/sampleRate) + 1;
    if record2d(ishot).ntr > 1
        record2d(ishot).dtr = mean(diff(recordStruct(ishot).gx));
    else
        record2d(ishot).dtr = 0; % or set to dtr if known
    end
    record2d(ishot).dt = sampleRate * dt;

    % Source info
    record2d(ishot).src = src;
    record2d(ishot).fpeak = calculate_peak_frequency(src, dt); % peak frequency of source wavelet
end
end
