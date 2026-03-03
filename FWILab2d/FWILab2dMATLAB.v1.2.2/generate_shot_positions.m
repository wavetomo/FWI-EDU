function [sx, sz, Nshot] = generate_shot_positions(x, z, fsx, ds, sdepth, Nshot)
%% GENERATE_SHOT_POSITIONS Generate shot coordinates and determine number of shots.
%
%   [sx, sz, Nshot] = numShots(x, z, fsx, ds, sdepth, Nshot)
%
%   Inputs:
%     x         - array of x-coordinates in the model domain (e.g., 0:dx:Lx)
%     z         - array of z-coordinates in the model domain (e.g., 0:dz:Lz)
%     fsx       - x-coordinate of the first shot
%     ds        - interval between adjacent shots in x-direction
%     sdepth    - fixed z-coordinate for all shots
%     Nshot     - desired number of shots (may be adjusted based on limits)
%
%   Outputs:
%     sx        - array of x-positions of shots (length = Nshot)
%     sz        - array of z-positions of shots (same size as sx)
%     Nshot     - actual number of shots used (may be reduced)
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Check if the first shot x-location is within model bounds
if fsx < x(1) || fsx > x(end)
    error('x-coordinate of first shot is out of range!');
end

% Check if the z-location for all shots is within model bounds
if sdepth < z(1) || sdepth > z(end)
    error('z-coordinate of shots is out of range!');
end

% Generate x-locations of shots starting from fsx
sx = fsx:ds:x(end);

% Adjust Nshot if necessary to prevent exceeding available positions
if Nshot > length(sx)
    Nshot = length(sx); % reduce Nshot to available number
elseif Nshot < length(sx)
    sx = sx(1:Nshot); % trim sx array if Nshot is smaller
end

% All shots share the same z-coordinate
sz = sdepth * ones(1, Nshot);
end
