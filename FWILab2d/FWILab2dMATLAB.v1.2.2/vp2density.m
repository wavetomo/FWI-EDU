function rho = vp2density(vp)
%% VP2DENSITY Convert P-wave velocity (vp) to density (rho) by Gardner relation
%
%   rho = VP2DENSITY(vp) computes the density rho (kg/m^3) from
%   the input P-wave velocity vp (m/s) using an empirical power-law.
%
%   The empirical relationship used is:
%       rho = 310 * vp^0.25
%
%   Input:
%       vp  - P-wave velocity in m/s (scalar, vector, or matrix)
%
%   Output:
%       rho - Estimated density in kg/m^3 (same size as vp)
%
%   Notes:
%   - This is a simplified empirical model and may not be valid
%     for all rock types or velocity ranges.
%   - Typical use is for sedimentary basins in geophysical modeling.
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Input validation
if nargin < 1
    error('Function vp2density requires an input vp.');
end

if any(vp(:) <= 0)
    error('All input P-wave velocities must be positive.');
end

% Apply the empirical formula
rho = 310 * vp.^0.25;
end
