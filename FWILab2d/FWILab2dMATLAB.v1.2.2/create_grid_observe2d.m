function record2d = create_grid_observe2d(record2d, dx, dz, nx)
%% CREATE_GRID_OBSERVE2D Convert real-world coordinates to grid indices in 2D.
%
%   record2d = create_grid_observe2d(record2d, dx, dz, nx)
%
%   This function converts the real-world coordinates of sources and
%   receivers to integer grid indices based on the spatial sampling
%   intervals dx and dz. It also sets local subgrid indices for simulation.
%
%   Input:
%     record2d - structure array containing shot gathers
%     dx       - spatial grid spacing in x-direction
%     dz       - spatial grid spacing in z-direction
%     nx       - total number of x-grid points (for bounding ix_max)
%
%   Output:
%     record2d - updated with grid-based indices and local index ranges
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

Nshot = length(record2d);

for ishot = 1:Nshot
    % Validate gx/gz are vectors
    if ~isvector(record2d(ishot).gx)
        error('gx must be a vector');
    end
    if ~isvector(record2d(ishot).gz)
        error('gz must be a vector');
    end
    
    % Convert real coordinates (in meters) to global grid indices
    isx = floor(record2d(ishot).sx/dx) + 1;
    isz = floor(record2d(ishot).sz/dz) + 1;
    igx = floor(record2d(ishot).gx/dx) + 1;
    igz = floor(record2d(ishot).gz/dz) + 1;

    % Define local x-axis simulation window
    ix_min = max(igx(1)-30, 1); % pad 30 grid points to the left
    ix_max = min(igx(end)+30, nx); % pad 30 grid points to the right

    % Adjust global indices to local subgrid indices
    record2d(ishot).isx = isx - ix_min + 1;
    record2d(ishot).isz = isz;
    record2d(ishot).igx = igx - ix_min + 1;
    record2d(ishot).igz = igz;

    % Store subgrid x-boundaries
    record2d(ishot).ix_min = ix_min;
    record2d(ishot).ix_max = ix_max;
end

end
