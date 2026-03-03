function modp = pad_model2d(mod, nlayer)
%% PAD_MODEL2D Pads a 2D model with absorbing boundary layers.
%
%   modp = MODPAD2D(mod, nlayer) adds 'nlayer' layers of padding
%   around the input model 'mod' for use in absorbing boundary conditions.
%
%   INPUT:
%       mod     - Original 2D model (nz x nx)
%       nlayer  - Number of padding layers
%
%   OUTPUT:
%       modp    - Padded model ((nz+2*nlayer) x (nx+2*nlayer))
%
%   The padding replicates the nearest boundary values from 'mod' to fill
%   the edges, ensuring a smooth transition for absorbing boundaries.
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Validate input
if ndims(mod) ~= 2
    error('Input model must be a 2D matrix.');
end

[nz, nx] = size(mod); % Original model size
Nz = nz + 2 * nlayer;
Nx = nx + 2 * nlayer;

modp = zeros(Nz, Nx); % Initialize padded model

% Insert original model into center
modp(nlayer+1:nlayer+nz, nlayer+1:nlayer+nx) = mod;

% Top padding (replicate first row)
modp(1:nlayer, :) = repmat(modp(nlayer+1, :), nlayer, 1);

% Bottom padding (replicate last row)
modp(nz+nlayer+1:end, :) = repmat(modp(nz+nlayer, :), Nz-(nz + nlayer), 1);

% Left padding (replicate first column)
modp(:, 1:nlayer) = repmat(modp(:, nlayer+1), 1, nlayer);

% Right padding (replicate last column)
modp(:, nx+nlayer+1:end) = repmat(modp(:, nx+nlayer), 1, Nx-(nx + nlayer));

end
