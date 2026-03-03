
function [dfdx, dfdz] = calculate_gradient2d(f, dx, dz)
% CALCULATE_GRADIENT2D compute gradients of a 2D array (like Matlab gradient)
%
% Input:
%   f[nz][nx]   - input 2D array
%                 first dimension : z direction (rows)
%                 second dimension: x direction (columns)
%   dx          - grid spacing along x direction
%   dz          - grid spacing along z direction
%
% Output:
%   dfdx[nz][nx] - gradient along x direction
%   dfdz[nz][nx] - gradient along z direction
%
% Finite-difference scheme:
%   interior points : central difference
%   boundary points : forward / backward difference

[nz, nx] = size(f);

dfdx = zeros(nz, nx);
dfdz = zeros(nz, nx);

%% X-direction derivative (along columns)
for ix = 1:nx
    for iz = 1:nz
        if ix == 1
            % forward difference
            dfdx(iz, ix) = (f(iz, ix+1) - f(iz, ix)) / dx;
        elseif ix == nx
            % backward difference
            dfdx(iz, ix) = (f(iz, ix) - f(iz, ix-1)) / dx;
        else
            % central difference
            dfdx(iz, ix) = (f(iz, ix+1) - f(iz, ix-1)) / (2.0 * dx);
        end
    end
end

%% Z-direction derivative (along rows)
for ix = 1:nx
    for iz = 1:nz
        if iz == 1
            % forward difference
            dfdz(iz, ix) = (f(iz+1, ix) - f(iz, ix)) / dz;
        elseif iz == nz
            % backward difference
            dfdz(iz, ix) = (f(iz, ix) - f(iz-1, ix)) / dz;
        else
            % central difference
            dfdz(iz, ix) = (f(iz+1, ix) - f(iz-1, ix)) / (2.0 * dz);
        end
    end
end
end
