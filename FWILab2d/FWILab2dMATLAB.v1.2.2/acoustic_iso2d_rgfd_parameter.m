function [K, invrho, invrhox, invrhoz] = acoustic_iso2d_rgfd_parameter(vp, rho, dx, dz)
%% model's parameter at grid node for second-order acoustic wave equation
% Input:
%   vp    - 2D velocity model (nz x nx)
%   rho   - 2D density model (nz x nx)
%   dx    - grid spacing in x-direction
%   dz    - grid spacing in z-direction
%
% Output:
%   K         - bulk modulus, K = rho * vp^2
%   inv_rho   - reciprocal of density
%   invrhox   - first derivative of 1/rho along x (using optimized stencil)
%   invrhoz   - first derivative of 1/rho along z (using optimized stencil)
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

[nz, nx] = size(vp);

K = rho .* vp .* vp;
invrho = 1 ./ rho;
% [invrhox, invrhoz] = gradient(invrho, dx, dz);
[invrhox, invrhoz] = calculate_gradient2d(invrho, dx, dz);
C1 = 0.91067892; % omptimized cofficients for 1st order derivative based on L0 norm
C2 = -0.34187892;
C3 = 0.13833962;
C4 = -0.04880710;
C5 = 0.01302148;
C6 = -0.00199047;

invdx = 1.0 / dx;
invdz = 1.0 / dz;

for ix = 7:nx - 6
    for iz = 7:nz - 6
        invrhox(iz, ix) = ( ...
            C1 * (invrho(iz, ix+1) - invrho(iz, ix-1)) + ...
            C2 * (invrho(iz, ix+2) - invrho(iz, ix-2)) + ...
            C3 * (invrho(iz, ix+3) - invrho(iz, ix-3)) + ...
            C4 * (invrho(iz, ix+4) - invrho(iz, ix-4)) + ...
            C5 * (invrho(iz, ix+5) - invrho(iz, ix-5)) + ...
            C6 * (invrho(iz, ix+6) - invrho(iz, ix-6))) * invdx;
        invrhoz(iz, ix) = ( ...
            C1 * (invrho(iz+1, ix) - invrho(iz-1, ix)) + ...
            C2 * (invrho(iz+2, ix) - invrho(iz-2, ix)) + ...
            C3 * (invrho(iz+3, ix) - invrho(iz-3, ix)) + ...
            C4 * (invrho(iz+4, ix) - invrho(iz-4, ix)) + ...
            C5 * (invrho(iz+5, ix) - invrho(iz-5, ix)) + ...
            C6 * (invrho(iz+6, ix) - invrho(iz-6, ix))) * invdz;
    end
end

end