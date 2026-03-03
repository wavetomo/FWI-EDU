function [D, Beta] = mealCoefs2d(vp, vs, fpeak, alpha, dx, dz, nlayer)
%% MEALCOEFS2D Compute damping and scaling coefficients for PML (Perfectly Matched Layer)
%
% Inputs:
%   vp     - 2D matrix of P-wave velocity [Nz x Nx]
%   vs     - 2D matrix of S-wave velocity [Nz x Nx], can be empty
%   fpeak  - peak frequency of the source (Hz)
%   alpha  - taper power exponent (set -1 for default 2.5)
%   dx     - spatial grid spacing in x-direction
%   dz     - spatial grid spacing in z-direction
%   nlayer - number of PML absorbing layers
%
% Outputs:
%   D      - combined damping coefficient matrix [Nz x Nx]
%   Beta   - combined scaling coefficient matrix [Nz x Nx]
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Get model size
[Nz, Nx] = size(vp);

% PML profile power exponent (polynomial order of profile)
pd = 2;
if alpha == -1.0
    alpha = 2.5;  % default value
end

% Compute reflection coefficient target (based on standard formula)
Rcoef = 10^(-(log10(nlayer) - 1) / log10(2) - 3 - 0.5 * pd);

% Define desired minimum points per wavelength in PML
PPWfd = 10.0;

% Physical thickness of the PML in each direction
thicknessx = nlayer * dx;
thicknessz = nlayer * dz;

% Compute Vs from Vp if not provided (assume Poisson solid)
if isempty(vs)
    Vs = vp / sqrt(3.0);
else
    Vs = vs;
    % If vs == 0, approximate using vp
    zeroIdx = (vs == 0);
    Vs(zeroIdx) = vp(zeroIdx) / sqrt(3.0);
end

% Initialize intermediate damping and scaling coefficient matrices
Dx = zeros(Nz, Nx);
Dz = zeros(Nz, Nx);
Betax = ones(Nz, Nx);
Betaz = ones(Nz, Nx);

%% Apply damping and scaling in top and bottom boundaries
for ix = 1:Nx
    for l = 1:nlayer
        % Weighting function for polynomial profile
        weight = exp(log(2) * (l / nlayer)^alpha) - 1;

        % --- Top boundary
        iz = nlayer - l + 1;
        d0z = -(1 + pd) * vp(iz, ix) * log(Rcoef) / (2 * thicknessz);
        PPWfcs = Vs(iz, ix) / (dz * fpeak);
        betaz = 2 * PPWfcs / PPWfd;
        Dz(iz, ix) = d0z * weight;
        Betaz(iz, ix) = 1 + (betaz - 1) * weight;

        % --- Bottom boundary
        iz = Nz - nlayer + l;
        d0z = -(1 + pd) * vp(iz, ix) * log(Rcoef) / (2 * thicknessz);
        PPWfcs = Vs(iz, ix) / (dz * fpeak);
        betaz = 2 * PPWfcs / PPWfd;
        Dz(iz, ix) = d0z * weight;
        Betaz(iz, ix) = 1 + (betaz - 1) * weight;
    end
end

%% Apply damping and scaling in left and right boundaries
for l = 1:nlayer
    weight = exp(log(2) * (l / nlayer)^alpha) - 1;
    for iz = 1:Nz
        % --- Left boundary
        ix = nlayer - l + 1;
        d0x = -(1 + pd) * vp(iz, ix) * log(Rcoef) / (2 * thicknessx);
        PPWfcs = Vs(iz, ix) / (dx * fpeak);
        betax = 2 * PPWfcs / PPWfd;
        Dx(iz, ix) = d0x * weight;
        Betax(iz, ix) = 1 + (betax - 1) * weight;

        % --- Right boundary
        ix = Nx - nlayer + l;
        d0x = -(1 + pd) * vp(iz, ix) * log(Rcoef) / (2 * thicknessx);
        PPWfcs = Vs(iz, ix) / (dx * fpeak);
        betax = 2 * PPWfcs / PPWfd;
        Dx(iz, ix) = d0x * weight;
        Betax(iz, ix) = 1 + (betax - 1) * weight;
    end
end

%% Combine x- and z-components to get final coefficients
D = sqrt(Dx.^2 + Dz.^2);
Beta = sqrt((Betax - 1.0).^2 + (Betaz - 1.0).^2) + 1.0;

end
