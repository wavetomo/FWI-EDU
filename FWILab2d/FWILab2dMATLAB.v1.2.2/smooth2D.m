function out = smooth2D(in, N, iszdirect)
%% SMOOTH2D
%  One-dimensional cosine-squared smoothing for 2D data
%
%  out = smooth2D(in, nx, nz, N, iszdirect)
%
%  This function applies a cosine-squared smoothing operator to a 2D array
%  along one specified direction (z or x). The smoothing window is truncated
%  near boundaries and re-normalized to preserve amplitude.
%
%  INPUTS:
%      in         - (nz x nx) input 2D data array, indexed as in(iz,ix)
%
%      N          - half window length of smoothing
%                   Total window size = 2*N + 1
%
%      iszdirect  - smoothing direction flag:
%                   = 1 : smooth along z-direction
%                   ~=1 : smooth along x-direction
%
%  OUTPUTS:
%      out        - (nz x nx) smoothed 2D data array
%
%  NOTES:
%      - The smoothing kernel is defined as:
%            w(k) = cos^2( pi * k / (2*N) ),   k = -N : N
%
%      - Near model boundaries, the smoothing window is truncated and
%        the remaining weights are re-normalized.
%
%      - This implementation is equivalent to a spatial low-pass filter.
%        It is NOT equivalent to convolution with zero-padding.
%
%      - Common applications include model smoothing and gradient
%        preconditioning in FWI.
%
%  Author: Chen QuYang
%  Date: 2026-01-23
%  Last Modified: 2026-01-23
%  Copyright (c) 2026 WaveTomo. All rights reserved.

[nz,nx] =size(in);
out = zeros(nz, nx);

%% cosine-squared weights
idx = -N:N;
weight = cos(pi * idx / (2*N)).^2;   % 1 x (2N+1)

%% -------- smooth along z direction --------
if iszdirect == 1

    for iz = 1:nz
        % window range in data
        z1 = max(1, iz - N);
        z2 = min(nz, iz + N);

        % corresponding weight range
        w1 = 1 + (z1 - (iz - N));
        w2 = 1 + (z2 - (iz - N));

        w = weight(w1:w2).';          % column vector
        wsum = sum(w);

        % vectorized over x
        out(iz, :) = (w.' * in(z1:z2, :)) / (wsum + eps);
    end

%% -------- smooth along x direction --------
else

    for ix = 1:nx
        % window range in data
        x1 = max(1, ix - N);
        x2 = min(nx, ix + N);

        % corresponding weight range
        w1 = 1 + (x1 - (ix - N));
        w2 = 1 + (x2 - (ix - N));

        w = weight(w1:w2);            % row vector
        wsum = sum(w);

        % vectorized over z
        out(:, ix) = (in(:, x1:x2) * w.') / (wsum + eps);
    end
end

end