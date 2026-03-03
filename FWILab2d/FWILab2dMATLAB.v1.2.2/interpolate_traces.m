function trs = interpolate_traces(tr, r)
%% INTERPOLSTE_TRACES -Interpolate traces in the time direction by a factor of r using zero-stuffing in time domain.
% Inputs:
%   tr  - [nt × ntr] original trace matrix
%   r   - interpolation factor (integer > 1)
% Output:
%   trs - [ns × ntr] interpolated traces
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

[nt, ntr] = size(tr);
ns = (nt - 1) * r + 1; % Number of samples after interpolation
trs = zeros(ns, ntr); % Allocate output matrix
% Step 1: Insert zeros between samples (zero-stuffing in time domain)
trs(1:r:ns, :) = tr;
% Step 3: Construct a frequency domain low-pass window (simple rectangular window)
% It keeps only the frequencies up to the original Nyquist frequency
n = floor(ns/(2 * r));
win = zeros(ns, 1);
win(1:n) = r;
win(ns-n+2:ns) = r;
trs = real(ifft(win.*fft(trs)));

end