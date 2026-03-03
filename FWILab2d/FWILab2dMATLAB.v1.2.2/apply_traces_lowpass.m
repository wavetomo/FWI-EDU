function trs = apply_traces_lowpass(tr, fpass, fcut, dt)
%% APPLY_TRACES_LOWPASS Applies a low-pass filter to 1D or 2D trace data in frequency domain
%
% Inputs:
%   tr     - input trace(s), size (nt × ntr)
%   fpass  - passband frequency (Hz)
%   fcut   - stopband frequency (Hz), must be > fpass
%   dt     - sampling interval (s)
%
% Output:
%   trs    - filtered trace(s), same size as input
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Validate
if fcut <= fpass
    error('fcut must be greater than fpass');
end

nt = size(tr, 1);
f = generate_frequency(dt, nt) / (2 * pi); % Should return in Hz
f = abs(f); % Symmetric magnitude

% Create frequency-domain window
win = zeros(nt, 1);
taperLen = fcut - fpass;

idx_pass = f < fpass;
idx_taper = (f >= fpass) & (f <= fcut);

win(idx_pass) = 1;
win(idx_taper) = 0.5 * (1 + cos(pi * (f(idx_taper) - fpass) / taperLen));
% f > fcut automatically zero

trs = real(ifft(fft(tr).*win));

end
