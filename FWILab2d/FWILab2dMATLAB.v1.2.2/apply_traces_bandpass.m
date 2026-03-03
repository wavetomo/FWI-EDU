function trs = apply_traces_bandpass(tr, fpass, fcut, dt)
%% APPLY_TRACES_BANDPASS Applies a band-pass filter to 1D or 2D trace data in frequency domain
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
%   Author: Chen Quyang
%   Date: 2025-11-4
%   Last Modified: 2025-11-4
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
lowTaperLen = fpass - 0;
highTaperLen = fcut - fpass;

idx_low_taper = f < fpass;
idx_high_taper = (f >= fpass) & (f <= fcut);

win(idx_low_taper) = 0.5 * (1 + cos(pi * (fpass - f(idx_low_taper)) / lowTaperLen));
win(idx_high_taper) = 0.5 * (1 + cos(pi * (f(idx_high_taper) - fpass) / highTaperLen));
% f > fcut automatically zero

trs = real(ifft(fft(tr).*win));

end
