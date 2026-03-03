function fpeak = calculate_peak_frequency(s, dt)
%% CALCULATE_PEAK_FREQUENCY Estimate the peak frequency of a 1D time-domain signal.
%
%   fpeak = calculate_peak_frequency(s, dt) returns the dominant frequency (Hz)
%   of a time-domain signal 's' sampled at time interval 'dt'.
%
%   Inputs:
%       s  - 1D signal vector (time-domain)
%       dt - Sampling interval in seconds
%
%   Output:
%       fpeak - Peak (dominant) frequency in Hz
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Input validation
if nargin < 2
    error('Both signal vector s and sampling interval dt are required.');
end
if ~isvector(s)
    error('Input signal s must be a 1D vector.');
end

% Ensure signal is a column vector
s = s(:);

nt = length(s); % Number of time samples
freq = generate_frequency(dt, nt); % Angular frequency [rad/s]
f_Hz = abs(freq) / (2 * pi); % Convert to Hz

% Compute amplitude spectrum
spec = abs(fft(s)); % Magnitude spectrum

% Limit to positive frequencies only (real signal symmetry)
half = floor(nt/2) + 1;
spec = spec(1:half);
f_Hz = f_Hz(1:half);

% Find the frequency with maximum amplitude
[~, idx] = max(spec);
fpeak = f_Hz(idx);
end
