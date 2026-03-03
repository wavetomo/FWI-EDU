function f = generate_frequency(dt, nt)
%% GENERATE_FREQUENCY Generate angular frequency vector for uniformly sampled signals.
%
% Inputs:
%   dt - sampling interval in time domain (seconds)
%   nt - number of time samples
%
% Output:
%   f  - frequency vector (angular frequency, rad/s), size [nt x 1]
%
% Notes:
%   - This generates the angular frequencies corresponding to the FFT
%     output ordering in MATLAB for a real-valued time series.
%   - Frequencies run from 0 up to Nyquist, then negative frequencies.
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

f = zeros(nt, 1); % preallocate frequency vector
df = 2 * pi / (nt * dt); % frequency resolution (rad/s)
nf = floor(nt/2) + 1; % index of Nyquist frequency

% Set zero frequency explicitly
f(1) = 0;

% Fill positive frequencies
for it = 2:nf
    f(it) = (it - 1) * df;
    % Assign corresponding negative frequencies symmetrically
    if nt - it + 2 > nf
        f(nt-it+2) = -f(it);
    end
end

end
