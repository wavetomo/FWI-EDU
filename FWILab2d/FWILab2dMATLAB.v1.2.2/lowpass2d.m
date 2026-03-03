function s_filter = lowpass2d(s, dx, dz, kpass, kcut)
%% LOWPASS2D 2D low-pass filter in the wavenumber domain
%
% Inputs:
%   s      - 2D scalar field (nz x nx)
%   dx     - spatial grid spacing in x-direction
%   dz     - spatial grid spacing in z-direction
%   kpass  - passband cutoff wavenumber (no attenuation for k < kpass)
%   kcut   - stopband cutoff wavenumber (full attenuation for k > kcut)
%
% Output:
%   s      - filtered 2D field
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Validate cutoff settings
if kcut <= kpass
    error('kcut must be greater than kpass');
end

% Get size of input
[nz, nx] = size(s);

% taper data
taperLenX = round(4/(kpass * dx));
taperLenZ = round(4/(kpass * dz));
s = taper_traces(s, taperLenZ, taperLenX);

% padding zeros to data
padLen = 10 * max([taperLenX, taperLenZ]);
spad = pad_model2d(s, padLen);
[Nz, Nx] = size(spad);

% Frequency vectors (cycles/m) in x and z
kx = generate_frequency(dx, Nx)' / (2 * pi); % Convert to spatial wavenumber [1/m]
kz = generate_frequency(dz, Nz) / (2 * pi);

% Create 2D meshgrid of frequencies
k = sqrt(kx.^2+kz.^2); % Magnitude of 2D wavenumber

% Define smooth taper window
win = zeros(Nz, Nx);
taperLen = kcut - kpass;

% Build cosine taper in annulus between kpass and kcut
win(k < kpass) = 1.0; % fully pass
inTaper = (k >= kpass) & (k <= kcut);
win(inTaper) = 0.5 * (1 + cos((k(inTaper) - kpass)/taperLen*pi));
% beyond kcut is already zero

% Apply filter in frequency domain
s_filter = real(ifft2(win.*fft2(spad)));
s_filter = s_filter(padLen+1:padLen+nz, padLen+1:padLen+nx);

end
