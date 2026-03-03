function tr = taper_traces(tr, taperLenTime, taperLenTrace)
%% TRACESTAPER Apply cosine taper to the edges of seismic traces.
%
% Inputs:
%   tr              - Seismic traces (nt x ntr matrix)
%   taperLenTime    - Taper length in the time direction (number of samples)
%   taperLenTrace   - Taper length in the trace/spatial direction (number of traces)
%
% Output:
%   tr              - Tapered seismic traces
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

[nt, ntr] = size(tr);

% Time-domain taper (top and bottom of each trace)
if taperLenTime > 0
    taperLenTime = min(taperLenTime, round(nt/3));
    x = linspace(0, pi, taperLenTime+1)'; % Cosine window
    taperCoef = 0.5 * (1.0 + cos(x)); % From 1 to 0
    tr(1:taperLenTime+1, :) = tr(1:taperLenTime+1, :) .* flipud(taperCoef); % Top
    tr(nt-taperLenTime:nt, :) = tr(nt-taperLenTime:nt, :) .* taperCoef; % Bottom
end

% Trace-domain taper (left and right of each time sample)
if taperLenTrace > 0
    taperLenTrace = min(taperLenTrace, round(ntr/3));
    x = linspace(0, pi, taperLenTrace+1); % Cosine window
    taperCoef = 0.5 * (1.0 + cos(x)); % From 1 to 0
    tr(:, 1:taperLenTrace+1) = tr(:, 1:taperLenTrace+1) .* fliplr(taperCoef); % Left
    tr(:, ntr-taperLenTrace:ntr) = tr(:, ntr-taperLenTrace:ntr) .* taperCoef; % Right
end

end
