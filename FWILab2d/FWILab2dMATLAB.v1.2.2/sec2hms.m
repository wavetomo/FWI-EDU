function [h, m, s] = sec2hms(sec)
%% SEC2HMS Convert seconds to hours, minutes and seconds.
%
%   [h, m, s] = SEC2HMS(sec) converts elapsed time given in seconds into
%   hours, minutes and seconds.
%
%   Input:
%       sec  - elapsed time in seconds (scalar, can be non-integer)
%
%   Output:
%       h    - hours (integer)
%       m    - minutes (integer, 0–59)
%       s    - seconds (double, 0–60)
%
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% ---------- input check ----------
if sec < 0
    error('Input time (sec) must be non-negative.');
end

% ---------- calculate h:m:s ----------
h = floor(sec / 3600);
m = floor(mod(sec, 3600) / 60);
s = mod(sec, 60);

end
