function mod1 = cap_model(mod0, mod1, cap)
%% CAP_MODEL Limit the relative update of model parameters.
%
%   mod1 = cap2dmodel(mod0, mod1, cap) ensures that the relative update
%   between mod1 and mod0 is less than cap * mod0 (in absolute value).
%
%   INPUTS:
%       mod0 - original model (nz x nx matrix)
%       mod1 - updated model (nz x nx matrix)
%       cap  - scalar, 0 < cap < 1.0
%
%   OUTPUT:
%       mod1 - capped updated model (nz x nx matrix)
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

if cap <= 0 || cap >= 1
    error('cap must be in range (0, 1)');
end

delta = mod1 - mod0;
delta_abs = abs(delta);
delta_cap = cap * abs(mod0);

% Find where update exceeds the allowed cap
mask = delta_abs > delta_cap;

% Normalize direction and scale update
mod1(mask) = mod0(mask) + sign(delta(mask)) .* delta_cap(mask);

end
