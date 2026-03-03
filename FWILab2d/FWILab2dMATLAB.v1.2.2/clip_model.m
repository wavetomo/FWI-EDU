function mod = clip_model(mod, mod_min, mod_max)
%% CLIP_MODEL Clip 2D model parameter values within specified range.
%
%   mod = clip2dmodel(mod, mod_min, mod_max) clips the values in the
%   2D matrix 'mod' so that no value is smaller than mod_min or greater than mod_max.
%
%   INPUTS:
%       mod      - 2D matrix (nz x nx), the model parameter
%       mod_min  - scalar, minimum allowed value
%       mod_max  - scalar, maximum allowed value
%
%   OUTPUT:
%       mod      - clipped 2D matrix
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

mod(mod > mod_max) = mod_max;
mod(mod < mod_min) = mod_min;

end
