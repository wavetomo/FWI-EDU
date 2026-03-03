function mod_out = slow2vp(mod_in)
%% SLOW2VP Invert model parameters element-wise: mod_out = 1.0 ./ mod_in
%
% Inputs:
%   mod_in  - input model parameters (array)
%
% Outputs:
%   mod_out - inverted model parameters, with zeros preserved (same size as mod_in)
%
% Note:
%   Elements in mod_in that are zero remain zero to avoid division by zero.
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

mod_out = zeros(size(mod_in));
nonzero_idx = mod_in ~= 0;
mod_out(nonzero_idx) = 1.0 ./ mod_in(nonzero_idx);

end
