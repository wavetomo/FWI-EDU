function scalar = compute_model_perturbation(fraction, slowness, grad)
%% MODEL_PERTURBATION Compute model perturbation based on gradient statistics
%
%   [out, scalar] = model_perturbation(fraction, slowness, gradient)
%   calculates a scaled perturbation of the slowness model such that the
%   perturbation follows the direction of the gradient but is capped to
%   include 99% of gradient variation.
%
%   INPUTS:
%       fraction  - relative scale factor (e.g., 0.01)
%       slowness  - slowness model, 1D array (length N)
%       gradient  - gradient corresponding to the model, 1D array
%
%   OUTPUTS:
%       scalar    - computed scaling factor for the gradient
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.


len = length(slowness);
cap = 0.99;

% Compute histogram range
MaxG = max(grad);
MinG = min(grad);
nbins = 1001;
interval = (MaxG - MinG) / (nbins - 1);

% Initialize histogram
histogram = zeros(1, nbins);

% Fill histogram
idx = round((grad - MinG)/interval) + 1;
idx(idx < 1) = 1;
idx(idx > nbins) = nbins;
for i = 1:len
    histogram(idx(i)) = histogram(idx(i)) + 1;
end

% Find peak bin
[~, pos] = max(histogram);
icur = pos;
icur2 = pos;

% Accumulate probability mass until it exceeds cap
accumulator = histogram(pos) / len;
while accumulator < cap
    if icur < nbins
        icur = icur + 1;
        accumulator = accumulator + histogram(icur) / len;
    end
    if icur2 > 1
        icur2 = icur2 - 1;
        accumulator = accumulator + histogram(icur2) / len;
    end
end

% Compute scaling range
g1 = abs(MinG+(icur - 1)*interval);
g2 = abs(MinG+(icur2 - 1)*interval);
scalar_max = max(g1, g2);

% Compute average slowness
avg_slowness = mean(slowness);

% Final scaling factor
scalar = (avg_slowness * fraction) / scalar_max;

end
