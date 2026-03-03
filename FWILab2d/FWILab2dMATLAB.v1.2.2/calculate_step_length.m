function steplen = calculate_step_length(record2d_pre0, record2d_pre1)
%% CALCULATE_STEP_LENGTH Estimate actual step length based on tentative update
% Inputs:
%   record2d_pre - predicted data before model update
%   record2d_per - predicted data after tentative model update (perturbation)
% Output:
%   steplen      - estimated optimal step length based on misfit reduction
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

Nshot = length(record2d_pre0);
numerator = 0;
denominator = 0;

for ishot = 1:Nshot
    residual = record2d_pre0(ishot).residual_p;
    delta = record2d_pre0(ishot).residual_p - record2d_pre1(ishot).residual_p;

    numerator = numerator + sum(residual.*delta, 'all');
    denominator = denominator + sum(delta .^ 2, 'all');
end

if denominator == 0
    warning('Denominator is zero in step length calculation. Returning step length = 0.');
    steplen = 0;
else
    steplen = numerator / denominator;
end

end
