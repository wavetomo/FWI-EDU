function fun = calculate_objective_value(record2d)
%% CALCULATE_OBJECTIVE_VALUE Compute the objective function value for FWI.
%
%   fun = OBJECTIVE_VALUE(record2d) calculates the misfit function
%   defined as 0.5 * sum of squared residuals for all shots.
%
%   Input:
%       record2d - struct array with field 'residual_p'
%                      residual_p is the difference between observed and
%                      predicted pressure wave data (2D matrix per shot)
%
%   Output:
%       fun      - scalar objective function value (least-squares misfit)
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

if ~isstruct(record2d) || ~isfield(record2d, 'residual_p')
    error('Input must be a struct array with field ''residual_p''.');
end

Nshot = numel(record2d); % Number of shots
fun = 0; % Initialize misfit value

for ishot = 1:Nshot
    residual = record2d(ishot).residual_p;
    fun = fun + sum(residual(:).^2); % Efficient summation
end

fun = 0.5 * fun; % Least-squares cost function
end
