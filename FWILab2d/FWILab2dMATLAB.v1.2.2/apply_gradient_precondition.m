function grad = apply_gradient_precondition(grad, precond, alpha)
%% APPLY_GRADIENT_PRECONDITION Apply preconditioning to the gradient vector.
%
% Inputs:
%   grad    - gradient vector (row or column vector)
%   precond - preconditioning vector or matrix (same size as grad)
%   alpha   - stabilization factor (scalar, optional, default: 0.005)
%
% Output:
%   grad    - preconditioned gradient vector
%
% Description:
%   Applies a preconditioning to the input gradient by dividing element-wise
%   by the sum of the preconditioning matrix and a scaled RMS value of it,
%   improving numerical stability during optimization.
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

if nargin < 3
    alpha = 0.005; % default stabilization factor
end

% Compute RMS of the preconditioning matrix/vector
precond_rms = sqrt(mean(precond(:).^2));

% Apply preconditioning with stabilization term
grad = grad ./ (precond + alpha * precond_rms);

end
