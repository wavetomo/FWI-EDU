function [grad, precond] = apply_gradient_shot_precondition(grad, precond)
%% APPLY_GRADIENT_SHOT_PRECONDITION Apply preconditioning to a single shot gradient.
%
% Inputs:
%   grad    - gradient vector (same size as precond)
%   precond - preconditioning matrix or vector (same size as grad)
%
% Outputs:
%   grad    - preconditioned gradient
%   precond - updated preconditioning matrix
%
% Description:
%   Applies a stabilization-based preconditioning to the input gradient and
%   updates the preconditioning matrix accordingly. The operation normalizes
%   gradient by a factor related to the average preconditioning value to
%   improve convergence behavior.
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

avg = mean(precond(:)); % compute average value of preconditioning matrix

% Apply preconditioning to the gradient
grad = grad ./ (1.0 + precond ./ (7.0 * avg));

% Update the preconditioning matrix for next iteration
precond = precond ./ (1.0 + precond ./ (7.0 * avg));
end
