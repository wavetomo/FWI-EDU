function set_all_axes_labels(hFig, xLabelStr, yLabelStr)
%% SET_ALL_AXES_LABELS - Set xlabel and ylabel for all axes in a figure
%
% Usage:
%   set_all_axes_labels(hFig, 'X (km)', 'Z (km)')
%
% Inputs:
%   hFig        - Handle to figure (e.g. gcf)
%   xLabelStr   - Label for x-axis
%   yLabelStr   - Label for y-axis
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

if nargin < 1 || isempty(hFig)
    hFig = gcf; % Use current figure if not specified
end
if nargin < 2
    xLabelStr = '';
end
if nargin < 3
    yLabelStr = '';
end

% Find all axes in the figure (excluding legends/colorbars)
axList = findall(hFig, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar');

for k = 1:length(axList)
    ax = axList(k);
    set(ax, 'Units', 'centimeters');
    set(ax, 'TickDir', 'out');
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
    % set(ax, 'XAxisLocation', 'top');
    set(ax, 'FontName', 'Times New Roman');
    set(ax, 'FontSize', 9);
    xlabel(ax, xLabelStr);
    ylabel(ax, yLabelStr);
end
end
