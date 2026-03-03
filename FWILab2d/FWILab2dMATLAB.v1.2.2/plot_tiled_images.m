function hFig = plot_tiled_images(row, col, hFig, x, z, isEqual, scalar, varargin)
%% PLOT_TILED_IMAGES - Plot multiple 2D images in a tiled layout with optional color control
%
% Syntax:
%   hFig = plot_tiled_images(row, col, hFig, x, z, isEqual, scalar, img1, img2, ...)
%
% Inputs:
%   row     - Number of rows in the tiled layout
%   col     - Number of columns in the tiled layout
%   hFig    - Figure handle or [] to create a new figure
%   x, z    - Coordinate vectors for X and Z axes
%   isEqual - Logical flag, if true uses axis image (default: true)
%   scalar  - Scalar for symmetric clipping range (default: 1)
%   varargin - List of 2D matrices or cells of 2D matrices to be plotted
%
% Output:
%   hFig    - Handle to the created or used figure
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Default value handling
if nargin < 1 || isempty(row), row = 1;
end
if nargin < 2 || isempty(col), col = 1;
end
if nargin < 3 || isempty(hFig) || ~ishandle(hFig) || ~strcmp(get(hFig, 'Type'), 'figure')
    hFig = figure;
else
    figure(hFig);
    clf(hFig);
end
if nargin < 6 || isempty(isEqual), isEqual = 1;
end
if nargin < 7 || isempty(scalar), scalar = 1;
end

% Validate number of plots
nImg = length(varargin);
if nImg ~= row * col
    error('Number of input images (%d) exceeds row * col = %d.', nImg, row*col);
end

% Create layout
tiledlayout(row, col, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot each image
for k = 1:nImg
    nexttile;
    img = varargin{k};
    if isempty(img)
        text(0.5, 0.5, 'Empty', 'HorizontalAlignment', 'center');
        axis off;
        continue;
    end

    if iscell(img), img = img{1};
    end

    imagesc(x, z, img);
    colorbar;

    % Axis scaling
    if isEqual
        axis image;
    else
        axis tight;
    end

    % Color scaling
    clip_min = min(img(:));
    clip_max = max(img(:));

    if clip_min * clip_max < 0
        clip = scalar * max(abs([clip_min, clip_max]));
        clim([-clip, clip]);
        sucolormap(gca, 5); % Custom colormap, ensure it exists
    else
        sucolormap(gca, -5);
    end

    % Optional titles
    % title(sprintf('Image #%d', k));
end
end
