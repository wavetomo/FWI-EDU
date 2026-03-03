function set_subplot_titles(varargin)
%% SET_SUBPLOT_TITLES - Set titles for all subplots (axes) in current figure
%
% Syntax:
%   set_subplot_titles('Title1', 'Title2', ...)
%
% Example:
%   set_subplot_titles('(a)', '(b)', '(c)', '(d)')
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% 获取当前图窗中的所有 Axes 对象（排除 colorbar 和 legend）
axList = findall(gcf, 'Type', 'Axes', ...
    '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar');

% 翻转顺序使得从上到下、左到右排列（match subplot视觉顺序）
axList = flipud(axList);

nAxes = length(axList);
nTitles = length(varargin);

if nTitles < nAxes
    warning('标题数量少于子图数量，部分子图将无标题。');
elseif nTitles > nAxes
    warning('标题数量多于子图数量，多余的标题将被忽略。');
end

% 设置每个子图的标题
for i = 1:min(nAxes, nTitles)
    ax = axList(i);
    title(ax, varargin{i}, 'FontWeight', 'normal');
end
end
