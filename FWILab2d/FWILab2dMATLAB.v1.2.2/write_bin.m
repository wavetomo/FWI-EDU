function write_bin(filename, data, precision)
%% WRITE_BIN Write data to a binary file.
%
% Syntax:
%   write_bin(filename, data)
%   write_bin(filename, data, precision)
%
% Inputs:
%   filename  - Name of the binary file to write
%   data      - Numeric array to write
%   precision - (Optional) Precision of data type to write (e.g., 'float32', 'int16')
%               Default is 'float32'
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Set default precision
if nargin < 3 || isempty(precision)
    precision = 'float32';
end

% Validate inputs
if ~ischar(filename) && ~isstring(filename)
    error('filename must be a character vector or string scalar.');
end
if ~isnumeric(data)
    error('data must be numeric.');
end

% Open file for binary write
fileID = fopen(filename, 'wb');
if fileID == -1
    error('Failed to open file: %s', filename);
end

% Write data in column-major order
count = fwrite(fileID, data(:), precision);
if count ~= numel(data)
    warning('Only %d of %d elements were written.', count, numel(data));
end

fclose(fileID);
end
