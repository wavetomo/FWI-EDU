function data = read_bin1d(fileName, n, precision)
%% READ_BIN1D Read 1D binary data from file.
%
% Syntax:
%   data = read_bin1d(fileName)
%   data = read_bin1d(fileName, n)
%   data = read_bin1d(fileName, n, precision)
%
% Inputs:
%   fileName  - Path to binary file
%   n         - Number of elements to read (optional)
%               If omitted, it is inferred from file size
%   precision - Data precision, e.g., 'float32', 'double' (default: 'float32')
%
% Output:
%   data      - Column vector of read data
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-11-3
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Set default precision
if nargin < 3 || isempty(precision)
    precision = 'float32';
end

% Infer number of elements if not provided
if nargin < 2 || isempty(n)
    elemSize = get_precision_bytes(precision);
    n = count_file_elements(fileName, elemSize);
end

% Open file
fid = fopen(fileName, 'rb');
if fid == -1
    error('Failed to open file: %s', fileName);
end

% Read data
[data, count] = fread(fid, [n, 1], precision);
if count ~= n
    warning('Expected to read %d elements, but only read %d.', n, count);
end

fclose(fid);
end
