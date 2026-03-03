function data = read_bin2d(fileName, n1, n2, precision)
%% READ_BIN2D Read 2D binary data from file.
%
% Syntax:
%   data = read_bin2d(fileName, n1)
%   data = read_bin2d(fileName, n1, n2)
%   data = read_bin2d(fileName, n1, n2, precision)
%
% Inputs:
%   fileName   - Path to binary file
%   n1         - First dimension (e.g., number of rows)
%   n2           - Second dimension (optional, inferred if omitted)
%   precision  - Data precision (default: 'float32')
%
% Output:
%   data       - n1-by-n2 matrix of binary data
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-11-3
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Set default precision
if nargin < 4 || isempty(precision)
    precision = 'float32';
end

% Get bytes per element for precision
elemSize = get_precision_bytes(precision);

% Infer n2 if not provided
if nargin < 3 || isempty(n2)
    totalElements = count_file_elements(fileName, elemSize);
    if mod(totalElements, n1) ~= 0
        error('Cannot evenly divide file into [%d, ?]. File size mismatch.', n1);
    end
    n2 = totalElements / n1;
end

% Open file
fileID = fopen(fileName, 'rb');
if fileID == -1
    error('Failed to open file: %s', fileName);
end

% Read data
[data, count] = fread(fileID, [n1, n2], precision);
fclose(fileID);

% Check if data size matches expectation
if count ~= n1 * n2
    warning('Expected to read %d elements, but only read %d.', n1*n2, count);
end
end
