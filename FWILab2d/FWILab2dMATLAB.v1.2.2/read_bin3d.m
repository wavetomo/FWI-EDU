function data = read_bin3d(fileName, n1, n2, n3, precision)
%% READ_BIN3D Read 3D binary data from file into [n1, n2, n3] shape
%
% Usage:
%   data = read_bin3d(fileName, n1, n2, n3)
%   data = read_bin3d(fileName, n1, n2, [], 'float32')
%
% Inputs:
%   fileName  - Path to binary file
%   n1        - First dimension (e.g. z)
%   n2        - Second dimension (e.g. x)
%   n3        - Third dimension (e.g. shot or time)
%   precision - Data type, default is 'float32'
%
% Output:
%   data      - 3D matrix [n1, n2, n3]
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Default precision
if nargin < 5 || isempty(precision)
    precision = 'float32';
end

% Bytes per element
elemSize = get_precision_bytes(precision);

% Infer n3 if not provided
if nargin < 4 || isempty(n3)
    totalElems = count_file_elements(fileName, elemSize);
    if mod(totalElems, n1*n2) ~= 0
        error('Data size mismatch: cannot reshape to [%d, %d, ?]', n1, n2);
    end
    n3 = totalElems / (n1 * n2);
end

% Open file
fileID = fopen(fileName, 'rb');
if fileID == -1
    error('Failed to open file: %s', fileName);
end

% Read and reshape
[rawData, count] = fread(fileID, n1*n2*n3, precision);
fclose(fileID);

if count ~= n1 * n2 * n3
    warning('Expected %d elements, but read %d.', n1*n2*n3, count);
end

% Reshape to [n1, n2, n3]
data = reshape(rawData, [n1, n2, n3]);
end
