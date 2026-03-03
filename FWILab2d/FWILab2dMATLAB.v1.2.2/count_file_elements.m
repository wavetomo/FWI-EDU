function numElements = count_file_elements(fileName, elementSize)
%% COUNT_FILE_ELEMENTS  Calculate number of elements in a binary file.
%
% Syntax:
%   numElements = count_file_elements(fileName)
%   numElements = count_file_elements(fileName, elementSize)
%
% Inputs:
%   fileName     - Path to the binary file
%   elementSize  - Size of each element in bytes (default is 4, i.e., float32)
%
% Output:
%   numElements  - Number of elements in the file
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Default element size if not provided
if nargin < 2
    elementSize = 4;  % float32 default
end

% Validate elementSize
if elementSize <= 0
    error('Element size must be a positive integer.');
end

% Get total file size in bytes
fileSizeInBytes = get_file_size(fileName);

% Calculate number of elements
numElements = fileSizeInBytes / elementSize;

% Ensure result is an integer
if mod(numElements, 1) ~= 0
    warning('File size is not a multiple of element size. Result is not an integer.');
end

end
