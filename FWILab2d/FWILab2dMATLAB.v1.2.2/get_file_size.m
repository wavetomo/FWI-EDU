function fileSize = get_file_size(fileName)
%% GET_FILE_SIZE Get the size of a file in bytes.
%
% Syntax:
%   fileSize = get_file_size(fileName)
%
% Input:
%   fileName - Path to the file as a string or character array
%
% Output:
%   fileSize - Size of the file in bytes (integer)
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

% Try to open the file in read-only binary mode
fileID = fopen(fileName, 'rb');

% Check if file opened successfully
if fileID == -1
    error('Could not open file: %s', fileName);
end

% Move to the end of the file and get its position (i.e., file size)
fseek(fileID, 0, 'eof');
fileSize = ftell(fileID);

% Close the file
fclose(fileID);
end
