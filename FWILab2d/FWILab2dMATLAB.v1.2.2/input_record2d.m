function [record2d, Nshot] = input_record2d(recordfile, srcfile, nt,dt)
%% INPUT_RECORD2D - Load shot information and attach source wavelet
%
%   [record2d, Nshot] = input_record2d(recordfile, srcfile, nt, dt)
%
%   This function reads in a structure array of shot records, loads the 
%   source wavelet, and attaches it to each shot along with its estimated
%   peak frequency.
%
%   INPUTS:
%       recordfile - path to a .mat file or ASCII file that contains 
%                    a structure array of shots. Each element must include:
%                        .sx, .sz : source position indices
%                        .rx, .rz : receiver position indices
%
%       srcfile    - path to binary file storing the 1D source wavelet
%                    (assumed to be in `float32` format with length = nt)
%
%       nt         - number of time samples in the source wavelet
%       dt         - temporal sampling interval [s]
%
%   OUTPUTS:
%       record2d   - structure array containing per-shot information.
%                    Each element includes:
%                        .sx, .sz       : source indices
%                        .rx, .rz       : receiver indices
%                        .src           : 1D source wavelet (nt x 1)
%                        .fpeak         : estimated peak frequency [Hz]
%
%       Nshot      - total number of shots
%
%   Note:
%       - The peak frequency is estimated using the magnitude spectrum.
%       - This function is typically used before wave propagation or FWI.
%
%   Author: Zhang PingMin
%   Date: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

src = read_bin1d(srcfile,nt);
record2d = importdata(recordfile);
Nshot = length(record2d);

for ishot = 1:Nshot
    record2d(ishot).src = src;
    record2d(ishot).fpeak = calculate_peak_frequency(src, dt);
end

end