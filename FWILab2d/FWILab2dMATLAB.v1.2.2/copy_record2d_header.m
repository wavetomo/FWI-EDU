function record2d_dest = copy_record2d_header(record2d_src)
%% COPY_RECORD2D_HEADER Copy metadata fields from a source record2d array.
%
%   record2d_dest = copy_record2d_header(record2d_src)
%
%   This function copies only the header (metadata) information from
%   record2d_src to a new structure array record2d_dest. It excludes
%   wavefield data fields such as traces_vx, traces_vz, and traces_p.
%
%   Input:
%     record2d_src  - Input structure array (1×Nshot) with source records
%
%   Output:
%     record2d_dest - Output structure array with copied header fields
% 
%   Author: Zhang PingMin
%   Date: 2025-07-22
%   Last Modified: 2025-07-29
%   Copyright (c) 2025 WaveTomo. All rights reserved.

Nshot = length(record2d_src);
record2d_dest = alloc_record2d(Nshot);  % Allocate target structure

for ishot = 1:Nshot
    % Copy source information
    record2d_dest(ishot).shotID    = record2d_src(ishot).shotID;
    record2d_dest(ishot).sx        = record2d_src(ishot).sx;
    record2d_dest(ishot).sz        = record2d_src(ishot).sz;
    record2d_dest(ishot).isx       = record2d_src(ishot).isx;
    record2d_dest(ishot).isz       = record2d_src(ishot).isz;
    record2d_dest(ishot).fpeak     = record2d_src(ishot).fpeak;
    record2d_dest(ishot).src       = record2d_src(ishot).src;

    % Copy receiver geometry
    record2d_dest(ishot).gx        = record2d_src(ishot).gx;
    record2d_dest(ishot).gz        = record2d_src(ishot).gz;
    record2d_dest(ishot).igx       = record2d_src(ishot).igx;
    record2d_dest(ishot).igz       = record2d_src(ishot).igz;

    % Copy sampling information
    record2d_dest(ishot).sampleRate = record2d_src(ishot).sampleRate;
    record2d_dest(ishot).ntr        = record2d_src(ishot).ntr;
    record2d_dest(ishot).ns         = record2d_src(ishot).ns;
    record2d_dest(ishot).dt         = record2d_src(ishot).dt;
    record2d_dest(ishot).dtr        = record2d_src(ishot).dtr;

    % Copy recording window
    record2d_dest(ishot).ix_min     = record2d_src(ishot).ix_min;
    record2d_dest(ishot).ix_max     = record2d_src(ishot).ix_max;
end

end
