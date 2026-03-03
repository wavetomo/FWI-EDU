clc, 
clear, 
close all;
%% Used to observed data

%% path and file
work_path = './marm-ii-modeling-freesurface';
vpfile = './model_geometry/vp_00221_00601_12.5m.bin';
rhofile = '';
% srcfile = './model_geometry/src_ricker_06001_1.0ms.bin';
srcfile = './model_geometry/ricker_6001_1ms_10Hz_delay0.15ms.bin';
jsonfile = './model_geometry/acquisition.json';
%% parameters
problem = 0; % 0:forward 1:inversion
dx = 12.5;
dz = 12.5;
dt = 1e-3;
nx = 601;
nz = 221;
nt = 6001;
pmlThick = 30;
freesurface = 1;

Nshot = 30;
fsx = 0; % first sx
ds = 60 * dx;
sdepth = 2 * dz;

sampleRate = 2;
offsmin = dx;
offsmax = 600 * dx;
dtr = dx;
gdepth = 2 * dz;

Gardner = 1;
vp_water = 1500;
rho_water = 1010;
waterdepth = 36; % unit: grid

parallel = 0; % 1: turn on 
maxNumCompThreads(1);
verbose = 1; 

%% coordinate define
x = (0:nx - 1) * dx;
z = (0:nz - 1) * dz;
t = (0:nt - 1) * dt;

%% Input file
src = read_bin1d(srcfile, nt);
vp = read_bin2d(vpfile, nz, nx);
if isempty(rhofile) || Gardner == 1
    rho = vp2density(vp);
else
    rho = read_bin2d(rhofile, nz, nx);
end
vp(1:waterdepth, :) = vp_water;
rho(1:waterdepth, :) = rho_water;
slow = vp2slow(vp);

if verbose > 0
    hFig = plot_tiled_images(2, 1, [], x/1000, z/1000, 1, 1, vp, rho);
    set_all_axes_labels(hFig, 'Distance / km', 'Depth / km');
    set_subplot_titles('V_P', 'Density');
    drawnow;
end

%% observation define
if strcmp(jsonfile,'')
    [sx, sz, Nshot] = generate_shot_positionsx(x, z, fsx, ds, sdepth, Nshot);
    record2d = init_observe2d(x, z, dt, nt, sx, sz, src, offsmin, offsmax, dtr, gdepth, sampleRate);
else
    [record2d, Nshot] = init_observe2d_from_json(dt, nt, src, jsonfile, sampleRate);
end

%% griding observation
record2d = create_grid_observe2d(record2d, dx, dz, nx);

%% modeling
tic;
[record2d, ~] = iso_acoustic2d_propagation(vp, rho, record2d, ...
    dx, dz, dt, nt, pmlThick, freesurface, 0, 0, 0, problem, parallel, verbose);
toc;

%% gather output
fileName = sprintf('%s/shots_gather_%03d.mat', work_path, Nshot);
folder = fileparts(fileName);
if ~exist(folder, 'dir')
    mkdir(folder);
end
save(fileName, 'record2d');
