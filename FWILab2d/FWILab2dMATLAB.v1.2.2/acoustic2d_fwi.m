clc, clear, close all;

%% path and file
work_path = './marm-ii-fwi-freesurface';
vpfile = './model_geometry/vp_smooth_00221_00601_12.5m.bin';
rhofile = '';
% srcfile = './model_geometry/src_ricker_06001_1.0ms.bin';
srcfile = './model_geometry/ricker_6001_1ms_10Hz_delay0.15ms.bin';
recordfile = './marm-ii-modeling-freesurface/shots_gather_030.mat';

%% parameters
problem = 1;
dx = 12.5;
dz = 12.5;
dt = 1e-3;
nx = 601;
nz = 221;
nt = 6001;
pmlThick = 30;
freesurface = 1;

% freq = [2, 4, 6, 8, 10, 12];
% filter = [0, 0, 0, 0, 0, 2];
% iteration = [10, 5, 5, 5, 5, 5];
% rx = [2, 2, 2, 1, 1, 1];
% rz = [2, 2, 2, 1, 1, 1];
freq = [2, 4, 6, 8, 10];
filter = [0, 0, 0, 0, 2];
iteration = [10, 5, 5, 5, 5];
rx = [2, 2, 2, 1, 1];
rz = [2, 2, 2, 1, 1];

Gardner = 1;
vpmax = 4700;
vpmin = 1500;
rhomax = 2700;
rhomin = 1000;
vp_water = 1500;
rho_water = 1010;
waterdepth = 36;
precondition = 1;

parallel = 1;
verbose = 1;

%% coordinate define
x = (0:nx - 1) * dx;
z = (0:nz - 1) * dz;
t = (0:nt - 1) * dt;

%% Input file
vp = read_bin2d(vpfile, nz, nx);
% vp = importdata('vp_00221_00601_039iter.mat');
if isempty(rhofile) || Gardner == 1
    rho = vp2density(vp);
else
    rho = read_bin2d(rhofile, nz, nx);
end
vp(1:waterdepth, :) = vp_water;
rho(1:waterdepth, :) = rho_water;
slow = vp2slow(vp);

[record2d_obs, Nshot] = input_record2d(recordfile, srcfile, nt, dt);

if verbose > 0
    hFig1 = plot_tiled_images(2, 1, [], x/1000, z/1000, 1, 1, vp, rho);
    set_all_axes_labels(hFig1, 'Distance / km', 'Depth / km');
    set_subplot_titles('V_P', 'Density');
end

if ~exist(work_path, 'dir')
    mkdir(work_path);
end
%% griding observation
record2d_obs = create_grid_observe2d(record2d_obs, dx, dz, nx);

%% inversion
fun1 = zeros(sum(iteration), 1);
fun2 = zeros(sum(iteration), 1);
timeCost = zeros(sum(iteration), 1);
it = 1;
hFig = [];
for k = 1:length(iteration)
    for iter = 1:iteration(k)
        tic;
        % GRADIENT
        [record2d_pre0, grad] = iso_acoustic2d_propagation(vp, rho, record2d_obs, ...
            dx, dz, dt, nt, pmlThick, freesurface, freq(k), filter(k), precondition, problem, parallel, verbose);
        %%% If the built-in smoothing function is used, the inversion results show slight differences compared with the C implementation.
        grad = smooth2D(grad, rz(k), 1);
        grad = smooth2D(grad, rx(k), 0);
        %%% grad = smoothdata(grad, 1, 'gaussian', rz(k));
        %%% grad = smoothdata(grad, 2, 'gaussian', rx(k));
        grad(1:waterdepth, :) = 0;
        % STEP-LENGTH
        fun1(it) = calculate_objective_value(record2d_pre0);
        scalar = compute_model_perturbation(0.01, slow(:), grad(:));
        pertub = -scalar.*grad;
        vp = slow2vp(slow+pertub);
        vp(1:waterdepth, :) = vp_water;
        vp = clip_model(vp, vpmin, vpmax);
        if Gardner == 1
            rho = vp2density(vp); % Gardner relation
            rho(1:waterdepth, :) = rho_water;
            rho = clip_model(rho, rhomin, rhomax);
        end
        [record2d_pre1, ~] = iso_acoustic2d_propagation(vp, rho, record2d_obs, ...
            dx, dz, dt, nt, pmlThick, freesurface, freq(k), filter(k), precondition, 2, parallel, verbose);
        fun2(it) =  calculate_objective_value(record2d_pre1);
        steplen = calculate_step_length(record2d_pre0, record2d_pre1);
        % UPDATE
        slow = slow + steplen * pertub;
        %%% slow = smoothdata(slow, 1, 'gaussian', rz(k));
        %%% slow = smoothdata(slow, 2, 'gaussian', rx(k));
        vp = slow2vp(slow);
        vp(1:waterdepth, :) = vp_water;
        vp = clip_model(vp, vpmin, vpmax);
        if Gardner == 1
            rho = vp2density(vp); % Gardner relation  
            % rho = smoothdata(rho, 1, 'gaussian', rz(k));
            % rho = smoothdata(rho, 2, 'gaussian', rx(k));
            rho(1:waterdepth, :) = rho_water;
            rho = clip_model(rho, rhomin, rhomax);
        end
        slow = vp2slow(vp);
        % OUTPUT
        fileName = sprintf('%s/grad_sp_%05d_%05d_%03diter.mat', work_path, nz, nx, it);
        save(fileName, "grad");
        fileName = sprintf('%s/vp_%05d_%05d_%03diter.mat', work_path, nz, nx, it);
        save(fileName, "vp");
        fileName = sprintf('%s/rho_%05d_%05d_%03diter.mat', work_path, nz, nx, it);
        save(fileName, "rho");
        % PLOT
        fileName = sprintf('Iteration = %03d: Inversion Model', it);
        hFig = plot_tiled_images(3, 1, hFig, x/1000, z/1000, 1, 0.5, vp, rho, grad);
        set_all_axes_labels(hFig, 'Distance / km', 'Depth / km');
        set_subplot_titles('V_P', 'Density', 'Gradient of V_P');
        set(hFig, 'Name', fileName);
        drawnow;
        %
        timeCost(it) = toc;
        fprintf('Iteration %03d: freq = %3.1fHz, steplength = %e, fun1 = %e, fun2 = %e, takes %f seconds.\n', ...
            it, freq(k), steplen, fun1(it), fun2(it), timeCost(it));
        %
        it = it + 1;
    end
end
[h,m,s] = sec2hms(sum(timeCost));
fprintf('Inversion Takes: %d h %d min %.3f sec\n', h, m, s);

[record2d_pre0, ~] = iso_acoustic2d_propagation(vp, rho, record2d_pre0, ...
    dx, dz, dt, nt, pmlThick, freesurface, 0, 0, 0, 0, parallel, verbose);
%% OUTPUT
fileName = sprintf('%s/shots_predicted_gather_%03d.mat', work_path, Nshot);
save(fileName, 'record2d_pre0');
fileName = sprintf('%s/objective_function.mat', work_path);
save(fileName, 'fun1');

%% Shutting down Parfor
if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
    fprintf('parallel is shutting down!\n');
end
