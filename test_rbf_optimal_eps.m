% This script finds the optimal shape parameter for each RBF method
% it performs a scalar minimization on the average fit over a range of the
% shape parameters
% the results are somewhat dependent on the bounds of the optimization,
% so we had to hand tune the bounds a bit to get what we wanted
% The resulting optimal shape parameters are saved into data/best_eps with
% a file corresponding to the model name
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

clear all;
load_settings;

grids = [3,4,5,6,7,8];
best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));

fprintf('Starting method RBF-G-3D\n');
for i=1:length(grids)
    
    nodes_dataset = sprintf(options.base_dataset, grids(i));
    [best_eps(i), best_scores(i)] = fminbnd(@(x) test_single_eps_score(nodes_dataset, x, 0, options), 0.1, 4);
end

save(sprintf('%s/best_eps/RBF-G-3D', options.data_base_path), 'grids', 'best_eps', 'best_scores');

bounds = [
    0.1,1;
    0.1,1;
    0.1,2;
    0.1,4;
    0.1,8;
    0.1,16];

fprintf('Starting method RBF-MQ-3D\n');
for i=1:length(grids)
    nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
    fprintf('finding best eps for grid:  %d with bounds between[%.1f, %.1f]\n', grids(i), bounds(i,1), bounds(i,2));
    [best_eps(i), best_scores(i)] = fminbnd(@(x) test_single_eps_score(nodes_dataset, x, 1, options), bounds(i,1), bounds(i,2));
end

save(sprint('%s/best_eps/RBF-MQ-3D', options.data_base_path), 'grids', 'best_eps', 'best_scores');

plot(grids, best_eps);

best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));
bounds = [
    0.1,1;
    0.1,1;
    0.1,3;
    0.1,4;
    1,5;
    2,6];
fprintf('Starting method RBF-G-DF\n');
for i=1:length(grids)
    nodes_dataset = sprintf(options.base_dataset, grids(i));
    [best_eps(i), best_scores(i)] = fminbnd(@(x) test_single_eps_score(nodes_dataset, x, 2, options), bounds(i,1), bounds(i,2));
end

save(sprint('%s/best_eps/RBF-G-DF', options.data_base_path), 'grids', 'best_eps', 'best_scores');

best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));

bounds = [
    0.1,1;
    0.1,1;
    0.1,2;
    0.1,4;
    0.1,4;
    0.1,4];

fprintf('Starting method RBF-MQ-DF\n');
for i=1:length(grids)
    nodes_dataset = sprintf(options.base_dataset, grids(i));
    [best_eps(i), best_scores(i)] = fminbnd(@(x) test_single_eps_score(nodes_dataset, x, 3, options), bounds(i,1), bounds(i,2));
end

save(sprint('%s/best_eps/RBF-MQ-DF', options.data_base_path), 'grids', 'best_eps', 'best_scores');
