clear all;

%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_5_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_6_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_7_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/');

grids = [3,4,5,6,7,8];
best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));

fprintf('Starting method RBF-G-3D\n');
for i=1:length(grids)

nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 0), 0.1, 4);
end

save('data/best_eps/RBF-G-3D', 'grids', 'best_eps', 'best_scores');


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
    [best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 1), bounds(i,1), bounds(i,2));
end

save('data/best_eps/RBF-MQ-3D', 'grids', 'best_eps', 'best_scores');

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
nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
%[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 2), 0.08, 2);
[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 2), bounds(i,1), bounds(i,2));
end

save('data/best_eps/RBF-G-DF', 'grids', 'best_eps', 'best_scores');

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

nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 3), bounds(i,1), bounds(i,2));
end

save('data/best_eps/RBF-MQ-DF', 'grids', 'best_eps', 'best_scores');
 
plot_best_eps;