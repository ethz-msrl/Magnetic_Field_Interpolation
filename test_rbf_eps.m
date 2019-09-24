clear all;

%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_5_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_6_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_7_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/');

grids = [3,4,5,6,7,8];
% best_eps = zeros(1, length(grids));
% best_scores = zeros(1, length(grids));

% for i=1:length(grids)
% 
% nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
% [best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 0), 1, 50);
% end
% 
% save('data/best_eps/RBF-3D', 'grids', 'best_eps', 'best_scores');
% 
% for i=1:length(grids)
% 
% nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
% [best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 1), 1, 50);
% end
% 
% save('data/best_eps/RBF-MQ-3D', 'grids', 'best_eps', 'best_scores');
% 

% best_eps = zeros(1, length(grids));
% best_scores = zeros(1, length(grids));
% for i=1:length(grids)
% 
% nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
% [best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 2), 1, 100);
% end
% 
% save('data/best_eps/RBF-DF', 'grids', 'best_eps', 'best_scores');

best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));

for i=1:length(grids)

nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x, 3), 3, 100);
end

save('data/best_eps/RBF-MQ-DF', 'grids', 'best_eps', 'best_scores');

plot_best_eps;