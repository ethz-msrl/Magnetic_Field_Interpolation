clear all;

%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_5_h5/');
% plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_6_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_7_h5/');
%plot_eps_score('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/');

grids = [4,5,6,7,8];
best_eps = zeros(1, length(grids));
best_scores = zeros(1, length(grids));

for i=1:length(grids)

nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grids(i));
[best_eps(i), best_scores(i)] = fminbnd(@(x) rbf_eps_score(nodes_dataset, x), 1, 50);
end

save('best_eps_size_rbf', 'grids', 'best_eps', 'best_scores');
