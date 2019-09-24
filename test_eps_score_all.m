clear all;
nodes_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/';

% [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 0);
% save('data/eps_score/RBF-3D_4.mat', 'eps_v', 'scores', 'cond_numbers');
% 
% [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 1);
% save('data/eps_score/RBF-MQ-3D_4.mat', 'eps_v', 'scores', 'cond_numbers');
% 
% [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 2);
% save('data/eps_score/RBF-DF_4.mat', 'eps_v', 'scores', 'cond_numbers');
% 
% [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 3);
% save('data/eps_score/RBF-MQ-DF_4.mat', 'eps_v', 'scores', 'cond_numbers');

colors = [0.882,0.416,0.525;
          0.722,0.541,0.000;
          0.314,0.639,0.082;
          0.000,0.678,0.604;
          0.000,0.604,0.871;
          0.784,0.427,0.843;
          0.882,0.416,0.525];

load('data/eps_score/RBF-MQ-3D_4.mat');

fh = figure;
plot(eps_v, scores, '-o', 'LineWidth', 2.0);
ylabel('N-RMSE (%)');
yyaxis 'right';
ylabel('condition number');
plot(eps_v, cond_numbers);
xlabel('Eps');

% opt.XLabel = 'Grid Size';
% opt.YLabel = 'Best Eps';
opt.Markers = {'+', '*', 'x'};
opt.BoxDim = [3.45/2, 2.1];
opt.FontName = 'Helvetica';
opt.AxisLineWidth = 1.5;
opt.FontSize = 8;
opt.Colors = colors;

%opt.FileName = 'figures/best_eps.eps';

setPlotProp(opt, fh);

export_fig(fh, 'figures/eps_score_RBF-MQ-3D.pdf');