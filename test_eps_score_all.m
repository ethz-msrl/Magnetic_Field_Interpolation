clear all;
nodes_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/';

recompute = 0;

if recompute ~= 0
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 0);
    save('data/eps_score/RBF-G-3D_4.mat', 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 1);
    save('data/eps_score/RBF-MQ-3D_4.mat', 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 2);
    save('data/eps_score/RBF-G-DF_4.mat', 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 3);
    save('data/eps_score/RBF-MQ-DF_4.mat', 'eps_v', 'scores', 'cond_numbers');
end

% colors = [0.882,0.416,0.525;
%           0.722,0.541,0.000;
%           0.314,0.639,0.082;
%           0.000,0.678,0.604;
%           0.000,0.604,0.871;
%           0.784,0.427,0.843;
%           0.882,0.416,0.525];

load('data/eps_score/RBF-MQ-3D_4.mat');
cmap = cbrewer('qual', 'Set1', 3);

fh = figure;
colormap(cmap);
plot(eps_v, scores, '-o', 'LineWidth', 2.0);
ylabel('N-RMSE (%)');
yyaxis 'right';
ylabel('condition number');
semilogy(eps_v, cond_numbers, 'LineWidth', 2.0);
xlabel('Shape Parameter');

opt.XLabel = 'Shape Parameter';
opt.YLabel = 'Condition Number';
opt.Markers = {'+', '*', 'x'};
opt.BoxDim = [4.6, 3.];
opt.FontName = 'Helvetica';
opt.AxisLineWidth = 1.5;
opt.FontSize = 11;
opt.Colors = cmap;
opt.Legend = '';

%opt.FileName = 'figures/best_eps.eps';

setPlotProp(opt, fh);

export_fig(fh, 'figures/eps_score_RBF-MQ-3D.pdf');