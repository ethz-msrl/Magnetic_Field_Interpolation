% Generates the plot showing how the interpolation accuracy and condition
% number vary with the RBF shape parameter 
%
% This is useful to show the well-known uncertainty relation in RBF
% methods. In short, lower values of the shape parameter or "flatter" basis
% functions should lead to a better interpolation, but they also lead to
% higher condition numbers of the interpolation matrix. There is therefore
% a tradeoff between the conditioning and the interpolation power.
% Ultimatly this leads to an optimal shape parameter for that problem.
%
% Here we only plot the relationship for one of the RBF method and for a 
% single grid size. You can set that with options.rbfs_eps_method and
% options.rbf_eps_grid
%
% The positions are normalized using normalize_positions_maxmin so that the
% optimal parameter does not depend on the scale of the positions but
% rather on the grid size
%
%   Copyright 2020, Samuel Charreyron

clear variables;
load_settings;

data_fn = sprintf('%s/eps_score/%s_%d.mat', options.data_base_path, ...
    options.rbf_eps_method, options.rbf_eps_grid);
load(data_fn);
cmap = cbrewer('qual', 'Set1', 3);

fh = figure;
colormap(cmap);
plot(eps_v, scores, '-o', 'LineWidth', 1.2, 'DisplayName', 'N-RMSE');
ylabel('N-RMSE (%)');
yyaxis 'right';
ylabel('condition number');
semilogy(eps_v, cond_numbers, 'LineWidth', 1.2, 'DisplayName', 'Condition number');
xlabel('Shape Parameter');

% The following options are for the thesis
if strcmp(options.plot_mode,'thesis')
    opt.XLabel = 'Shape Parameter';
    opt.YLabel = 'Condition Number';
    opt.Markers = {'+', '*', 'x'};
    opt.BoxDim = [4.6, 3.];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = options.thesis_font_size;
    opt.Colors = cmap;
else
    % The following options are for IEEE
    opt.XLabel = 'Shape Parameter';
    opt.YLabel = 'Condition Number';
    opt.Markers = {'+', '*', 'x'};
    opt.BoxDim = [1.75, 2.0];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = options.ieee_font_size;
    opt.Colors = cmap;
    opt.LineWidth = [1.2,1.2];

end

setPlotProp(opt, fh);
export_fig(fh, sprintf('figures/eps_score_%s_%d_%s.pdf', options.rbf_eps_method, ...
    options.rbf_eps_grid, options.plot_mode));