clear all;
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