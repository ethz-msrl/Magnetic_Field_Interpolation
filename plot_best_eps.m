load_settings;

myDir = sprintf('%s/best_eps', options.data_base_path); %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct

load('data/colors.mat');

fh = figure;
hold on;
% because setPlotProp stupidly overwrites colors this stores the assigned
% colors
c = [];
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  [~,name,~] = fileparts(baseFileName);
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  data = load(fullFileName);
  plot(data.grids, data.best_eps, 'DisplayName', name, 'Color', colors(name));
  c = [c; colors(name)'];
end

legend('show');
hold off;

if strcmp(options.plot_mode, 'thesis')
    opt = {};
    opt.XLabel = 'Grid Size';
    opt.XMinorTick = 'off';
    opt.YLabel = 'Best Eps';
    opt.Markers = {'o', '+', '*', 'x'};
    opt.BoxDim = [4.6, 3.];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = options.thesis_font_size;
    opt.LegendLoc = 'NorthWest';
    opt.Colors = c;
    opt.YGrid = 'on';

else
    opt = {};
    opt.XLabel = 'Grid Size';
    opt.XMinorTick = 'off';
    opt.YLabel = 'Best Eps';
    opt.Markers = {'o', '+', '*', 'x'};
    opt.BoxDim = [1.75, 2.0];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = options.ieee_font_size;
    opt.LegendLoc = 'NorthWest';
    opt.Colors = c;
    opt.YGrid = 'on';
    opt.LineWidth = [1.2, 1.2, 1.2, 1.2];
end

setPlotProp(opt);

export_fig(fh, sprintf('figures/best_eps_%s.pdf', options.plot_mode));