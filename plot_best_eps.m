myDir = 'data/best_eps'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct

cmap = cbrewer('qual', 'Set1', length(myFiles));

fh = figure;
hold on;
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  [~,name,~] = fileparts(baseFileName);
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  data = load(fullFileName);
  plot(data.grids, data.best_eps, 'DisplayName', name);
end
legend('show');
hold off;

opt.XLabel = 'Grid Size';
opt.XMinorTick = 'off';
opt.YLabel = 'Best Eps';
opt.Markers = {'o', '+', '*', 'x'};
opt.BoxDim = [4.6, 3.];
opt.FontName = 'Helvetica';
opt.AxisLineWidth = 1.5;
opt.FontSize = 12;
opt.LegendLoc = 'NorthWest';
opt.Colors = cmap;
opt.YGrid = 'on';

setPlotProp(opt);

export_fig(fh, 'figures/best_eps.pdf');