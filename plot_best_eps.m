myDir = 'data/best_eps'; %gets directory
myFiles = dir(fullfile(myDir,'*.mat')); %gets all wav files in struct
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
opt.BoxDim = [3.45/2, 2.1];
opt.FontName = 'Helvetica';
opt.AxisLineWidth = 1.5;
opt.FontSize = 8;
opt.LegendLoc = 'NorthWest';
opt.Colors = [0.882,0.416,0.525;
          0.722,0.541,0.000;
          0.314,0.639,0.082;
          0.000,0.678,0.604;
          0.000,0.604,0.871;
          0.784,0.427,0.843;
          0.882,0.416,0.525];

setPlotProp(opt);

export_fig(fh, 'figures/best_eps.pdf');