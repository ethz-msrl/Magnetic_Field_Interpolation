MODE = 'ieee';

myDir = 'data/best_eps'; %gets directory
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

if strcmp(MODE, 'thesis')
    opt = {};
    opt.XLabel = 'Grid Size';
    opt.XMinorTick = 'off';
    opt.YLabel = 'Best Eps';
    opt.Markers = {'o', '+', '*', 'x'};
    opt.BoxDim = [4.6, 3.];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = 12;
    opt.LegendLoc = 'NorthWest';
    opt.Colors = c;
    opt.YGrid = 'on';

else
    opt = {};
    opt.XLabel = 'Grid Size';
    opt.XMinorTick = 'off';
    opt.YLabel = 'Best Eps';
    opt.Markers = {'o', '+', '*', 'x'};
    opt.BoxDim = [1.75, 1.4];
    opt.FontName = 'Helvetica';
    opt.AxisLineWidth = 1.5;
    opt.FontSize = 6;
    opt.LegendLoc = 'NorthWest';
    opt.Colors = c;
    opt.YGrid = 'on';


end

setPlotProp(opt);

export_fig(fh, sprintf('figures/best_eps_%s.pdf', MODE));