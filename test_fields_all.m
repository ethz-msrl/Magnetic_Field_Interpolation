clear all;

RECOMPUTE = 0;

grid_sizes = {3,4,5,6};
noise_std = 0;

if RECOMPUTE ~= 0
     disp('Testing RBF 3D');
    load('data/best_eps/RBF-G-3D', 'best_eps');
    test_fields('RBF-G-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);
   
    disp('Testing RBF Multiquadric 3D');
    load('data/best_eps/RBF-MQ-3D', 'best_eps');
    test_fields('RBF-MQ-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing RBF Div-free');
    load('data/best_eps/RBF-G-DF', 'best_eps');
    test_fields('RBF-G-DF', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing RBF Multiquadric Div-free');
    load('data/best_eps/RBF-MQ-DF');
    test_fields('RBF-MQ-DF', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing Tricubic 3D');
    test_fields('TRI-3D', struct('size', grid_sizes), noise_std);

    disp('Testing Scalar Field Tricubic');
    test_fields('TRI-LPL', struct('size', grid_sizes), noise_std);

%     disp('Testing Divergence Free Tricubic');
%     test_fields('TRI-DF', struct('size', grid_sizes), noise_std);

    disp('Testing 3D BSpline');
    test_fields('SPL-3D', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);

    disp('Testing Laplacian BSpline');
    test_fields('SPL-LPL', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);
end

output_files = dir('data/fields/*.mat');

close all;

Nf = length(output_files);
cmap = cbrewer('qual', 'Set1', Nf);

mae = zeros(Nf, length(grid_sizes));
nmae = zeros(Nf, length(grid_sizes));
rmse = zeros(Nf, length(grid_sizes));
nrmse = zeros(Nf, length(grid_sizes));
r2 = zeros(Nf, length(grid_sizes));
mean_div = zeros(Nf, length(grid_sizes));
mean_curl = zeros(Nf, length(grid_sizes));
model_names = {};
for i=1:Nf
    filename = fullfile(output_files(i).folder, output_files(i).name);
    
    temp = strsplit(output_files(i).name, '_');
    model_name = temp{1};
    model_names{i} = model_name;
    noise_std = temp{2} * 1e-6;
    
    results = importdata(filename);
    mae(i,:) = mean([results.mae]);
    nmae(i,:) = mean([results.nmae]);
    rmse(i,:) = mean([results.rmse]);
    nrmse(i,:) = mean([results.nrmse]);
    r2(i,:) = mean([results.r2]);
end

%% nmae
fh_mae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);
colormap(cmap);
[~, idx] = sort(nmae(:,1), 1);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = cmap(i,:);
end
ax = fh_mae.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
legend(model_names(idx));

%set(fh_mae, 'PaperUnits', 'inches');
%set(fh_mae, 'PaperSize', [3.45/2, 2.1]);

export_fig(fh_mae, 'figures/interp_field_nmae.pdf');

%% r2

fh_r2 = figure('Name', 'Mean R2', 'units', 'inch', ...
     'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 8);
colormap(cmap);
% we sort the by the r2 in the lowest grid resolution
%[~, idx] = sort(r2(:,1), 1, 'descend');
b = bar([results.grid_size], r2(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = cmap(i,:);
end
ax = fh_r2.CurrentAxes;
xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, '$R^2$ ', 'Interpreter', 'latex');
ylim(ax, [min(r2(:))-0.1, 1])
legend(model_names(idx));

export_fig(fh_r2, 'figures/interp_field_r2.pdf');

%% Table generation
input.data = [1000*mae(:,3)'; 100*nmae(:,3)'; 1000*rmse(:,3)'; 100*nrmse(:,3)'; r2(:,3)'];
input.tableRowLabels = {'MAE (mT)', 'N-MAE (\%)', 'RMSE (mT)', 'N-RMSE (\%)', '$R^2$'};
input.dataFormatMode = 'row';
input.dataFormat = {'%.1f', 4, '%.3f', 1};
input.tableBorders = 0;
input.tableColLabels = model_names;
input.tableCaption = 'Field interpolation performance results with $n_g = 5$';
input.tableLabel = 'tab:interp_performance';
table_mae = latexTable(input);

%% Color Generation
colors = containers.Map(model_names(idx), num2cell(cmap', [1,3]));
save('data/colors', 'colors');
save('data/idx', 'idx');