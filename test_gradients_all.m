clear all;

RECOMPUTE = 0;

grid_sizes = {3,4,5,6};

noise_std = 0;

if RECOMPUTE ~= 0
    disp('Testing RBF 3D');
    load('data/best_eps/RBF-3D', 'best_eps');
    test_gradients('RBF-G-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);
   
    disp('Testing RBF Multiquadric 3D');
    load('data/best_eps/RBF-MQ-3D', 'best_eps');
    test_gradients('RBF-MQ-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing RBF Div-free');
    load('data/best_eps/RBF-DF', 'best_eps');
    test_gradients('RBF-G-DF', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing RBF Multiquadric Div-free');
    load('data/best_eps/RBF-MQ-DF');
    test_gradients('RBF-MQ-DF', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing Tricubic 3D');
    test_gradients('TRI-3D', struct('size', grid_sizes), noise_std);

    disp('Testing Scalar Field Tricubic');
    test_gradients('TRI-LPL', struct('size', grid_sizes), noise_std);
    
    disp('Testing 3D BSpline');
    test_gradients('SPL-3D', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);

    disp('Testing Laplacian BSpline');
    test_gradients('SPL-LPL', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);
end

output_files = dir('data/gradients/*.mat');

close all;

load('data/colors');
load('data/idx');

Nf = length(output_files);
%cmap = cbrewer('qual', 'Set1', Nf);

nmae = zeros(Nf, length(grid_sizes));
mae = zeros(Nf, length(grid_sizes));
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
%     bar(fh_mae.CurrentAxes, [results.grid_size], mean([results.mae],1), ...
%         'DisplayName', model_name, 'FaceColor', cmap(i,:));
    r2(i,:) = mean([results.r2]);
    mean_div(i,:) = [results.mean_div];
    mean_curl(i,:) = sqrt(sum(reshape([results.mean_curl],3,4)'.^2,2));
end
fh_nmae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);
%colormap(cmap);
ax = gca;
%ax.ColorOrder = cmap;
%[~, idx] = sort(nmae(:,1), 1);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    %b(i).FaceColor = cmap(idx(i),:);
    b(i).FaceColor = colors(model_names{idx(i)});
end
ax = fh_nmae.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
legend(model_names(idx));

%set(fh_nmae, 'PaperUnits', 'inches');
%set(fh_nmae, 'PaperSize', [3.45/2, 2.1]);

export_fig(fh_nmae, 'figures/interp_gradient_nmae.pdf');

fh_r2 = figure('Name', 'Mean R2', 'units', 'inch', ...
     'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);

% we sort the by the r2 in the lowest grid resolution
%[~, idx] = sort(r2(:,1), 1, 'descend');
% hold on;
% for i=1:Nf
%     b = bar([results.grid_size], r2(idx(i),:)');
%     b.FaceColor = cmap(idx(i),:);
% end

b = bar([results.grid_size], r2(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    %b(i).FaceColor = cmap(idx(i),:);
    b(i).FaceColor = colors(model_names{idx(i)});
end

ax = fh_r2.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, '$R^2$ ', 'Interpreter', 'latex');
ylim(ax, [min(r2(:))-0.1, 1])
legend(model_names(idx));

export_fig(fh_r2, 'figures/interp_gradient_r2.pdf');

%% div
idx_r = find(abs(mean_div(:,1)) > 1e-12);

fh_md = figure('Name', 'Mean Divergence', 'units', 'inch', ...
     'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);
%colormap(cmap);
b = bar([results.grid_size], 1000*mean_div(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax = fh_md.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'Mean Divergence (mT/m)', 'Interpreter', 'latex');
legend(model_names(idx_r));

export_fig(fh_md, 'figures/interp_divergence.pdf');

%% curl
idx_r = find(abs(mean_curl(:,1)) > 1e-12);

fh_mc = figure('Name', 'Mean Curl Magnitude', 'units', 'inch', ...
     'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);
%colormap(cmap);
b = bar([results.grid_size], 1000*mean_curl(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax = fh_mc.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');

ylabel(ax, 'Mean Curl Magnitude (mT/m)', 'Interpreter', 'latex');
legend(model_names(idx_r));

export_fig(fh_mc, 'figures/interp_curl.pdf');

%% Table generation
input.data = [1000*mae(:,3)'; 100*nmae(:,3)'; 1000*rmse(:,3)'; 100*nrmse(:,3)'; ...
    r2(:,3)'; 1000*mean_div(:,3)'; 1000*mean_curl(:,3)'];
input.tableRowLabels = {'MAE (mT/m)', 'N-MAE (\%)', 'RMSE (mT/m)', 'N-RMSE (\%)', ...
    '$R^2$', '$\nabla \cdot \mathbf{b}$ ($\SI{}{\milli\tesla/\meter}$)', '$||\nabla \times \mathbf{b}||$ ($\SI{}{\milli\tesla/\meter}$)'};
input.dataFormatMode = 'row';
input.dataFormat = {'%.1f', 4, '%.3f', 1, '%.1f', 2 };
input.tableBorders = 0;
input.tableColLabels = model_names;
input.tableCaption = 'Gradient interpolation performance results with $n_g = 5$';
input.tableLabel = 'tab:interp_performance';
table_mae = latexTable(input);
