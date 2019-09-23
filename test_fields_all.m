clear all;

RECOMPUTE = 1;

grid_sizes = {3,4,5,6};
noise_std = 0;

if RECOMPUTE ~= 0
%     disp('Testing RBF 3D');
%     test_fields('RBF-3D', struct('size', grid_sizes, 'eps', {1., 10., 22., 50.}), noise_std);

%     disp('Testing RBF Multiquadric 3D');
%     test_fields('RBF-MQ-3D', struct('size', grid_sizes, 'eps', {10.4525, 14.6804, 24.3048,  24.6842}), noise_std);
% 
%     disp('Testing RBF Div-free');
%     test_fields('RBF-DF', struct('size', grid_sizes, 'eps', {35, 35, 50, 50}), noise_std);
    disp('Testing RBF Multiquadric Div-free');
    test_fields('RBF-MQ-DF', struct('size', grid_sizes, 'eps', {35, 35, 50, 50}), noise_std);
% 
%     disp('Testing Tricubic 3D');
%     test_fields('TRI-3D', struct('size', grid_sizes), noise_std);
% 
%     disp('Testing Scalar Field Tricubic');
%     test_fields('TRI-LPL', struct('size', grid_sizes), noise_std);

%     disp('Testing Divergence Free Tricubic');
%     test_fields('TRI-DF', struct('size', grid_sizes), noise_std);
% 
%     disp('Testing 3D BSpline');
%     test_fields('SPL-3D', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);

%     disp('Testing Laplacian BSpline');
%     test_fields('SPL-LPL', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);
end

output_files = dir('data/fields/*.mat');

close all;

cmap = [0.105882352941176 0.619607843137255 0.466666666666667;
    0.850980392156863 0.372549019607843 0.00784313725490196;
    0.458823529411765 0.43921568627451 0.701960784313725;
    0.905882352941176 0.16078431372549 0.541176470588235;
    0.4 0.650980392156863 0.117647058823529;
    0.901960784313726 0.670588235294118 0.00784313725490196];

Nf = length(output_files);
mae = zeros(Nf, length(grid_sizes));
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
    r2(i,:) = mean([results.r2]);
end
fh_mae = figure('Name', 'Mean MAE');
colormap(cmap);
bar([results.grid_size], mae', 'grouped');
ax = fh_mae.CurrentAxes;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'Mean Absolute Error (mT)', 'Interpreter', 'latex');
legend(model_names);

fh_r2 = figure('Name', 'Mean R2');
colormap(cmap);
bar([results.grid_size], r2', 'grouped');
ax = fh_r2.CurrentAxes;
xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'R2', 'Interpreter', 'latex');
ylim(ax, [min(r2(:))-0.1, 1])
legend(model_names);