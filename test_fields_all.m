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
%     
    % Using a degree equal to the grid size appears to work best for some
    % reason
    disp('Testing 3D BSpline');
    test_fields('SPL-3D', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);

    disp('Testing Laplacian BSpline');
    test_fields('SPL-LPL', struct('size', grid_sizes, 'degree', {3, 4, 5, 6}), noise_std);
end

%% Color Generation
% colors = containers.Map(model_names(idx), num2cell(cmap', [1,3]));
% save('data/colors', 'colors');
% save('data/idx', 'idx');

load_field_data;

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

