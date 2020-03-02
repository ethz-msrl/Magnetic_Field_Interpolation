clear all;

RECOMPUTE = 1;

grid_sizes = {3,4,5,6};

save('data/grid_sizes', 'grid_sizes');

noise_std = 0;

if RECOMPUTE == 1
    disp('Testing RBF 3D');
    load('data/best_eps/RBF-G-3D', 'best_eps');
    test_gradients('RBF-G-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);

    disp('Testing RBF Multiquadric 3D');
    load('data/best_eps/RBF-MQ-3D', 'best_eps');
    test_gradients('RBF-MQ-3D', struct('size', grid_sizes, 'eps', num2cell(best_eps(1:length(grid_sizes)))), noise_std);
    
    disp('Testing RBF Div-free');
    load('data/best_eps/RBF-G-DF', 'best_eps');
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

load_gradient_data;

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
