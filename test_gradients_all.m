% This script tests all interpolation methods on their ability to
% interpolate magnetic field gradients.
%
% For the RBF values optimal shape parameters that are obtained from the
% test_optimal_eps experiment are used. The positions are normalized using
% normalize_positions_maxmin to be between 0 and 1.
%
% For the B-Spline methods, we use a degree that matches the grid size.
% That seemed to work best for us.
%
% The data is saved to data/gradients with a file corresponding to the method
% name.
% This script also generates latex code with the results from the Ng=5 grid
% that can be pasted into the paper.
%
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

clear variables;
load_settings;

if options.recompute ~= 0
    disp('Testing RBF Gaussian 3D');
    load('data/best_eps/RBF-G-3D', 'best_eps');
    test_gradients('RBF-G-3D', struct('size', options.grid_sizes, 'eps', num2cell(best_eps(1:length(options.grid_sizes)))), options);

    disp('Testing RBF Multiquadric 3D');
    load('data/best_eps/RBF-MQ-3D', 'best_eps');
    test_gradients('RBF-MQ-3D', struct('size', options.grid_sizes, 'eps', num2cell(best_eps(1:length(options.grid_sizes)))), options);
    
    disp('Testing RBF Gaussian Div-free');
    load('data/best_eps/RBF-G-DF', 'best_eps');
    test_gradients('RBF-G-DF', struct('size', options.grid_sizes, 'eps', num2cell(best_eps(1:length(options.grid_sizes)))), options);
    
    disp('Testing RBF Multiquadric Div-free');
    load('data/best_eps/RBF-MQ-DF');
    test_gradients('RBF-MQ-DF', struct('size', options.grid_sizes, 'eps', num2cell(best_eps(1:length(options.grid_sizes)))), options);
    
    disp('Testing Tricubic 3D');
    test_gradients('TRI-3D', struct('size', options.grid_sizes), options);
    
    disp('Testing Scalar Field Tricubic');
    test_gradients('TRI-LPL', struct('size', options.grid_sizes), options);
    
    disp('Testing 3D BSpline');
    test_gradients('SPL-3D', struct('size', options.grid_sizes, 'degree', options.bspline_degrees), options);
    
    disp('Testing Laplacian BSpline');
    test_gradients('SPL-LPL', struct('size', options.grid_sizes, 'degree', options.bspline_degrees), options);
end

load_gradient_data;

%% Table generation (for paper)
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
