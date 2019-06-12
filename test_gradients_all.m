clear all;

RECOMPUTE = 1;

grid_sizes = {3,4,5,6};

if RECOMPUTE ~= 0
%     disp('Testing RBF 3D');
%     test_gradients('RBF3D', struct('size', grid_sizes, 'eps', {1., 10., 22., 50.}), 0);

%     disp('Testing RBF Div-free');
%     test_gradients('RBFD', struct('size', grid_sizes, 'eps', {35, 35, 50, 50}), 0);
% 
    disp('Testing Tricubic 3D');
    test_gradients('TRI3D', struct('size', grid_sizes), 0);

%     disp('Testing Scalar Field Tricubic');
%     test_gradients('TRILPL', struct('size', grid_sizes), 0);
% 
%     disp('Testing 3D BSpline');
%     test_gradients('SPL3D', struct('size', grid_sizes, 'degree', {3, 4, 4, 4}), 0);
% 
%     disp('Testing Laplacian BSpline');
%     test_gradients('SPLLPL', struct('size', grid_sizes, 'degree', {3, 4, 4, 4}), 0);
end

output_files = dir('data/gradients');

for i=1:length(output_files)
    filename = output_files(i);
end