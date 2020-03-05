function  plot_field_streamslice( model_name, grid_size, current_index, noise_std, eps)
% PLOT_FIELD_STREAMSLICE plots the streamslice
%   model_name: a string defining the interpolation model type
%   grid: a struct containing both the grid size and params that are needed
%   by the model
%   noise_std: the standard deviation or Normally distributed noise to add
%   to the data (set to 0 for no noise)
%
%   Copyright 2020, Samuel Charreyron

i = current_index;

nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grid_size);
EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

% load positions
nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
xg = h5read(nodes_pos_fn, '/xg');
yg = h5read(nodes_pos_fn, '/yg');
zg = h5read(nodes_pos_fn, '/zg');

xg = permute(xg, [3, 2, 1]);
yg = permute(yg, [3, 2, 1]);
zg = permute(zg, [3, 2, 1]);

nodes = cat(4, xg, yg, zg);

% load positions
eval_pos_fn = fullfile(EVAL_DATASET,'/positions.h5');

xg_ev = h5read(eval_pos_fn, '/xg');
yg_ev = h5read(eval_pos_fn, '/yg');
zg_ev = h5read(eval_pos_fn, '/zg');

xg_ev = permute(xg_ev, [3, 2, 1]);
yg_ev = permute(yg_ev, [3, 2, 1]);
zg_ev = permute(zg_ev, [3, 2, 1]);

positions_ev = cat(4, xg_ev, yg_ev, zg_ev);

fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
fields = permute(fields, [4, 3, 2, 1]);
% adding random noise
fields = fields + noise_std * randn(size(fields));

% fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
% fields_ev = permute(fields_ev, [4, 3, 2, 1]);

if strcmp(model_name, 'RBF-G-3D')
    model = RBF3DInterpolator(nodes, fields, eps);
elseif strcmp(model_name, 'RBF-MQ-3D')
    model = RBF3DMultiquadricInterpolator(nodes, fields, eps);
elseif strcmp(model_name, 'RBF-G-DF')
    model = RBFDivFreeInterpolator(nodes, fields, eps);
elseif strcmp(model_name, 'RBF-MQ-DF')
    model = RBFDivFreeMultiquadricInterpolator(nodes, fields, eps);
elseif strcmp(model_name, 'SPL-3D')
    degree = 4;
    model = BSpline3DInterpolator(nodes, fields, degree);
elseif strcmp(model_name, 'SPL-LPL')
    degree = 4;
    model = BSplineLaplacianInterpolator(nodes, fields, degree);
elseif strcmp(model_name, 'TRI-3D')
    load('tricubic_3D_M.mat');
    model = Tricubic3DInterpolator(nodes, fields, M);
elseif strcmp(model_name, 'TRI-LPL')
    load('tricubic_scalar_field_M.mat');
    model = TricubicScalarFieldInterpolator(nodes, fields, M);
elseif strcmp(model_name, 'TRI-DF')
    load('tricubic_divfree_M.mat');
    model = TricubicDivFreeInterpolator(nodes, fields, M);
else
    error('invalid model name');
end

[xge, yge] = meshgrid(linspace(-0.1, 0.1, 16), linspace(-0.1, 0.1, 16));
nzg = floor(size(zg, 1)/2)+1;
positions_ev = cat(4, xge, yge, zg(nzg,1,1) * ones(size(xge)));

fields_interp = model.getFieldsAtPositions(positions_ev);
ss = streamslice(1000*xge, 1000*yge, fields_interp(:,:,:,1), fields_interp(:,:,:,2)); 
set(ss, 'Color', [0, 0.4470, 0.7410]);
set(ss, 'LineWidth', 1.0);
hold on;

quiver(1000*xg(nzg,:,:), 1000*yg(nzg,:,:), fields(nzg, :, :, 1), fields(nzg, :, :, 2), 'r', ...
'AutoScaleFactor', 0.3, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.5);

xlabel('x (mm)');
ylabel('y (mm)');
axis('equal');
%axis 'tight';

hold off;

end

