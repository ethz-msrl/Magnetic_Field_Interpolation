%% This script outputs interpolation data using the TRI-3D method to be
% checked with the C++ implementation

clear variables;
load_settings;

grid_size = 5;        
nodes_dataset = sprintf(options.base_dataset, 5);

% the are the nodes that are used by the interpolant
nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
xg = h5read(nodes_pos_fn, '/xg');
yg = h5read(nodes_pos_fn, '/yg');
zg = h5read(nodes_pos_fn, '/zg');

% the positions are in the z,y,x format because of deep-fluids
%xg = permute(xg, [3, 2, 1]);
%yg = permute(yg, [3, 2, 1]);
%zg = permute(zg, [3, 2, 1]);

nodes = cat(4, xg, yg, zg);

% these are the positions at which the interpolation is evaluated
% usually it's a fine grid with Ng=16
eval_pos_fn = fullfile(options.eval_dataset,'/positions.h5');

xg_ev = h5read(eval_pos_fn, '/xg');
yg_ev = h5read(eval_pos_fn, '/yg');
zg_ev = h5read(eval_pos_fn, '/zg');

%xg_ev = permute(xg_ev, [3, 2, 1]);
%yg_ev = permute(yg_ev, [3, 2, 1]);
%zg_ev = permute(zg_ev, [3, 2, 1]);

% the positions are in the z,y,x format because of deep-fluids
positions_ev = cat(4, xg_ev, yg_ev, zg_ev);

fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', 0)), '/fields');
%fields = permute(fields, [4, 3, 2, 1]);
fields = permute(fields, [2,3,4,1]);

% these are the field values on the evaluation grid
% note that they are obtained from synthetic data generated
% with the MPEM so they are naturally curl and divergence free
fields_ev = h5read(fullfile(options.eval_dataset, 'v', sprintf('%04d.h5', 0)), '/fields');
%fields_ev = permute(fields_ev, [4, 3, 2, 1]);
fields_ev = permute(fields_ev, [2,3,4,1]);
load('tricubic_3D_M.mat');
model = Tricubic3DInterpolator(nodes, fields, M);
fields_interp = model.getFieldsAtPositions(positions_ev);

ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);

writematrix(reshape(permute(nodes, [3,2,1,4]), [], 3), ...
    'interp_Ng16_0_nodes.txt', 'Delimiter', ' ');
writematrix(reshape(permute(fields, [3,2,1,4]), [], 3), ...
    'interp_Ng16_0_fields.txt', 'Delimiter', ' ');
writematrix(reshape(permute(positions_ev, [3,2,1,4]), [], 3), ...
    'interp_Ng16_0_pos_ev.txt', 'Delimiter', ' ');
writematrix(reshape(permute(fields_interp, [3,2,1,4]), [], 3), ...
    'interp_Ng16_0_tri3d_ev.txt', 'Delimiter', ' ');