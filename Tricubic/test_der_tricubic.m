%% This script is to test that the derivatives of the tribubic interpolator
% are calculated correctly. We compute the derivatives using an analytical
% expression and calculate them with numerical derivatives and compare
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

A = rand(4,4,4);

eps = 1e-8;

x = [0.1, 0.2, 0.3];
f0 = tricubic(A, x(1) , x(2), x(3));

df_dx = (tricubic(A, x(1) + eps, x(2), x(3)) - f0) / eps;
df_dy = (tricubic(A, x(1), x(2) + eps, x(3)) - f0) / eps;
df_dz = (tricubic(A, x(1), x(2), x(3) + eps) - f0) / eps;


df_dv = tricubic_grad(A, x(1), x(2), x(3));
df_dv_ = [df_dx; df_dy; df_dz] ;

assert(all(isalmost(df_dv, df_dv_, 1e-4)));

df0 = tricubic_grad(A, x(1), x(2), x(3));
df2_dx2 = (tricubic_grad(A, x(1) + eps, x(2), x(3)) - df0 ) /eps;
df2_dy2 = (tricubic_grad(A, x(1), x(2) + eps, x(3)) - df0 ) /eps;
df2_dz2 = (tricubic_grad(A, x(1), x(2), x(3) + eps) - df0 ) /eps;
df2_dv2_ = [df2_dx2, df2_dy2, df2_dz2];

df2_dv2 = tricubic_hess(A, x(1), x(2), x(3));

assert(all(isalmost(df2_dv2, df2_dv2_, 1e-4), 'all'));

%% Testing derivatives of Tricubic 
grid_size = 5;
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

fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', 0)), '/fields');
            fields = permute(fields, [4, 3, 2, 1]);

load('tricubic_scalar_field_M.mat');
model = TricubicScalarFieldInterpolator(nodes, fields, M);

pev = [0,0,0];
b = model.getFieldAtPosition(pev);
b_ = model.getFieldAtPostionNumeric(pev);

assert(all(isalmost(b, b_, 1e-4)));

g_ = model.getGradientAtPositionNumeric(pev);
g = model.getGradientAtPosition(pev);
assert(all(isalmost(g, g_, 1e-4), 'all'));

