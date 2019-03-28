% function [ output_args ] = load_mag_data( filepath )
%LOAD_MAG_DATA Summary of this function goes here
%   Detailed explanation goes here

% load currents file
currents = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/currents_3787.h5', '/currents')';

% load positions
xg = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/positions.h5', '/xg');
yg = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/positions.h5', '/yg');
zg = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/positions.h5', '/zg');

xg = permute(xg, [3, 2, 1]);
yg = permute(yg, [3, 2, 1]);
zg = permute(zg, [3, 2, 1]);

xd = xg(:);
yd = yg(:);
zd = zg(:);

nodes = [xd, yd, zd];

fields = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/v/0000.h5', '/fields');

fields = permute(fields, [4, 3, 2, 1]);

% load positions
xg_ev = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/positions.h5', '/xg');
yg_ev = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/positions.h5', '/yg');
zg_ev = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/positions.h5', '/zg');

%xd_ev = xg_ev(:);
%yd_ev = yg_ev(:);
%zd_ev = zg_ev(:);

positions_ev = cat(4, xg_ev, yg_ev, zg_ev);

fields_ev = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_8_h5/v/0000.h5', '/fields');
fields_ev = permute(fields_ev, [4, 3, 2, 1]);

% end

