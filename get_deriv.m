function [ dvg_dy, dvg_dx, dvg_dz ] = get_deriv( vg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dvg_dx_fd = diff(vg, 1, 1);
dvg_dx = cat(1, dvg_dx_fd, vg(end,:,:,:) - vg(end-1,:,:,:));

dvg_dy_fd = diff(vg, 1, 2);
dvg_dy = cat(2, dvg_dy_fd, vg(:,end,:,:) - vg(1,end-1,:,:));

dvg_dz_fd = diff(vg, 1, 3);
dvg_dz = cat(3, dvg_dz_fd, vg(:,:,end,:) - vg(:,:,end-1,:));

end

