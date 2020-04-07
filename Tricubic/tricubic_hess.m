function [ hess ] = tricubic_hess( A, x, y, z )
%TRICUBICHESS Calculates the Hessian of a tricubic function
%   Arguments:
%               A: as in tricubic, it should be 4x4x4 matrix
%               x: should be a scalar
%               y: should be a scalar
%               z: should be a scalar
%   Returns:
%               hess: a 3x3 symmetric Hessian matrix
%
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

xv = [x^3, x^2, x, 1]';
xv1 = [3*x^2, 2*x, 1, 0]';
xv2 = [6*x, 2, 0, 0]';

yv = [y^3, y^2, y, 1]';
yv1 = [3*y^2, 2*y, 1, 0]';
yv2 = [6*y, 2, 0, 0]';

zv = [z^3, z^2, z, 1]';
zv1 = [3*z^2, 2*z, 1, 0]';
zv2 = [6*z, 2, 0, 0]';

% av is the vector arrangement of A which is 4x4x4
% this ensures that the ordering of the parameters of av matches xv
av = flip(reshape(permute(A, [3,2,1]), [], 1));

dgx = [
    kron(kron(xv2, yv), zv)';
    kron(kron(xv1, yv1), zv)'; 
    kron(kron(xv1, yv), zv1)'; 
    ]';
dgy = [
    kron(kron(xv1, yv1), zv)';
    kron(kron(xv, yv2), zv)';
    kron(kron(xv, yv1), zv1)';
    ]';
dgz = [
    kron(kron(xv1, yv), zv1)';
    kron(kron(xv, yv1), zv1)';
    kron(kron(xv, yv), zv2)';
    ]';

H = zeros(64, 3, 3);
H(:,1,:) = dgx;
H(:,2,:) = dgy;
H(:,3,:) = dgz;

hess = squeeze(sum(repmat(av, [1,3,3]) .* H, 1));

end
