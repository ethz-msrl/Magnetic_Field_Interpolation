function [ curl ] = tricubic_curl( A, x, y, z )
%TRICUBIC_CURL Calculates the curl of a 3D tricubic function
%   Arguments:
%               A: as in tricubic, it should be 4x4x4 matrix
%               x: should be a scalar
%               y: should be a scalar
%               z: should be a scalar
%   Returns:
%               curl: a vector containing the curl
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

xv = kron(kron([3*x^2, 2*x, 1, 0], [y^3, y^2, y, 1]), [z^3, z^2, z, 1])';
yv = kron(kron([x^3, x^2, x, 1], [3*y^2, 2*y, 1, 0]), [z^2, z^2, z, 1])';
zv = kron(kron([x^3, x^2, x, 1], [y^3, y^2, y, 1]), [3*z^2, 2*z, 1, 0])';

% av is the vector arrangement of A which is 4x4x4
% this ensures that the ordering of the parameters of av matches xv
avx = flip(reshape(permute(A(:,:,:,1), [3,2,1]), [], 1));
avy = flip(reshape(permute(A(:,:,:,2), [3,2,1]), [], 1));
avz = flip(reshape(permute(A(:,:,:,3), [3,2,1]), [], 1));

curl = [
    sum(avz .* yv - avy .* zv );
    sum(avx .* zv - avz .* xv );
    sum(avy .* xv - avx .* yv );
    ];
end

