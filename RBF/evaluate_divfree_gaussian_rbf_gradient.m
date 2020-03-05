function [ output ] = evaluate_divfree_gaussian_rbf_gradient( position, nodes, eps, C )
%EVALUTE_DIVFREE_GAUSSIAN_RBF_GRADIENT Evaluates the gradient of a 3D divergence-free radial basis
%function with a Gaussian kernel
%
%   Evaluates the function:
%   del(y(position)) = sum_j ( Psi(position - node_j) * C_j)
%
%   where Psi is the 3D matrix-valued divergence-free RBF with a Gaussian
%   kernl
%
%   position: 3D vector at which to evaluate the RBF
%   nodes: Nx3 matrix of node positions
%   eps: scalar value of the decay param of the RBF
%   C: Nx3 coefficient matrix
%
%   output 3x3 symmetric gradient matrix
%
%   Copyright 2020, Samuel Charreyron

x = position(1) - nodes(:,1)';
y = position(2) - nodes(:,2)';
z = position(3) - nodes(:,3)';

K = exp(-eps * (x.^2 + y.^2 + z.^2));

del_psi_11 = [
    (-8*eps^2 * x + 8*eps^3 * (x .* (y.^2) + x .* (z.^2))) .* K;
    (-16*eps^2 * y + 8*eps^3 * (y.^3 + y .* (z.^2))) .* K; 
    (-16*eps^2 * z + 8*eps^3 * (z.^3 + (y.^2) .* z)) .* K;
    ];

del_psi_12 = [
    (-8*eps^3 * (x.^2) .* y + 4*eps^2 * y) .* K;
    (-8*eps^3 * x .* (y.^2) + 4*eps^2 * x) .* K; 
    (-8*eps^3 * x .* y .* z) .* K;
    ];

del_psi_13 = [
    (-8*eps^3 * (x.^2) .* z + 4*eps^2 * z) .* K;
    (-8*eps^3 * x .* y .* z) .* K; 
    (-8*eps^3 * x .* (z.^2) + 4 * eps^2 * x) .* K;
    ];

del_psi_22 = [
    (-16*eps^2 * x + 8*eps^3 * (x.^3 + x.*(z.^2))) .* K;
    (-8*eps^2 * y + 8*eps^3 * ((x.^2) .* y + y .* (z.^2))) .* K; 
    (-16*eps^2 * z + 8*eps^3 * (z.^3 + (x.^2) .* z)) .* K;
    ];

del_psi_23 = [
    (-8*eps^3 * x .* y .* z) .* K;
    (-8*eps^3 * (y.^2) .* z + 4*eps^2 * z) .* K; 
    (-8*eps^3 * y .* (z.^2) + 4 * eps^2 * y) .* K;
    ];

del_psi_33 = [
    (-16*eps^2 * x + 8*eps^3 * (x.^3 + x .* (y.^2))) .* K;
    (-16*eps^2 * y + 8*eps^3 * (y.^3 + (x.^2) .* y)) .* K; 
    (-8*eps^2 * z + 8*eps^3 * ((x.^2 .* z) + (y.^2) .* z)) .* K;
    ];

output = [
    del_psi_11 * C(:,1) + del_psi_12 * C(:,2) + del_psi_13 * C(:,3), ...
    del_psi_12 * C(:,1) + del_psi_22 * C(:,2) + del_psi_23 * C(:,3), ...
    del_psi_13 * C(:,1) + del_psi_23 * C(:,2) + del_psi_33 * C(:,3)];

end

