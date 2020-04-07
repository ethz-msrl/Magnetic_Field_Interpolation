function [ output ] = evaluate_divfree_multiquadric_rbf_gradient( position, nodes, eps, C )
%EVALUTE_RBF Evaluates the gradient of a 3D divergence-free radial basis
%function with a multiquadric kernel
%
%   Evaluates the function:
%   del(y(position)) = sum_j ( Psi(position - node_j) * C_j)
%
%   where Psi is the 3D matrix-valued divergence-free RBF with a
%   multiquadric
%   kernel
%
%   position: 3D vector at which to evaluate the RBF
%   nodes: Nx3 matrix of node positions
%   eps: scalar value of the decay param of the RBF
%   C: Nx3 coefficient matrix
%
%   output 3x3 symmetric gradient matrix
%
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

x = position(1) - nodes(:,1)';
y = position(2) - nodes(:,2)';
z = position(3) - nodes(:,3)';
x2 = x.^2;
y2 = y.^2;
z2 = z.^2;

x3 = x.^3;
y3 = y.^3;
z3 = z.^3;

K = sqrt(1 + eps * (x2 + y2 + z2));
K3 = K.^3;
K5 = K.^5;

del_psi_11 = [
    (3*eps^3 * x .* ( y2 + z2)) ./ K5 - (2*eps^2*x) ./ K3;
    (3*eps^3*(y3 + y .* z2)) ./ K5 - (4*eps^2*y) ./ K3; 
    (3*eps^3* (z3 + y2 .*z)) ./ K5  - (4*eps^2*z) ./ K3;
    ];

del_psi_12 = [
    (eps^2*y) ./ K3 - (3*eps^3*x2 .*y ) ./ K5;
    (eps^2*x) ./ K3 - (3*eps^3*x.*y2) ./ K5;
    -(3*eps^3*x.*y.*z) ./ K5
    ];

del_psi_13 = [
    (eps^2*z) ./ K3 - (3*eps^3* x2 .*z) ./ K5;
    -(3*eps^3*x.*y.*z) ./ K5;
    (eps^2*x) ./ K3 - (3*eps^3*x.*z2) ./ K5;
    ];

del_psi_22 = [
    (3*eps^3* (x3 + x.*z2)) ./ K5 - (4*eps^2*x) ./ K3;
    (3*eps^3*(x2.*y + y.*z2)) ./ K5 - (2*eps^2*y) ./ K3;
    (3*eps^3*(z3 + x2.*z)) ./ K5 - (4*eps^2*z) ./ K3;
    ];

del_psi_23 = [
    -(3*eps^3*x.*y.*z) ./ K5;
    (eps^2*z) ./ K3 - (3*eps^3*y2.*z) ./ K5;
    (eps^2*y) ./ K3 - (3*eps^3*y.*z2) ./ K5;
    ];

del_psi_33 = [
    (3*eps^3*(x3 + x.*y2)) ./ K5 - (4*eps^2*x) ./ K3;
    (3*eps^3*(y3 + x2.*y)) ./ K5 - (4*eps^2*y) ./ K3;
    (3*eps^3*(x2.*z + y2.*z)) ./ K5 - (2*eps^2*z) ./ K3;
    ];

output = [
    del_psi_11 * C(:,1) + del_psi_12 * C(:,2) + del_psi_13 * C(:,3), ...
    del_psi_12 * C(:,1) + del_psi_22 * C(:,2) + del_psi_23 * C(:,3), ...
    del_psi_13 * C(:,1) + del_psi_23 * C(:,2) + del_psi_33 * C(:,3)];

end

