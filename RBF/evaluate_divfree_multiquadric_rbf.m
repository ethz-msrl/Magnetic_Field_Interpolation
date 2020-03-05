function [ output ] = evaluate_divfree_multiquadric_rbf( position, nodes, eps, C )
%EVALUTE_DIVFREE_RBF Evaluates the value of a 3D divergence-free radial basis
%function with a multiquadric kernel
%
%   Evaluates the function:
%   y(position) = sum_j ( Psi(position - node_j) * C_j)
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
%   output: 3x1 vector interpolated value
%
%   Copyright 2020, Samuel Charreyron


    xd = nodes(:,1);
    yd = nodes(:,2);
    zd = nodes(:,3);
    
    K = sqrt(1 + eps * ((position(1) - xd).^2 + (position(2) - yd).^2 + ...
        (position(3) - zd).^2));
    
    K3 = K.^3;

    psi_11 = 2*eps ./ K - eps^2 * ((position(2) - yd).^2 + (position(3) - zd).^2) ./ K3 ;
    psi_22 = 2*eps ./ K - eps^2 * ((position(1) - xd).^2 + (position(3) - zd).^2) ./ K3 ;
    psi_33 = 2*eps ./ K - eps^2 * ((position(1) - xd).^2 + (position(2) - yd).^2) ./ K3 ;

    psi_12 = eps^2 * (position(1) - xd) .* (position(2) - yd) ./ K3;
    psi_13 = eps^2 * (position(1) - xd) .* (position(3) - zd) ./ K3;
    psi_23 = eps^2 * (position(2) - yd) .* (position(3) - zd) ./ K3;
        
    output = [
        dot(psi_11, C(:,1)) + dot(psi_12, C(:,2)) + dot(psi_13, C(:,3));
        dot(psi_12, C(:,1)) + dot(psi_22, C(:,2)) + dot(psi_23, C(:,3));
        dot(psi_13, C(:,1)) + dot(psi_23, C(:,2)) + dot(psi_33, C(:,3));
        ];

end