function [ output, gradient ] = evaluate_gaussian_rbf( position, nodes, eps, C )
%EVALUTE_RBF Evaluates the value of a 3D radial basis function
%with a Gaussian kernel
%
%   Evaluates the function:
%   y(position) = sum_j ( Psi(position - node_j) * C_j)
%
%   where Psi is the Gaussian kernel
%
%   position: 3D vector at which to evaluate the RBF
%   nodes: Nx3 matrix of node positions
%   eps: scalar value of the decay param of the RBF
%   C: Nx3 coefficient matrix
%
%   Returns:
%   output: 3x1 interpolated value
%   gradient: 3x3 symmetric gradient of interpolant
%
%   Copyright 2020, Samuel Charreyron

    xd = nodes(:,1);
    yd = nodes(:,2);
    zd = nodes(:,3);
    K = exp(-eps * ((position(1) - xd).^2 + (position(2) - yd).^2 + ...
        (position(3) - zd).^2));

    output = [
        dot(K, C(:,1));
        dot(K, C(:,2));
        dot(K, C(:,3));
        ];
    
    pdiff = repmat(position, [size(nodes, 1), 1]) - nodes;
    gradient = -2*eps*[ ...
        (K .* C(:,1))' * pdiff ; ...
        (K .* C(:,2))' * pdiff; ...
        (K .* C(:,3))' * pdiff];
        

end