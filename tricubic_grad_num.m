function [ grad ] = tricubic_grad_num( A, x, y, z )
%TRICUBICGRADNUM Calculates the gradient of a tricubic function
%   Computes the gradient of tricubic numerically.
xv = kron(kron([3*x^2, 2*x, 1, 0], [y^3, y^2, y, 1]), [z^3, z^2, z, 1])';
yv = kron(kron([x^3, x^2, x, 1], [3*y^2, 2*y, 1, 0]), [z^2, z^2, z, 1])';
zv = kron(kron([x^3, x^2, x, 1], [y^3, y^2, y, 1]), [3*z^2, 2*z, 1, 0])';

% av is the vector arrangement of A which is 4x4x4
% this ensures that the ordering of the parameters of av matches xv
av = flip(reshape(permute(A, [3,2,1]), [], 1));

grad = [
    sum(av .* xv);
    sum(av .* yv);
    sum(av .* zv);
    ];
end

