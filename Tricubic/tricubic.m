function [ y ] = tricubic( A, x, y, z )
%TRICUBIC Calculates the tricubic function using numerical computations
%instead of symbolically in tricubic
%
%  This should allow you to compute things much faster than using symbolic
%  substitutions
%
%   Arguments:
%               A: as in tricubic, it should be 4x4x4 matrix
%               x: should be a scalar
%               y: should be a scalar
%               z: should be a scalar
%   Returns:
%               y: a scalar function corresponding to the tricubic value
%   Copyright 2020, Samuel Charreyron

xv = kron(kron([x^3, x^2, x, 1], [y^3, y^2, y, 1]), [z^3, z^2, z, 1])';
% av is the vector arrangement of A which is 4x4x4
% this ensures that the ordering of the parameters of av matches xv
av = flip(reshape(permute(A, [3,2,1]), [], 1));

y = sum(av .* xv);

end

