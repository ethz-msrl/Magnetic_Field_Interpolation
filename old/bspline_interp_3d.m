function [ output_args ] = bspline_interp_3d( knots, v, eval, m )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

xv = knots(:,1);
yv = knots(:,2);
zv = knots(:,3);
Ni = spcol(padarray(xv, [2,0], 'replicate'), m, eval(:,1));
Mj = spcol(padarray(yv, [2,0], 'replicate'), m, eval(:,2));
Pk = spcol(padarray(zv, [2,0], 'replicate'), m, eval(:,2));


end

