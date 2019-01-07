function [ B ] = evaluate_tricubic_div_free( fields, x, y, z, M, B_fun )
%EVALUATE_TRICUBIC_DIV_FREE Evalutes the divergence free tricubic
%interpolation
%   Takes as inputs:
%       - fields: a 8x3 matrix of field values at the 8 corners of the
%       interpolation cube. The order of the corners is where the order of
%       the columns is x, y, z
%           0,0,0
%           0,0,1
%           0,1,0
%           0,1,1
%           1,0,0
%           1,0,1
%           1,1,0
%           1,1,1
%       - x: the x position between 0 and 1
%       - y: the y position between 0 and 1
%       - z: the z position between 0 and 1
%       - the interpolation matrix (see get_tribubic_div_free_matrix)
%       - the magnetic field function handle (see get_tribuic_div_free_matrix)
%   Returns:
%       - the interpolated magnetic field at (x, y, z)
    
    D = zeros(48,1);
    D(1:24) = [
        reshape(fields', 8*3, 1);
        ];
    %     Bx_x';
    %     Bx_y';
    %     Bx_z';
    %     By_y';
    %     By_z'];


    a_sol = M \ D;
    B = B_fun(x, y, z, a_sol);
end

