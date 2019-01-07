function [M, B_fun, G_fun] = get_tricubic_div_free_matrix()
%GET_TRICUBIC_DIV_FREE_MATRIX Derives the divergence free tricubic
%   interpolation matrix and function handle
%   The matrix allows you to solve the linear system 
%   M * a = d
%   where d is the constraint vector where the first 24 rows are the
%   magnetic fields at the 8 corners and a is the solution vector of the
%   tribubic coefficients
%
%   Divergence free is obtained by interpolating a 3D vector field A and
%   taking its curl. The vector field is constrained such that its curl
%   equals the field values at the 8 corners, and that the curl of the curl
%   is also zero at the eight corners (Ampère's law with no  enclosed
%   currents)
%   Returns:
%       - the interpolation matrix M
%       - the interpolation function 
%       - the gradient function
    
    syms x y z
    Ax = sym('ax%d%d%d', [4, 4, 4]);
    axv = reshape(Ax, [], 1);

    Ay = sym('ay%d%d%d', [4, 4, 4]);
    ayv = reshape(Ay, [], 1);

    Az = sym('az%d%d%d', [4, 4, 4]);
    azv = reshape(Az, [], 1);

    Af = [tricubic(Ax, x, y, z);
    tricubic(Ay, x, y, z);
    tricubic(Az, x, y, z)];

    curlAf = symfun(curl(Af, [x,y,z]), [x,y,z]);
    B_fun = matlabFunction(curlAf, 'Vars', {x, y, z, vertcat(axv, ayv, azv)});
    G_x = diff(curlAf, x)';
    G_y = diff(curlAf, y)';
    G_z = diff(curlAf, z)';
    G_fun = matlabFunction([G_x; G_y; G_z], 'Vars', {x, y, z, vertcat(axv, ayv, azv)});

    % enforcing field at control points
    curl_equations = [
        curlAf(0,0,0) == 0,
        curlAf(0,0,1) == 0,
        curlAf(0,1,0) == 0,
        curlAf(0,1,1) == 0,
        curlAf(1,0,0) == 0,
        curlAf(1,0,1) == 0,
        curlAf(1,1,0) == 0,
        curlAf(1,1,1) == 0
        ];

    dB_dx = diff(curlAf(x,y,z),x);
    dB_dy = diff(curlAf(x,y,z),y);
    dB_dz = diff(curlAf(x,y,z),z);

    dBx_dx(x,y,z) = dB_dx(1);
    dBy_dy(x,y,z) = dB_dy(2);
    dBz_dz(x,y,z) = dB_dz(3);
    dBy_dx(x,y,z) = dB_dx(3);
    dBx_dy(x,y,z) = dB_dy(1);
    dBx_dz(x,y,z) = dB_dz(1);
    dBz_dx(x,y,z) = dB_dx(3);
    dBy_dz(x,y,z) = dB_dz(2);
    dBz_dy(x,y,z) = dB_dy(3);

    grad_equations = [];
    for i = 1:2
        for j = 1:2
            for k = 1:2
                grad_equations = [
                    grad_equations;
                    % enforcing gradient measurements
                    dBx_dx(i-1, j-1, k-1) == 0;
                    dBx_dy(i-1, j-1, k-1) == 0;
                    dBx_dz(i-1, j-1, k-1) == 0;
                    dBy_dx(i-1, j-1, k-1) == 0;
                    dBy_dy(i-1, j-1, k-1) == 0;
                    dBy_dz(i-1, j-1, k-1) == 0;
                    dBz_dx(i-1, j-1, k-1) == 0;
                    dBz_dy(i-1, j-1, k-1) == 0;
                    dBz_dz(i-1, j-1, k-1) == 0;
                    % enforcing zero curl at the control points
%                     dBx_dy(i-1, j-1, k-1) - dBy_dx(i-1, j-1, k-1) == 0;
%                     dBx_dz(i-1, j-1, k-1) - dBz_dx(i-1, j-1, k-1) == 0;
%                     dBy_dz(i-1, j-1, k-1) - dBz_dy(i-1, j-1, k-1) == 0;
                    ];
            end
        end
    end

    M = double(equationsToMatrix(vertcat(curl_equations, grad_equations), ...
        vertcat(axv, ayv, azv)));

end

