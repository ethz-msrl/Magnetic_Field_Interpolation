function [M, B_fun, G_fun] = get_tricubic_3d_matrix()
    % GETTRICUBIC3DMATRIX
    % Here we setup the linear system of constraints for doing tricubic
    % interpolation of a scalar field f(x,y,z) using field and derivative constraints 

    % We use the remaining constraints on derivatives of the fields that are
    % calculated using finite differences

    % the following derivatives are used
    % v = [f, df_dx, df_dy, df_dz, d2f_dxy, d2f_dxz, d2f_dyz, d3f_dxdydz]
    %
    % returns the matrix representing the system 
    % M * av = d
    % where d are stacked versions of v for all eight corners
    % d = [f(1), ..., f(8), df_dx(1) ... df_dx(8) ... ]
    % and av are the coefficients of the tricubic interpolation
    % sum_i=1^4 sum_j=1^4 sum_k=1^4 a_ijk x^i y^j z^k 
    %   M: 64x64 matrix to solve for coefficients
    %   B_fun: handle to matlab function for performing 1D interpolation
    %   with a vector of coefficients
    %   G_fun: handle to matlab function for performing gradient of 1D
    %   interpolation with a vector of coefficients
    %
    %   Copyright 2020, Samuel Charreyron

    syms x y z
    A = sym('a%d%d%d', [4, 4, 4], 'real');
    av = reshape(A, 1, []);

    tri_fun = symfun(tricubic(A, x, y, z), [x y z]);
    order_0_eqn = [
        tri_fun(0,0,0)==0,
        tri_fun(0,0,1)==0,
        tri_fun(0,1,0)==0,
        tri_fun(0,1,1)==0,
        tri_fun(1,0,0)==0,
        tri_fun(1,0,1)==0,
        tri_fun(1,1,0)==0,
        tri_fun(1,1,1)==0
        ];

    tri_fun_dx = symfun(diff(tricubic(A,x,y,z),x), [x y z]);
    tri_fun_dy = symfun(diff(tricubic(A,x,y,z),y), [x y z]);
    tri_fun_dz = symfun(diff(tricubic(A,x,y,z),z), [x y z]);

    order_1_eqn = [
        tri_fun_dx(0,0,0)==0,
        tri_fun_dx(0,0,1)==0,
        tri_fun_dx(0,1,0)==0,
        tri_fun_dx(0,1,1)==0,
        tri_fun_dx(1,0,0)==0,
        tri_fun_dx(1,0,1)==0,
        tri_fun_dx(1,1,0)==0,
        tri_fun_dx(1,1,1)==0,
        tri_fun_dy(0,0,0)==0,
        tri_fun_dy(0,0,1)==0,
        tri_fun_dy(0,1,0)==0,
        tri_fun_dy(0,1,1)==0,
        tri_fun_dy(1,0,0)==0,
        tri_fun_dy(1,0,1)==0,
        tri_fun_dy(1,1,0)==0,
        tri_fun_dy(1,1,1)==0,
        tri_fun_dz(0,0,0)==0,
        tri_fun_dz(0,0,1)==0,
        tri_fun_dz(0,1,0)==0,
        tri_fun_dz(0,1,1)==0,
        tri_fun_dz(1,0,0)==0,
        tri_fun_dz(1,0,1)==0,
        tri_fun_dz(1,1,0)==0,
        tri_fun_dz(1,1,1)==0
        ];

    trif_fun_dxy = symfun(diff(tri_fun_dx(x,y,z),y), [x y z]);
    trif_fun_dxz = symfun(diff(tri_fun_dx(x,y,z),z), [x y z]);
    trif_fun_dyz = symfun(diff(tri_fun_dy(x,y,z),z), [x y z]);

    order_2_eqn = [
        trif_fun_dxy(0,0,0)==0,
        trif_fun_dxy(0,0,1)==0,
        trif_fun_dxy(0,1,0)==0,
        trif_fun_dxy(0,1,1)==0,
        trif_fun_dxy(1,0,0)==0,
        trif_fun_dxy(1,0,1)==0,
        trif_fun_dxy(1,1,0)==0,
        trif_fun_dxy(1,1,1)==0,
        trif_fun_dxz(0,0,0)==0,
        trif_fun_dxz(0,0,1)==0,
        trif_fun_dxz(0,1,0)==0,
        trif_fun_dxz(0,1,1)==0,
        trif_fun_dxz(1,0,0)==0,
        trif_fun_dxz(1,0,1)==0,
        trif_fun_dxz(1,1,0)==0,
        trif_fun_dxz(1,1,1)==0,
        trif_fun_dyz(0,0,0)==0,
        trif_fun_dyz(0,0,1)==0,
        trif_fun_dyz(0,1,0)==0,
        trif_fun_dyz(0,1,1)==0,
        trif_fun_dyz(1,0,0)==0,
        trif_fun_dyz(1,0,1)==0,
        trif_fun_dyz(1,1,0)==0,
        trif_fun_dyz(1,1,1)==0
        ];

    tri_fun_dxyz = symfun(diff(trif_fun_dxy(x,y,z),z), [x y z]);

    order_3_eqn = [
        tri_fun_dxyz(0,0,0)==0,
        tri_fun_dxyz(0,0,1)==0,
        tri_fun_dxyz(0,1,0)==0,
        tri_fun_dxyz(0,1,1)==0,
        tri_fun_dxyz(1,0,0)==0,
        tri_fun_dxyz(1,0,1)==0,
        tri_fun_dxyz(1,1,0)==0,
        tri_fun_dxyz(1,1,1)==0,
        ];

    M = double(equationsToMatrix(vertcat(order_0_eqn, order_1_eqn, order_2_eqn, order_3_eqn), av));
    B_fun = matlabFunction(tri_fun, 'Vars', {x, y, z, av});
    
    G_fun = matlabFunction(gradient(tri_fun, [x,y,z]), 'Vars', {x, y, z, av});
end