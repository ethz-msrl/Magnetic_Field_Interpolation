function [M] = get_tricubic_scalar_field_M()

    % Here we setup the linear system of constraints for doing tricubic
    % interpolation of a scalar potential that is defined using magnetic 
    % field values. We enforce the following constaint at the cube corners

    % -grad(phi) = B
    % lapl(phi) = 0 (from 0 divergence property)

    % We use the remaining constraints on derivatives of the fields that are
    % calculated using finite differences

    % the following derivatives are used
    % dBx/dy, dBx/dz, dBy/dz d2Bx/dydz

    % returns the matrix representing the system 
    % M * av = d
    % where d are stacked versions of v for all eight corners
    % d = [lapl(1), ..., lapl(8), -Bx(1) -By(1) - Bz(1) ... 
    % dBx_dy(1) ... dBx_dy(8) ... d2Bx_dydz(1) ... d2Bx_dydz(8) ]
    % and av are the coefficients of the tricubic interpolation
    % sum_i=1^4 sum_j=1^4 sum_k=1^4 a_ijk x^i y^j z^k 

    syms x y z
    A = sym('a%d%d%d', [4, 4, 4]);
    av = reshape(A, 1, []);

    % we remove the constant element of the scalar potential
    av = av(2:end);

    % get equations for laplace equation
    lapl = symfun(diff(tricubic(A,x,y,z),x,2) + diff(tricubic(A,x,y,z),y,2) + diff(tricubic(A,x,y,z),z,2), [x y z]);

    lapl_equations = [
        lapl(0,0,0)==0, 
        lapl(0,0,1)==0, 
        lapl(0,1,0)==0, 
        lapl(0,1,1)==0,
        lapl(1,0,0)==0,
        lapl(1,0,1)==0,
        lapl(1,1,0)==0,
        lapl(1,1,1)==0];

    %grad = symfun(gradient(tricubic(A,x,y,z), [x y z]), [x y z]);
    
    gradx = symfun(diff(tricubic(A,x,y,z),x), [x y z]);
    grady = symfun(diff(tricubic(A,x,y,z),y), [x y z]);
    gradz = symfun(diff(tricubic(A,x,y,z),z), [x y z]);
    
    grad_equations = [
        -gradx(0,0,0) == 0,
        -gradx(0,0,1) == 0,
        -gradx(0,1,0) == 0,
        -gradx(0,1,1) == 0,
        -gradx(1,0,0) == 0,
        -gradx(1,0,1) == 0,
        -gradx(1,1,0) == 0,
        -gradx(1,1,1) == 0,
        -grady(0,0,0) == 0,
        -grady(0,0,1) == 0,
        -grady(0,1,0) == 0,
        -grady(0,1,1) == 0,
        -grady(1,0,0) == 0,
        -grady(1,0,1) == 0,
        -grady(1,1,0) == 0,
        -grady(1,1,1) == 0,
        -gradz(0,0,0) == 0,
        -gradz(0,0,1) == 0,
        -gradz(0,1,0) == 0,
        -gradz(0,1,1) == 0,
        -gradz(1,0,0) == 0,
        -gradz(1,0,1) == 0,
        -gradz(1,1,0) == 0,
        -gradz(1,1,1) == 0,
        ];

    Bx_y = symfun(diff(-gradx, y), [x y z]);
    Bx_z = symfun(diff(-gradx, z), [x y z]);
    By_z = symfun(diff(-grady, z), [x y z]);

    Gxyz = symfun(diff(Bx_y, z), [x y z]);

    hess_equations = [
        Bx_y(0,0,0) == 0,
        Bx_y(0,0,1) == 0,
        Bx_y(0,1,0) == 0,
        Bx_y(0,1,1) == 0,
        Bx_y(1,0,0) == 0,
        Bx_y(1,0,1) == 0,
        Bx_y(1,1,0) == 0,
        Bx_y(1,1,1) == 0,
        Bx_z(0,0,0) == 0,
        Bx_z(0,0,1) == 0,
        Bx_z(0,1,0) == 0,
        Bx_z(0,1,1) == 0,
        Bx_z(1,0,0) == 0,
        Bx_z(1,0,1) == 0,
        Bx_z(1,1,0) == 0,
        Bx_z(1,1,1) == 0,
        By_z(0,0,0) == 0,
        By_z(0,0,1) == 0,
        By_z(0,1,0) == 0,
        By_z(0,1,1) == 0,
        By_z(1,0,0) == 0,
        By_z(1,0,1) == 0,
        By_z(1,1,0) == 0,
        By_z(1,1,1) == 0,
        ];

    hess_der_equations = [
        Gxyz(0,0,0) == 0,
        Gxyz(0,0,1) == 0,
        Gxyz(0,1,0) == 0,
        Gxyz(0,1,1) == 0,
        Gxyz(1,0,0) == 0,
        Gxyz(1,0,1) == 0,
        Gxyz(1,1,0) == 0,
        Gxyz(1,1,1) == 0
        ];

    % There are 64 equations but we get a matrix of rank 63 since there
    % is no constant term to the scalar potential
    M = double(equationsToMatrix(vertcat(lapl_equations, grad_equations, ...
    hess_equations, hess_der_equations), av));

end
