function [C, I] = get_monomial_coeffs(BFun)
%GET_MONOMIAL_COEFFS gets the scalar coefficients

%% Will be deleted since it doensnt really work
% given a scalar polynomial function of the type
% f(av, x,y,z) = sum (c_ijk * a_ijk * x^l * y^m * z^n)
% where av = [a_ijk]i=1..p,j=1..r,k=1..s are symbolic quantities and c_ijk are scalars
% we factor out as 
% f = dv .* xv
% where dv is (cv .* av) permuted to match the indices of xv
%  
%   Args: BFun a symbolic function with av and x, y, z as variables
%   Returns: C: the vector of scalar coefficients
%            I: a vector of indices to permuate av before multiplying by cv
%            and xv
%
% 
    syms x y z real;
    A = sym('a%d%d%d', [4, 4, 4], 'real');
    %av = reshape(A, [], 1);
    av = flip(reshape(permute(A, [3,2,1]), [], 1));
    xv = kron(kron([x^3, x^2, x, 1], [y^3, y^2, y, 1]), [z^3, z^2, z, 1])';
    
    [C_, T_] = coeffs(BFun, [x y z]);

    [~, ia] = intersect(T_, xv);
    [~, ina] = setdiff(xv, T_);

    Mx = sym('m', [64, 1]);
    Mx(ia) = C_;
    Mx(ina) = 0;
    T = sym('t', [64, 1]);
    C = zeros(64,1);
    tz = [];
    for i =1:64
        [c,t] = coeffs(Mx(i), av);
        if ~ isempty(t)
            T(i) = t;
            C(i) = c;
        else
            T(i) = 0;
            tz = [tz, i];
        end
    end
    T(tz) = setdiff(av, T);
    [~,I] = sort(T);
end
