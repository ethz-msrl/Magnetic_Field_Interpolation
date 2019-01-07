function [ C ] = get_divfree_rbf_coefficients( positions, values, eps )
%GET_RBF_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    Nd = length(positions);
    
    xd = positions(:,1);
    yd = positions(:,2);
    zd = positions(:,3);

    [jj, ii] = meshgrid(1:Nd, 1:Nd);

    K = exp(-eps*((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2 + ...
        (zd(ii) - zd(jj)).^2));

    psi_11 = (4*eps - 4*eps^2 * ((yd(ii) - yd(jj)).^2 + (zd(ii) - zd(jj)).^2)) .* K;
    psi_22 = (4*eps - 4*eps^2 * ((xd(ii) - xd(jj)).^2 + (zd(ii) - zd(jj)).^2)) .* K;
    psi_33 = (4*eps - 4*eps^2 * ((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2)) .* K;

    psi_12 = 4 * eps^2 * (xd(ii) - xd(jj)) .* (yd(ii) - yd(jj)) .* K;
    psi_13 = 4 * eps^2 * (xd(ii) - xd(jj)) .* (zd(ii) - zd(jj)) .* K;
    psi_23 = 4 * eps^2 * (yd(ii) - yd(jj)) .* (zd(ii) - zd(jj)) .* K;

    A = [
        psi_11, psi_12, psi_13;
        psi_12, psi_22, psi_23;
        psi_13, psi_23, psi_33;
        ];
    
    D = reshape(values, [], 1);
    
    c = A \ D;
    C = reshape(c, Nd, 3);
end

