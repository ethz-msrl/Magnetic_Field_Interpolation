function [ C ] = get_divfree_multiquadric_rbf_coefficients( positions, values, eps )
%GET_RBF_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    Nd = length(positions);
    
    xd = positions(:,1);
    yd = positions(:,2);
    zd = positions(:,3);

    [jj, ii] = meshgrid(1:Nd, 1:Nd);

    K = sqrt(1 + eps * ((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2 + ...
        (zd(ii) - zd(jj)).^2));
    
    K3 = K.^3;

    psi_11 = 2*eps ./ K - eps^2 * ((yd(ii) - yd(jj)).^2 + (zd(ii) - zd(jj)).^2) ./ K3 ;
    psi_22 = 2*eps ./ K - eps^2 * ((xd(ii) - xd(jj)).^2 + (zd(ii) - zd(jj)).^2) ./ K3 ;
    psi_33 = 2*eps ./ K - eps^2 * ((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2) ./ K3 ;

    psi_12 = eps^2 * (xd(ii) - xd(jj)) .* (yd(ii) - yd(jj)) ./ K3;
    psi_13 = eps^2 * (xd(ii) - xd(jj)) .* (zd(ii) - zd(jj)) ./ K3;
    psi_23 = eps^2 * (yd(ii) - yd(jj)) .* (zd(ii) - zd(jj)) ./ K3;

    A = [
        psi_11, psi_12, psi_13;
        psi_12, psi_22, psi_23;
        psi_13, psi_23, psi_33;
        ];
    
    D = reshape(values, [], 1);
    
    c = A \ D;
    C = reshape(c, Nd, 3);
end

