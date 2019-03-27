function [ C ] = get_rbf_coefficients( positions, values, eps )
%GET_RBF_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    Nd = length(positions);
    
    xd = positions(:,1);
    yd = positions(:,2);
    zd = positions(:,3);

    [jj, ii] = meshgrid(1:Nd, 1:Nd);

    K = exp(-eps*((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2 + ...
        (zd(ii) - zd(jj)).^2));
    
    C = K \ values;
    
end
