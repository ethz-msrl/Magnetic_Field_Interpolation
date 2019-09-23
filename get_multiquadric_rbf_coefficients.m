function [ C ] = get_multiquadric_rbf_coefficients( positions, values, eps )
%GET_RBF_COEFFICIENTS Calculates the RBF coefficients for a 3D vector field

    Nd = length(positions);
    
    xd = positions(:,1);
    yd = positions(:,2);
    zd = positions(:,3);

    [jj, ii] = meshgrid(1:Nd, 1:Nd);

    K = sqrt(1 + eps * ((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2 + ...
        (zd(ii) - zd(jj)).^2));
    
    C = K \ values;
    
end
