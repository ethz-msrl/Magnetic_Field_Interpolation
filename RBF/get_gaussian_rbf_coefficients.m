function [ C, cond_number ] = get_gaussian_rbf_coefficients( positions, values, eps )
%GET_GAUSSIAN_RBF_COEFFICIENTS Calculates the RBF coefficients for three
%   scalar RBF interpolants with a Gaussian kernel
%
%   Arguments:
%   position: 3D vector at which to evaluate the RBF
%   nodes: Nx3 matrix of node positions
%   eps: scalar value of the decay param of the RBF
%   
%   Returns:
%   C: Nx3 matrix with coefficients
%   cond_number: float with the condition number of the interpolation
%   matrix
%
%   Copyright 2020, Samuel Charreyron

    Nd = length(positions);
    
    xd = positions(:,1);
    yd = positions(:,2);
    zd = positions(:,3);

    [jj, ii] = meshgrid(1:Nd, 1:Nd);

    K = exp(-eps*((xd(ii) - xd(jj)).^2 + (yd(ii) - yd(jj)).^2 + ...
        (zd(ii) - zd(jj)).^2));
    
    cond_number = cond(K);
    
    C = K \ values;
    
end
