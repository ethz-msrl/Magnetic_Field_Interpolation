function y = gradrmse(observed, predicted )
% GRADRMSE: computes the root mean square error between observed and 
% predicted gradient matrices
%  
%   observed and prediction are assumed to have shape [Nx, Ny, Nz, 3, 3]
%   the returned MAE has shape [3,3]
% 
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

    y = sqrt(mean(reshape((observed - predicted).^2,[],9), 1));
    y = reshape(y, 3, 3);

end
