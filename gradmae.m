function y = gradmae(observed, predicted )
% MEANABSERROR: computes the mean absolute deviation between observed and 
% predicted gradient matrices
%  
%   observed and prediction are assumed to have shape [Nx, Ny, Nz, 3, 3]
%   the returned MAE has shape [3,3]

    y = mean(reshape(abs(observed - predicted),[],9), 1);
    y = reshape(y, 3, 3);

end