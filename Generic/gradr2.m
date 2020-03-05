function r2 = gradr2(observed, predicted)
% GRADR2 computes the R2 coefficient of determination between observed and
% predicted gradient matrices
%
%   observed and prediction are assumed to have shape [Nx, Ny, Nz, 3, 3]
%   the returned MAE has shape [3,3]
    ss_res = sum(reshape(observed - predicted,[],9).^2, 1);
    ss_tot = sum(reshape(observed, [], 9).^2, 1);
    r2 = reshape(1 - ss_res ./ ss_tot, 3, 3);
end

