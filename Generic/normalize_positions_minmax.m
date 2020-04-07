function [xgn, ygn, zgn, maxp, minp] = normalize_positions_minmax(xg, yg, zg, varargin)
% NORMALIZE_POSITIONS_MINMAX normalizes the positions to be between 0 and 1
%  xg, yg, zg are arrays to be normalized
%
%  normalize_positions_minmax(xg, yg, zg) will normalize based on the
%  maximums and miniums of xg, yg,zg
%
%  normalize_positions_minmax(xg, yg, zg, maxp, minp) will normalized based
%  on maxp and minp
%
%  if minp and maxp are optional 3D arrays with the min/max values used for
%  normalization. If they are supplied, they will be used. If not, the min
%  and max values will be determined from the xg, yg, zg arrays
%  the maxp and minp are returned for use in normalization later
%
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

    if nargin == 5
        maxp = varargin{1};
        minp = varargin{2};
    elseif nargin == 3
        maxp = [max(xg(:)), max(yg(:)), max(zg(:))];
        minp = [min(xg(:)), min(yg(:)), min(zg(:))];
    else
        error('wrong number of arguments');
    end

    xgn = (xg - minp(1)) / (maxp(1) - minp(1));
    ygn = (yg - minp(2)) / (maxp(2) - minp(2));
    zgn = (zg - minp(3)) / (maxp(3) - minp(3));
end
