classdef FieldInterpolator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        NodePositions
        NodeValues
    end
    
    methods (Abstract)
        field = getFieldAtPosition(position)
    end
    
    methods
        % non vectorized function to get fields at several positions
        % override with better version if you can vectorize
        function fields = getFieldsAtPositions(obj, positions)
            sz = size(positions);
            if sz(end) ~= 3
                error('last dimension of positions must be 3')
            end

            if ~ismatrix(positions)
                positions = reshape(positions, [], 3);
            end
            
            N = size(positions, 1);
            fields = zeros(size(positions));
            for i =1:N
                fields(i,:) = getFieldAtPosition(obj, positions(i,:));
            end
            
            if numel(sz) ~= numel(size(positions))
                fields = reshape(fields, sz);
            end
        end
    end
    
end

