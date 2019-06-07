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
%         function obj = FieldInterpolator(nodes, values)
%             obj.NodePositions = nodes;
%             obj.NodeValues = values;
%         end
        
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
        
        function gradient = getGradientAtPosition(obj, position)
            error('Not yet implemented');
        end
        
        function gradients = getGradientsAtPositions(obj, positions)
            sz = size(positions);
            if sz(end) ~= 3
                error('last dimension of positions must be 3')
            end

            if ~ismatrix(positions)
                positions = reshape(positions, [], 3);
            end
            
            N = size(positions, 1);
            gradients = zeros(N, 3, 3);
            for i =1:N
                gradients(i,:,:) = getGradientAtPosition(obj, positions(i,:));
            end
            
            % making sure that the shape of the gradients matches the
            % positions shape
            if numel(sz) ~= numel(size(positions))
                gradients = reshape(gradients, [sz, 3]);
            end
        end
    end
    
end

