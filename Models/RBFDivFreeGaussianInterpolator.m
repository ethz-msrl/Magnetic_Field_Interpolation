classdef RBFDivFreeGaussianInterpolator < FieldInterpolator
    %RBFDIVFREEGAUSSIANINTERPOLATOR Interpolates a 3D vector field by using
    % a Gaussian divergence-free matrix-valued kernel
    %   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron
    
    properties (SetAccess = private)
        Eps
        Coefs
        CondNumber
        Maxp
        Minp
    end
    
    methods
        function obj = RBFDivFreeGaussianInterpolator(nodes, values, eps)
            %   Args:
            %           nodes (4D array): the node positions.
            %           the dimensions are Nx,Ny,Nz,3
            %           values (4D array): the field values at the node
            %           positions. The dimensions are Nx,Ny,Nz,3
            %           eps (float): the weighting term of the Gaussian RBF
            obj.NodePositions = reshape(nodes, [], 3);
            % normalize node positions
            [obj.NodePositions(:,1), obj.NodePositions(:,2), obj.NodePositions(:,3), ...
                obj.Maxp, obj.Minp] = normalize_positions_minmax(obj.NodePositions(:,1),...
                obj.NodePositions(:,2), obj.NodePositions(:,3));
            obj.NodeValues = reshape(values, [], 3);
            obj.Eps = eps;
            [obj.Coefs, obj.CondNumber] = get_divfree_gaussian_rbf_coefficients(obj.NodePositions, ...
                obj.NodeValues, eps);
        end
        
        function field = getFieldAtPosition(obj, position)
            [position(1), position(2), position(3)] = ...
                normalize_positions_minmax(position(1), position(2), position(3), ...
                obj.Maxp, obj. Minp);
            field = evaluate_divfree_gaussian_rbf(position, obj.NodePositions, ...
                obj.Eps, obj.Coefs);
        end
        
        function gradient = getGradientAtPosition(obj, position)
            [position(1), position(2), position(3)] = ...
                normalize_positions_minmax(position(1), position(2), position(3), ...
                obj.Maxp, obj. Minp);
            gradient = evaluate_divfree_gaussian_rbf_gradient(position, ...
                obj.NodePositions, obj.Eps, obj.Coefs);
            % since we scaled positions by max - min, we need to also scale
            % the gradient
            gradient = gradient ./ repmat(obj.Maxp - obj.Minp, 3, 1);
        end
        
    end
    
end

