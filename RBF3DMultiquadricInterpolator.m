classdef RBF3DMultiquadricInterpolator < FieldInterpolator
    %SIMPLERBFINTERPOLATOR Interpolates a 3D vector field by creating a multiquadric RBF
    %   for each dimension
    %   B_i = sum_i^N c_i sqrt(1 + eps * ||pos - node_i||^2)
    %   for N nodes
    
    properties (SetAccess = private)
        Eps
        Coefs
        CondNumber
    end
    
    methods
        function obj = RBF3DMultiquadricInterpolator(nodes, values, eps)
            % Constructor
            %   Args:
            %           nodes (4D array): the node positions.
            %           the dimensions are Nx,Ny,Nz,3
            %           values (4D array): the field values at the node
            %           positions. The dimensions are Nx,Ny,Nz,3
            %           eps (float): the weighting term of the Gaussian RBF
            obj.NodePositions = reshape(nodes, [], 3);
            obj.NodeValues = reshape(values, [], 3);
            obj.Eps = eps;
            [obj.Coefs, obj.CondNumber] = get_multiquadric_rbf_coefficients(obj.NodePositions, obj.NodeValues, eps);
        end
        
        function field = getFieldAtPosition(obj, position)
            field = evaluate_multiquadric_rbf(position, obj.NodePositions, obj.Eps, obj.Coefs);
        end
        
        function gradient = getGradientAtPosition(obj, position)
            [~, gradient] = evaluate_multiquadric_rbf(position, obj.NodePositions, obj.Eps, obj.Coefs);
        end
    end
    
end

