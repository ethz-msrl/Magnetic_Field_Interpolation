classdef SimpleRBFInterpolator < FieldInterpolator
    %SIMPLERBFINTERPOLATOR Interpolates a 3D vector field by creating a Gaussian RBF
    %   for each dimension
    %   B_i = sum_i^N c_i * exp(-eps * ||pos - node_i||^2)
    %   for N nodes
    
    properties (SetAccess = private)
        Eps
        Coefs
    end
    
    methods
        function obj = SimpleRBFInterpolator(nodes, values, eps)
            % Constructor
            %   Args:
            %           nodes (4D array): the node positions.
            %           the dimensions are Nx,Ny,Nz,3
            %           values (4D array): the field values at the node
            %           positions. The dimensions are Nx,Ny,Nz,3
            %           eps (float): the weighting term of the Gaussian RBF
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            obj.Eps = eps;
            obj.Coefs = get_rbf_coefficients(nodes, values, eps);
        end
        
        function field = getFieldAtPosition(obj, position)
            field = evaluate_rbf(position, obj.NodePositions, obj.Eps, obj.Coefs);
        end
    end
    
end

