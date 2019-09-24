classdef RBFDivFreeMultiquadricInterpolator < FieldInterpolator
    %DIVFREERBFINTERPOLATOR Interpolates a 3D vector field by using a 
    % 3D divergence-free kernel
    
    properties (SetAccess = private)
        Eps
        Coefs
        CondNumber
    end
    
    methods
        function obj = RBFDivFreeMultiquadricInterpolator(nodes, values, eps)
            %   Args:
            %           nodes (4D array): the node positions.
            %           the dimensions are Nx,Ny,Nz,3
            %           values (4D array): the field values at the node
            %           positions. The dimensions are Nx,Ny,Nz,3
            %           eps (float): the weighting term of the Gaussian RBF
            obj.NodePositions = reshape(nodes, [], 3);
            obj.NodeValues = reshape(values, [], 3);
            obj.Eps = eps;
            [obj.Coefs, obj.CondNumber] = get_divfree_multiquadric_rbf_coefficients(obj.NodePositions, ...
                obj.NodeValues, eps);
        end
        
        function field = getFieldAtPosition(obj, position)
            field = evaluate_divfree_multiquadric_rbf(position, obj.NodePositions, ...
                obj.Eps, obj.Coefs);
        end
        
        function gradient = getGradientAtPosition(obj, position)
            error('not done yet');
%             gradient = evaluate_divfree_rbf_gradient(position, ...
%             obj.NodePositions, obj.Eps, obj.Coefs);
        end
    end
    
end

