classdef BSpline3DInterpolator < FieldInterpolator
    %BSPLINE3DINTERPOLATOR Evaluates a 3D B-Spline of degree D with no constraints 
    
    properties 
        D
        C
        k_n
        k_m
        k_p
    end
    
    methods
        function obj = BSpline3DInterpolator(nodes, values, degree)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            obj.D = degree;
            
            obj.getCoefficients();
        end
        
        function getCoefficients(obj)
            
            % gets back the grid vectors from gridded node data
            % assumes obj.NodePositions is array of size Nx, Ny, Nz, 3
            xv = unique(obj.NodePositions(:,:,:,1));
            yv = unique(obj.NodePositions(:,:,:,2));
            zv = unique(obj.NodePositions(:,:,:,3));

            % this gets the appropriate knot sequence that satisfies the
            % Schoenberg-Whitney condition
            obj.k_n = aptknt(xv, obj.D);
            obj.k_m = aptknt(yv, obj.D);
            obj.k_p = aptknt(zv, obj.D);

            N = spcol(obj.k_n, obj.D, xv);
            M = spcol(obj.k_m, obj.D, yv);
            P = spcol(obj.k_p, obj.D, zv);

            Z = kron(kron(N,M), P);
            obj.C = Z \ reshape(obj.NodeValues, [], 3);
        end
        
        function field = getFieldAtPosition(obj, position)
            N = spcol(obj.k_n, obj.D, position(1));
            M = spcol(obj.k_m, obj.D, position(2));
            P = spcol(obj.k_p, obj.D, position(3));
            Z = kron(kron(N, M), P);
            field = (Z * obj.C)';
        end
        
        function fields = getFieldsAtPositions(obj, positions)
            % get fields at grid positions
            % positions is assumed to be a 4D array of size Nx, Ny, Nz, 3
            % returns a 4D array of field values of size Nx, Ny, Nz, 3
            
            % gets back the grid vectors from gridded node data
            % assumes obj.NodePositions is array of size Nx, Ny, Nz, 3
            xve = unique(positions(:,:,:,1));
            yve = unique(positions(:,:,:,2));
            zve = unique(positions(:,:,:,3));
            
            N = spcol(obj.k_n, obj.D, xve);
            M = spcol(obj.k_m, obj.D, yve);
            P = spcol(obj.k_p, obj.D, zve);

            Z = kron(kron(N,M), P);
            fields = reshape(Z * obj.C, size(positions, 1), size(positions,2), ...
                size(positions,3), 3);
        end
    end
    
end

