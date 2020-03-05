classdef BSplineLaplacianInterpolator < FieldInterpolator
    %BSplineLaplacianInterpolator Interpolates using B-Splines with the
    %laplacian constrained to be zero on the grid of knot points
    %   Copyright 2020, Samuel Charreyron
    
    properties 
        D
        C
        k_n
        k_m
        k_p
    end
    
    methods
        function obj = BSplineLaplacianInterpolator(nodes, values, degree)
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
            
            %Zi = kron(P, kron(M, N));
            Zi = kron(N, kron(M, P));
            % we align the coefficients as [Cijkx Cijky; Cijkz]^T
            Z = blkdiag(Zi, Zi, Zi);
            
            % getting the first derivatives of the bases
            temp = spcol(obj.k_n, obj.D, brk2knt(xv, 2));
            Ndot = temp(2:2:end, :);

            temp = spcol(obj.k_m, obj.D, brk2knt(yv, 2));
            Mdot = temp(2:2:end, :);

            temp = spcol(obj.k_p, obj.D, brk2knt(zv, 2));
            Pdot = temp(2:2:end, :);
            
%             Zx = kron(P, kron(M, Ndot));
%             Zy = kron(P, kron(Mdot, N));
%             Zz = kron(Pdot, kron(M, N));
            Zx = kron(Ndot, kron(M, P));
            Zy = kron(N, kron(Mdot, P));
            Zz = kron(N, kron(M, Pdot));


            % This version minimizes the divergence and curl at the measurement
            % positions
            Q = [Zx, Zy, Zz; zeros(size(Zx)), -Zz, Zy; Zz, zeros(size(Zy)), -Zx; -Zy, Zx, zeros(size(Zz))];

            E = reshape(obj.NodeValues, [], 1);
            
            opt = optimoptions('quadprog', 'Algorithm', 'interior-point-convex',...
                'Display', 'off');
            
            coef = quadprog(Z'*Z, -E'*Z, [], [], Q, zeros(size(Q,1),1), ...
                [], [], [], opt);
            obj.C = reshape(coef, [], 3);

        end
        
        function field = getFieldAtPosition(obj, position)
            N = spcol(obj.k_n, obj.D, position(1));
            M = spcol(obj.k_m, obj.D, position(2));
            P = spcol(obj.k_p, obj.D, position(3));
            
            Z = kron(N, kron(M,P));

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
            
            Z = kron(N, kron(M,P));
            %Z = kron(P, kron(M, N));
            
%             temp = spcol(obj.k_n, obj.D, brk2knt(xve, 2));
%             Ndot = temp(2:2:end, :);
% 
%             temp = spcol(obj.k_m, obj.D, brk2knt(yve, 2));
%             Mdot = temp(2:2:end, :);
% 
%             temp = spcol(obj.k_p, obj.D, brk2knt(zve, 2));
%             Pdot = temp(2:2:end, :);
% 
%             Zx = kron(P, kron(M, Ndot));
%             Zy = kron(P, kron(Mdot, N));
%             Zz = kron(Pdot, kron(M, N));
% 
%             Q = [Zx, Zy, Zz];

            fields = reshape(Z * obj.C, size(positions, 1), size(positions,2), ...
                size(positions,3), 3);
        end
        
        function gradient = getGradientAtPosition(obj, position)
            gradient = zeros(3,3);
            N = spcol(obj.k_n, obj.D, brk2knt(position(1), 2));
            M = spcol(obj.k_m, obj.D, position(2));
            P = spcol(obj.k_p, obj.D, position(3));
            Z = kron(kron(N, M), P);
            temp = Z * obj.C;
            gradient(:,1) = temp(2:2:end,:);
            N = spcol(obj.k_n, obj.D, position(1));
            M = spcol(obj.k_m, obj.D, brk2knt(position(2), 2));
            P = spcol(obj.k_p, obj.D, position(3));
            Z = kron(kron(N, M), P);
            temp = Z * obj.C;
            gradient(:,2) = temp(2:2:end,:);
            N = spcol(obj.k_n, obj.D, position(1));
            M = spcol(obj.k_m, obj.D, position(2));
            P = spcol(obj.k_p, obj.D, brk2knt(position(3), 2));
            Z = kron(kron(N, M), P);
            temp = Z * obj.C;
            gradient(:,3) = temp(2:2:end,:);
        end
    end
    
end

