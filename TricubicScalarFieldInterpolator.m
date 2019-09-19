classdef TricubicScalarFieldInterpolator < FieldInterpolator
    %SIMPLETRICUBICINTERPOLATOR Performs tricubic interpolation separately
    % on the three field dimensions
    
    properties (SetAccess = private)
        M
        VG
        Nx
        Ny
        Nz
        dvg_dx
        dvg_dy
        dvg_dz
        dvg_dxy
        dvg_dxz
        dvg_dyz
        dvg_dxyz
        Coefs
        Steps
    end
    
    methods
        function obj = TricubicScalarFieldInterpolator(nodes, values, M)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            obj.M = M;
            obj.Nx = size(obj.NodePositions, 1);
            obj.Ny = size(obj.NodePositions, 2);
            obj.Nz = size(obj.NodePositions, 3);
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            obj.VG = permute(vg, [3,2,1,4]);
            
            [obj.dvg_dy, ~, obj.dvg_dz] = gradient(obj.VG);
            [~, ~, obj.dvg_dyz] = gradient(obj.dvg_dy);
            
            getAllCoefficients(obj);
            
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            obj.Steps = (pmax - pmin) ./ (size(obj.NodePositions(:,:,:,1)) - 1);

 
        end
        
       function [ix, iy, iz, x, y, z] = getIndices(obj, position)
            pn = getNormalizedPositions(obj, position);
            
            ix = min(floor(pn(1)) + 1, size(obj.NodePositions, 1) - 1);
            iy = min(floor(pn(2)) + 1, size(obj.NodePositions, 2) - 1);
            iz = min(floor(pn(3)) + 1, size(obj.NodePositions, 3) - 1);
            
            x = pn(1) + 1 - ix;
            y = pn(2) + 1 - iy;
            z = pn(3) + 1 - iz;
        end
        
        function getAllCoefficients(obj)
            obj.Coefs = zeros(obj.Nx-1, obj.Ny-1, obj.Nz-1, 64);
            for ix = 1:(obj.Nx-1)
                for iy = 1:(obj.Ny-1)
                    for iz = 1:(obj.Nz-1)
                        D = zeros(64,1);
                        D(9:end) = squeeze([
                        % f
                            reshape([
                            obj.VG(ix, iy, iz, :);
                            obj.VG(ix, iy, iz+1, :);
                            obj.VG(ix, iy+1, iz, :);
                            obj.VG(ix, iy+1, iz+1, :);
                            obj.VG(ix+1, iy, iz, :);
                            obj.VG(ix+1, iy, iz+1, :);
                            obj.VG(ix+1, iy+1, iz, :);
                            obj.VG(ix+1, iy+1, iz+1, :)], 24,1) ;
                       % dfx/dy
                            obj.dvg_dy(ix, iy, iz, 1);
                            obj.dvg_dy(ix, iy, iz+1, 1);
                            obj.dvg_dy(ix, iy+1, iz, 1);
                            obj.dvg_dy(ix, iy+1, iz+1, 1);
                            obj.dvg_dy(ix+1, iy, iz, 1);
                            obj.dvg_dy(ix+1, iy, iz+1, 1);
                            obj.dvg_dy(ix+1, iy+1, iz, 1);
                            obj.dvg_dy(ix+1, iy+1, iz+1, 1);
                        % dfx/dz
                            obj.dvg_dz(ix, iy, iz, 1);
                            obj.dvg_dz(ix, iy, iz+1, 1);
                            obj.dvg_dz(ix, iy+1, iz, 1);
                            obj.dvg_dz(ix, iy+1, iz+1, 1);
                            obj.dvg_dz(ix+1, iy, iz, 1);
                            obj.dvg_dz(ix+1, iy, iz+1, 1);
                            obj.dvg_dz(ix+1, iy+1, iz, 1);
                            obj.dvg_dz(ix+1, iy+1, iz+1, 1);
                        % dfy/dz
                            obj.dvg_dz(ix, iy, iz, 2);
                            obj.dvg_dz(ix, iy, iz+1, 2);
                            obj.dvg_dz(ix, iy+1, iz, 2);
                            obj.dvg_dz(ix, iy+1, iz+1, 2);
                            obj.dvg_dz(ix+1, iy, iz, 2);
                            obj.dvg_dz(ix+1, iy, iz+1, 2);
                            obj.dvg_dz(ix+1, iy+1, iz, 2);
                            obj.dvg_dz(ix+1, iy+1, iz+1, 2);
                        % d2fx/dydz
                            obj.dvg_dyz(ix, iy, iz, 1);
                            obj.dvg_dyz(ix, iy, iz+1, 1);
                            obj.dvg_dyz(ix, iy+1, iz, 1);
                            obj.dvg_dyz(ix, iy+1, iz+1, 1);
                            obj.dvg_dyz(ix+1, iy, iz, 1);
                            obj.dvg_dyz(ix+1, iy, iz+1, 1);
                            obj.dvg_dyz(ix+1, iy+1, iz, 1);
                            obj.dvg_dyz(ix+1, iy+1, iz+1, 1);
                            ]);

                        %a_sol = obj.M(9:32,:) \ D(9:32,:);
                        a_sol = obj.M \ D;
                        obj.Coefs(ix, iy, iz, :, :) = [0;a_sol];
                    end
                end
            end
        end
        
        function field = getFieldAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            field = -tricubic_grad(A_sol, xe, ye, ze);
        end
        
        function gradient = getGradientAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            H = -tricubic_hess(A_sol, xe, ye, ze);
%             gradient = H;
            gradient = [H(1,:)/obj.Steps(1); ...
                H(2,:)/obj.Steps(2); ...
                H(3,:)/obj.Steps(3)];
        end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            normalized = reshape((pm - pmin) ./ obj.Steps, size(positions));
        end
    end
    
end

