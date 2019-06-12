classdef Tricubic3DInterpolator < FieldInterpolator
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
        function obj = Tricubic3DInterpolator(nodes, values, M)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            obj.Nx = size(obj.NodePositions, 1);
            obj.Ny = size(obj.NodePositions, 2);
            obj.Nz = size(obj.NodePositions, 3);
            obj.M = M;
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            obj.VG = permute(vg, [3,2,1,4]);
            [obj.dvg_dy, obj.dvg_dx, obj.dvg_dz] = gradient(obj.VG);
            [obj.dvg_dxy, ~, obj.dvg_dxz] = gradient(obj.dvg_dx);
            [~,~,obj.dvg_dyz] = gradient(obj.dvg_dy);
            [~,~,obj.dvg_dxyz] = gradient(obj.dvg_dxy);
            
            getAllCoefficients(obj);
            
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            obj.Steps = (pmax - pmin) ./ (size(obj.NodePositions(:,:,:,1)) - 1);
        end
        
        function getAllCoefficients(obj)
            obj.Coefs = zeros(obj.Nx-1, obj.Ny-1, obj.Nz-1, 64, 3);
            for ix = 1:(obj.Nx-1)
                for iy = 1:(obj.Ny-1)
                    for iz = 1:(obj.Nz-1)
                        D = squeeze([
                            % f
                            obj.VG(ix, iy, iz, :);
                            obj.VG(ix, iy, iz+1, :);
                            obj.VG(ix, iy+1, iz, :);
                            obj.VG(ix, iy+1, iz+1, :);
                            obj.VG(ix+1, iy, iz, :);
                            obj.VG(ix+1, iy, iz+1, :);
                            obj.VG(ix+1, iy+1, iz, :);
                            obj.VG(ix+1, iy+1, iz+1, :);
                            % df/dx
                            obj.dvg_dx(ix, iy, iz, :);
                            obj.dvg_dx(ix, iy, iz+1, :);
                            obj.dvg_dx(ix, iy+1, iz, :);
                            obj.dvg_dx(ix, iy+1, iz+1, :);
                            obj.dvg_dx(ix+1, iy, iz, :);
                            obj.dvg_dx(ix+1, iy, iz+1, :);
                            obj.dvg_dx(ix+1, iy+1, iz, :);
                            obj.dvg_dx(ix+1, iy+1, iz+1, :);
                            % df/dy
                            obj.dvg_dy(ix, iy, iz, :);
                            obj.dvg_dy(ix, iy, iz+1, :);
                            obj.dvg_dy(ix, iy+1, iz, :);
                            obj.dvg_dy(ix, iy+1, iz+1, :);
                            obj.dvg_dy(ix+1, iy, iz, :);
                            obj.dvg_dy(ix+1, iy, iz+1, :);
                            obj.dvg_dy(ix+1, iy+1, iz, :);
                            obj.dvg_dy(ix+1, iy+1, iz+1, :);
                            % df/dz
                            obj.dvg_dz(ix, iy, iz, :);
                            obj.dvg_dz(ix, iy, iz+1, :);
                            obj.dvg_dz(ix, iy+1, iz, :);
                            obj.dvg_dz(ix, iy+1, iz+1, :);
                            obj.dvg_dz(ix+1, iy, iz, :);
                            obj.dvg_dz(ix+1, iy, iz+1, :);
                            obj.dvg_dz(ix+1, iy+1, iz, :);
                            obj.dvg_dz(ix+1, iy+1, iz+1, :);
                            % d2f/dxdy
                            obj.dvg_dxy(ix, iy, iz, :);
                            obj.dvg_dxy(ix, iy, iz+1, :);
                            obj.dvg_dxy(ix, iy+1, iz, :);
                            obj.dvg_dxy(ix, iy+1, iz+1, :);
                            obj.dvg_dxy(ix+1, iy, iz, :);
                            obj.dvg_dxy(ix+1, iy, iz+1, :);
                            obj.dvg_dxy(ix+1, iy+1, iz, :);
                            obj.dvg_dxy(ix+1, iy+1, iz+1, :);
                            % d2f/dxdz
                            obj.dvg_dxz(ix, iy, iz, :);
                            obj.dvg_dxz(ix, iy, iz+1, :);
                            obj.dvg_dxz(ix, iy+1, iz, :);
                            obj.dvg_dxz(ix, iy+1, iz+1, :);
                            obj.dvg_dxz(ix+1, iy, iz, :);
                            obj.dvg_dxz(ix+1, iy, iz+1, :);
                            obj.dvg_dxz(ix+1, iy+1, iz, :);
                            obj.dvg_dxz(ix+1, iy+1, iz+1, :);
                            % d2f/dydz
                            obj.dvg_dyz(ix, iy, iz, :);
                            obj.dvg_dyz(ix, iy, iz+1, :);
                            obj.dvg_dyz(ix, iy+1, iz, :);
                            obj.dvg_dyz(ix, iy+1, iz+1, :);
                            obj.dvg_dyz(ix+1, iy, iz, :);
                            obj.dvg_dyz(ix+1, iy, iz+1, :);
                            obj.dvg_dyz(ix+1, iy+1, iz, :);
                            obj.dvg_dyz(ix+1, iy+1, iz+1, :);
                            % d3f/dxdydz
                            obj.dvg_dxyz(ix, iy, iz, :);
                            obj.dvg_dxyz(ix, iy, iz+1, :);
                            obj.dvg_dxyz(ix, iy+1, iz, :);
                            obj.dvg_dxyz(ix, iy+1, iz+1, :);
                            obj.dvg_dxyz(ix+1, iy, iz, :);
                            obj.dvg_dxyz(ix+1, iy, iz+1, :);
                            obj.dvg_dxyz(ix+1, iy+1, iz, :);
                            obj.dvg_dxyz(ix+1, iy+1, iz+1, :);
                            ]);

                        %a_sol = obj.M(1:8, :) \ D(1:8,:);
                        a_sol = obj.M \ D;
                        obj.Coefs(ix, iy, iz, :, :) = a_sol;
                    end
                end
            end
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
        
        function field = getFieldAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4, 3]);
            field = [
                tricubic(A_sol(:,:,:,1), xe, ye, ze);
                tricubic(A_sol(:,:,:,2), xe, ye, ze);
                tricubic(A_sol(:,:,:,3), xe, ye, ze);
                ];

        end
        
        function gradient = getGradientAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4, 3]);
            dBx = tricubic_grad(A_sol(:,:,:,1), xe, ye, ze)';
            dBy = tricubic_grad(A_sol(:,:,:,2), xe, ye, ze)';
            dBz = tricubic_grad(A_sol(:,:,:,3), xe, ye, ze)';
            gradient = [dBx/obj.Steps(1); dBy/obj.Steps(2);...
                dBz/obj.Steps(3)];
        end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            normalized = reshape((pm - pmin) ./ obj.Steps, size(positions));
        end
    end
    
end

