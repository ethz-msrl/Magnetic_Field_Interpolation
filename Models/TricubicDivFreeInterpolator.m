classdef TricubicDivFreeInterpolator < FieldInterpolator
    %SIMPLETRICUBICINTERPOLATOR Performs tricubic interpolation separately
    % on the three field dimensions
    
    properties (SetAccess = private)
        M
        VG
        Nx
        Ny
        Nz
%         cvg_x
%         cvg_y
%         cvg_z
        gcvg_xx
        gcvg_xy
        gcvg_xz
        gcvg_yx
        gcvg_yy
        gcvg_yz
        gcvg_zx
        gcvg_zy
        gcvg_zz
        CVG
%         dvg_dx
%         dvg_dy
%         dvg_dz
%         dvg_dxy
%         dvg_dxz
%         dvg_dyz
%         dvg_dxyz
        Coefs
        Steps
    end
    
    methods
        function obj = TricubicDivFreeInterpolator(nodes, values, M)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            obj.M = M;
            obj.Nx = size(obj.NodePositions, 1);
            obj.Ny = size(obj.NodePositions, 2);
            obj.Nz = size(obj.NodePositions, 3);
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            obj.VG = permute(vg, [3,2,1,4]);
            
            [dvg_dy, dvg_dx, dvg_dz] = gradient(obj.VG);
            cvg_x = dvg_dy(:,:,:,3) - dvg_dz(:,:,:,2);
            cvg_y = dvg_dz(:,:,:,1) - dvg_dx(:,:,:,3);
            cvg_z = dvg_dx(:,:,:,2) - dvg_dy(:,:,:,1);
            obj.CVG = cat(4, cvg_x, cvg_y, cvg_z);
            
%             [obj.gcvg_xy, obj.gcvg_xx, obj.gcvg_xz] = gradient(obj.cvg_x);
%             [obj.gcvg_yx, obj.gcvg_yy, obj.gcvg_yz] = gradient(obj.cvg_y);
%             [obj.gcvg_zx, obj.gcvg_zy, obj.gcvg_zz] = gradient(obj.cvg_z);
            
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
            obj.Coefs = zeros(obj.Nx-1, obj.Ny-1, obj.Nz-1, 64, 3);
            for ix = 1:(obj.Nx-1)
                for iy = 1:(obj.Ny-1)
                    for iz = 1:(obj.Nz-1)
                        D = squeeze([
                       % B[
                            reshape([
                            obj.CVG(ix, iy, iz, :);
                            obj.CVG(ix, iy, iz+1, :);
                            obj.CVG(ix, iy+1, iz, :);
                            obj.CVG(ix, iy+1, iz+1, :);
                            obj.CVG(ix+1, iy, iz, :);
                            obj.CVG(ix+1, iy, iz+1, :);
                            obj.CVG(ix+1, iy+1, iz, :);
                            obj.CVG(ix+1, iy+1, iz+1, :)], 24,1) ;
                            ]);

                        %a_sol = obj.M(9:32,:) \ D(9:32,:);
                        a_sol = obj.M(1:24,:) \ D;
                        obj.Coefs(ix, iy, iz, :, :) = reshape(a_sol, 64, 3);
                    end
                end
            end
        end
        
        function field = getFieldAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4, 3]);
            field = tricubic_curl(A_sol, xe, ye, ze);
        end
        
%         function gradient = getGradientAtPosition(obj, position)
%             [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
%             A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
%             H = -tricubic_hess(A_sol, xe, ye, ze);
%             gradient = [H(1,:)/obj.Steps(1); ...
%                 H(2,:)/obj.Steps(2); ...
%                 H(3,:)/obj.Steps(3)];
%         end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            normalized = reshape((pm - pmin) ./ obj.Steps, size(positions));
        end
    end
    
end

