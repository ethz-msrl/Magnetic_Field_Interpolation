classdef TricubicScalarFieldInterpolator < FieldInterpolator
    %TRICUBICSCALARFIELDINTERPOLATOR Performs tricubic interpolation on
    % a magnetic scalar potential. 
    %   The magnetic field is the negative gradient of the scalar potential
    %   The magnetic field at the voxel corners, as well as dbx/dy dbx/dz
    %   dfy/dz d2fx/dydz are constrained at the corners. the laplacian is
    %   also constrained to be 0 at the voxel corners
    %
    %   Copyright 2020, Samuel Charreyron
    
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
            
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            obj.Steps = (pmax - pmin) ./ (size(obj.NodePositions(:,:,:,1)) - 1);
            
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            obj.VG = permute(vg, [3,2,1,4]);
            
            % here we prescale the interpolation values by the scale
            obj.VG = obj.VG .* permute(repmat(obj.Steps, [obj.Nx,1,obj.Ny, obj.Nz]), ...
                [1,3,4,2]);
            
            [obj.dvg_dy, ~, obj.dvg_dz] = gradient(obj.VG);
            [~, ~, obj.dvg_dyz] = gradient(obj.dvg_dy);
            
            getAllCoefficients(obj);

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

                        a_sol = obj.M \ D;
                        obj.Coefs(ix, iy, iz, :, :) = [0;a_sol];
                    end
                end
            end
        end
        
        function field = getFieldAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            % because the positions are normalized by the step sizes, we
            % need to scale all the derivatives
            field = -tricubic_grad(A_sol, xe, ye, ze) ./ obj.Steps';
        end
        
        function gradient = getGradientAtPosition(obj, position)
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            
            % the magnetic field gradient is the hessian of the scalar
            % potential
            H = -tricubic_hess(A_sol, xe, ye, ze);
            
            % because the positions are normalized by the step sizes, we
            % need to scale all the derivatives
            gradient =  H ./ (obj.Steps'.* obj.Steps);
        end
        
        function gradient = getGradientAtPositionNumeric(obj, position)
            % used this to make sure the analytical gradient is correct
            eps = 1e-6;
            f0 = getFieldAtPosition(obj, position);
            p_ = position + [eps, 0, 0];
            d_dx = (getFieldAtPosition(obj, p_) - f0) / eps;
            p_ = position + [0, eps, 0];
            d_dy = (getFieldAtPosition(obj, p_) - f0) / eps;
            p_ = position + [0, 0, eps];
            d_dz = (getFieldAtPosition(obj, p_) - f0) / eps;
            
            gradient = [d_dx, d_dy, d_dz];
            
        end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            normalized = reshape((pm - pmin) ./ obj.Steps, size(positions));
        end
        
        function field = getFieldAtPostionNumeric(obj, position)
            % used this to make sure the analytical gradient is correct
            eps = 1e-6;
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(position);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            f0 = tricubic(A_sol, xe, ye, ze);
            
            p_ = position + [eps, 0, 0];
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(p_);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            dphi_dx = tricubic(A_sol, xe, ye, ze);
            
            p_ = position + [0, eps, 0];
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(p_);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            dphi_dy = tricubic(A_sol, xe, ye, ze);
            
            p_ = position + [0, 0, eps];
            [ix, iy, iz, xe, ye, ze] = obj.getIndices(p_);
            A_sol = reshape(obj.Coefs(ix, iy, iz, :, :), [4, 4, 4]);
            dphi_dz = tricubic(A_sol, xe, ye, ze);
            
            field =  (-[dphi_dx; dphi_dy; dphi_dz] + f0) / eps;
        end
    end
end

