classdef TricubicScalarFieldInterpolator < FieldInterpolator
    %SIMPLETRICUBICINTERPOLATOR Performs tricubic interpolation separately
    % on the three field dimensions
    
    properties (SetAccess = private)
        M
        BFun
        GFun
        VG
    end
    
    methods
        function obj = TricubicScalarFieldInterpolator(nodes, values, M, BFun)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            %[obj.M, obj.BFun, obj.GFun] = get_tricubic_3d_matrix();
            obj.M = M;
            obj.BFun = BFun;
            %obj.GFun = GFun;
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            vg = permute(vg, [3,2,1,4]);
            % we pad the array so we can do i-1 and i+2 on the borders
            obj.VG = padarray(vg, [2, 2, 2, 0], 'circular', 'both');
        end
        
        function [a_sol, x, y, z] = getCoefficients(obj, position)
            pn = getNormalizedPositions(obj,position);
            ix = floor(pn(1)) + 3; 
            iy = floor(pn(2)) + 3; 
            iz = floor(pn(3)) + 3;
            
            x = pn(1) + 3 - ix;
            y = pn(2) + 3 - iy;
            z = pn(3) + 3 - iz;
            
            vg = obj.VG;
            D = zeros(64,1);
            D(9:end) = squeeze([
            % f
                reshape([vg(ix, iy, iz, :);
                vg(ix, iy, iz+1, :);
                vg(ix, iy+1, iz, :);
                vg(ix, iy+1, iz+1, :);
                vg(ix+1, iy, iz, :);
                vg(ix+1, iy, iz+1, :);
                vg(ix+1, iy+1, iz, :);
                vg(ix+1, iy+1, iz+1, :)], 24,1) ;
            % dfx/dy
                0.5 * (vg(ix, iy+1, iz, 1)       - vg(ix, iy-1, iz, 1));
                0.5 * (vg(ix, iy+1, iz+1, 1)     - vg(ix, iy-1, iz+1, 1));
                0.5 * (vg(ix, iy+2, iz, 1)       - vg(ix, iy, iz, 1));
                0.5 * (vg(ix, iy+2, iz+1, 1)     - vg(ix, iy, iz+1, 1));
                0.5 * (vg(ix+1, iy+1, iz, 1)     - vg(ix+1, iy-1, iz, 1));
                0.5 * (vg(ix+1, iy+1, iz+1, 1)   - vg(ix+1, iy-1, iz+1, 1));
                0.5 * (vg(ix+1, iy+2, iz, 1)     - vg(ix+1, iy, iz, 1));
                0.5 * (vg(ix+1, iy+2, iz+1, 1)   - vg(ix+1, iy, iz+1, 1));
            % dfx/dz
                0.5 * (vg(ix, iy, iz+1, 1)       - vg(ix, iy, iz-1, 1));
                0.5 * (vg(ix, iy, iz+2, 1)       - vg(ix, iy, iz, 1));
                0.5 * (vg(ix, iy+1, iz+1, 1)     - vg(ix, iy+1, iz-1, 1));
                0.5 * (vg(ix, iy+1, iz+2, 1)     - vg(ix, iy+1, iz, 1));
                0.5 * (vg(ix+1, iy, iz+1, 1)     - vg(ix+1, iy, iz-1, 1));
                0.5 * (vg(ix+1, iy, iz+2, 1)     - vg(ix+1, iy, iz, 1));
                0.5 * (vg(ix+1, iy+1, iz+1, 1)   - vg(ix+1, iy+1, iz-1, 1));
                0.5 * (vg(ix+1, iy+1, iz+2, 1)   - vg(ix+1, iy+1, iz, 1));
            % dfy/dz
                0.5 * (vg(ix, iy, iz+1, 2)       - vg(ix, iy, iz-1, 2));
                0.5 * (vg(ix, iy, iz+2, 2)       - vg(ix, iy, iz, 2));
                0.5 * (vg(ix, iy+1, iz+1, 2)     - vg(ix, iy+1, iz-1, 2));
                0.5 * (vg(ix, iy+1, iz+2, 2)     - vg(ix, iy+1, iz, 2));
                0.5 * (vg(ix+1, iy, iz+1, 2)     - vg(ix+1, iy, iz-1, 2));
                0.5 * (vg(ix+1, iy, iz+2, 2)     - vg(ix+1, iy, iz, 2));
                0.5 * (vg(ix+1, iy+1, iz+1, 2)   - vg(ix+1, iy+1, iz-1, 2));
                0.5 * (vg(ix+1, iy+1, iz+2, 2)   - vg(ix+1, iy+1, iz, 2));
            % d2fx/dydz
                0.25 * (vg(ix, iy+1, iz+1, 1)    - vg(ix, iy-1, iz+1, 1)    - vg(ix, iy+1, iz-1, 1)        + vg(ix, iy-1, iz-1, 1));
                0.25 * (vg(ix, iy+1, iz+2, 1)    - vg(ix, iy-1, iz+2, 1)    - vg(ix, iy+1, iz, 1)     + vg(ix, iy-1, iz, 1));
                0.25 * (vg(ix, iy+2, iz+1, 1)    - vg(ix, iy, iz+1, 1)      - vg(ix, iy+2, iz-1, 1)        + vg(ix, iy, iz-1, 1));
                0.25 * (vg(ix, iy+2, iz+2, 1)    - vg(ix, iy, iz+2, 1)      - vg(ix, iy+2, iz, 1)          + vg(ix, iy, iz, 1));
                0.25 * (vg(ix+1, iy+1, iz+1, 1)  - vg(ix+1, iy-1, iz+1, 1)  - vg(ix+1, iy+1, iz-1, 1)      + vg(ix+1, iy-1, iz-1, 1));
                0.25 * (vg(ix+1, iy+1, iz+2, 1)  - vg(ix+1, iy-1, iz+2, 1)  - vg(ix+1, iy+1, iz, 1)        + vg(ix+1, iy-1, iz, 1));
                0.25 * (vg(ix+1, iy+2, iz+1, 1)  - vg(ix+1, iy, iz+1, 1)    - vg(ix+1, iy+2, iz-1, 1)      + vg(ix+1, iy, iz-1, 1));
                0.25 * (vg(ix+1, iy+2, iz+2, 1)  - vg(ix+1, iy, iz+2, 1)    - vg(ix+1, iy+2, iz, 1)        + vg(ix+1, iy, iz, 1))
                ]);
             a_sol = obj.M \ D;
        end
        
        function field = getFieldAtPosition(obj, position)
           [a_sol, xe, ye, ze] = obj.getCoefficients(position);
            a_sol = [0; a_sol];
            A_sol = reshape(a_sol, [4 4 4]);
            field = -tricubic_grad_num(A_sol, xe, ye, ze);
        end
        
%         function gradient = getGradientAtPosition(obj, position)
%             [a_sol, x, y, z] = obj.getCoefficients(position);
%             gradient = [obj.GFun(x, y, z, a_sol(:,1)'), ...
%                 obj.GFun(x, y, z, a_sol(:,2)'), ...
%                 obj.GFun(x, y, z, a_sol(:,3)')];
%         end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            steps = (pmax - pmin) ./ (size(obj.NodePositions(:,:,:,1)) - 1);
            normalized = reshape((pm - pmin) ./ steps, size(positions));
        end
    end
    
end

