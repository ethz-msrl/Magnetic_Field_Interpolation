classdef SimpleTricubicInterpolator < FieldInterpolator
    %SIMPLETRICUBICINTERPOLATOR Performs tricubic interpolation separately
    % on the three field dimensions
    
    properties (SetAccess = private)
        M
        BFun
        GFun
        VG
    end
    
    methods
        function obj = SimpleTricubicInterpolator(nodes, values, M, BFun , GFun)
            obj.NodePositions = nodes;
            obj.NodeValues = values;
            %[obj.M, obj.BFun, obj.GFun] = get_tricubic_3d_matrix();
            obj.M = M;
            obj.BFun = BFun;
            obj.GFun = GFun;
            % we need to convert data axes from Z,Y,X ngrid to X,Y,Z
            vg = obj.NodeValues;
            vg = permute(vg, [3,2,1,4]);
            % we pad the array so we can do i-1 and i+2 on the borders
            obj.VG = padarray(vg, [2, 2, 2, 0], 'circular', 'both');
        end
        
        function field = getFieldAtPosition(obj, position)
            pn = getNormalizedPositions(obj,position);
            ix = floor(pn(1)) + 3; 
            iy = floor(pn(2)) + 3; 
            iz = floor(pn(3)) + 3;
            
            x = pn(1) + 3 - ix;
            y = pn(2) + 3 - iy;
            z = pn(3) + 3 - iz;
            
            vg = obj.VG;

            D = squeeze([
            % f
                vg(ix, iy, iz, :);
                vg(ix, iy, iz+1, :);
                vg(ix, iy+1, iz, :);
                vg(ix, iy+1, iz+1, :);
                vg(ix+1, iy, iz, :);
                vg(ix+1, iy, iz+1, :);
                vg(ix+1, iy+1, iz, :);
                vg(ix+1, iy+1, iz+1, :);
            % df/dx
                0.5 * (vg(ix+1, iy, iz, :)       - vg(ix-1, iy, iz, :));
                0.5 * (vg(ix+1, iy, iz+1, :)     - vg(ix-1, iy, iz+1, :));
                0.5 * (vg(ix+1, iy+1, iz, :)     - vg(ix-1, iy+1, iz, :));
                0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix-1, iy+1, iz+1, :));
                0.5 * (vg(ix+2, iy, iz, :)       - vg(ix, iy, iz, :));
                0.5 * (vg(ix+2, iy, iz+1, :)     - vg(ix, iy, iz+1, :));
                0.5 * (vg(ix+2, iy+1, iz, :)     - vg(ix, iy+1, iz, :));
                0.5 * (vg(ix+2, iy+1, iz+1, :)   - vg(ix, iy+1, iz+1, :));
            % df/dy
                0.5 * (vg(ix, iy+1, iz, :)       - vg(ix, iy-1, iz, :));
                0.5 * (vg(ix, iy+1, iz+1, :)     - vg(ix, iy-1, iz+1, :));
                0.5 * (vg(ix, iy+2, iz, :)       - vg(ix, iy, iz, :));
                0.5 * (vg(ix, iy+2, iz+1, :)     - vg(ix, iy, iz+1, :));
                0.5 * (vg(ix+1, iy+1, iz, :)     - vg(ix+1, iy-1, iz, :));
                0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix+1, iy-1, iz+1, :));
                0.5 * (vg(ix+1, iy+2, iz, :)     - vg(ix+1, iy, iz, :));
                0.5 * (vg(ix+1, iy+2, iz+1, :)   - vg(ix+1, iy, iz+1, :));
            % df/dz
                0.5 * (vg(ix, iy, iz+1, :)       - vg(ix, iy, iz-1, :));
                0.5 * (vg(ix, iy, iz+2, :)       - vg(ix, iy, iz, :));
                0.5 * (vg(ix, iy+1, iz+1, :)     - vg(ix, iy+1, iz-1, :));
                0.5 * (vg(ix, iy+1, iz+2, :)     - vg(ix, iy+1, iz, :));
                0.5 * (vg(ix+1, iy, iz+1, :)     - vg(ix+1, iy, iz-1, :));
                0.5 * (vg(ix+1, iy, iz+2, :)     - vg(ix+1, iy, iz, :));
                0.5 * (vg(ix+1, iy+1, iz+1, :)   - vg(ix+1, iy+1, iz-1, :));
                0.5 * (vg(ix+1, iy+1, iz+2, :)   - vg(ix+1, iy+1, iz, :));
            % d2f/dxdy
                0.25 * (vg(ix+1, iy+1, iz, :)    - vg(ix-1, iy+1, iz, :)    - vg(ix+1, iy - 1, iz)      + vg(ix-1, iy-1, iz));
                0.25 * (vg(ix+1, iy+1, iz+1, :)  - vg(ix-1, iy+1, iz+1, :)  - vg(ix+1, iy - 1, iz+1)    + vg(ix-1, iy-1, iz+1));
                0.25 * (vg(ix+1, iy+2, iz, :)    - vg(ix-1, iy+2, iz, :)    - vg(ix+1, iy, iz)          + vg(ix-1, iy, iz));
                0.25 * (vg(ix+1, iy+2, iz+1, :)  - vg(ix-1, iy+2, iz+1, :)  - vg(ix+1, iy, iz+1)        + vg(ix-1, iy, iz+1));
                0.25 * (vg(ix+2, iy+1, iz, :)    - vg(ix, iy+1, iz, :)      - vg(ix+2, iy - 1, iz)      + vg(ix, iy-1, iz));
                0.25 * (vg(ix+2, iy+1, iz+1, :)  - vg(ix, iy+1, iz+1, :)    - vg(ix+2, iy - 1, iz+1)    + vg(ix, iy-1, iz+1));
                0.25 * (vg(ix+2, iy+2, iz, :)    - vg(ix, iy+2, iz, :)      - vg(ix+2, iy, iz)          + vg(ix, iy, iz));
                0.25 * (vg(ix+2, iy+2, iz+1, :)  - vg(ix, iy+2, iz+1, :)    - vg(ix+2, iy, iz+1)        + vg(ix, iy, iz+1));
            % d2f/dxdz
                0.25 * (vg(ix+1, iy, iz+1, :)    - vg(ix-1, iy, iz+1, :)    - vg(ix+1, iy - 1, iz-1)    + vg(ix-1, iy, iz-1));
                0.25 * (vg(ix+1, iy, iz+2, :)    - vg(ix-1, iy, iz+2, :)    - vg(ix+1, iy - 1, iz)      + vg(ix-1, iy, iz));
                0.25 * (vg(ix+1, iy+1, iz+1, :)  - vg(ix-1, iy+1, iz+1, :)  - vg(ix+1, iy, iz-1)        + vg(ix-1, iy+1, iz-1));
                0.25 * (vg(ix+1, iy+1, iz+2, :)  - vg(ix-1, iy+1, iz+2, :)  - vg(ix+1, iy, iz)          + vg(ix-1, iy+1, iz));
                0.25 * (vg(ix+2, iy, iz+1, :)    - vg(ix, iy, iz+1, :)      - vg(ix+2, iy - 1, iz-1)    + vg(ix, iy, iz-1));
                0.25 * (vg(ix+2, iy, iz+2, :)    - vg(ix, iy, iz+2, :)      - vg(ix+2, iy - 1, iz)      + vg(ix, iy, iz));
                0.25 * (vg(ix+2, iy+1, iz+1, :)  - vg(ix, iy+1, iz+1, :)    - vg(ix+2, iy, iz-1)        + vg(ix, iy+1, iz-1));
                0.25 * (vg(ix+2, iy+1, iz+2, :)  - vg(ix, iy+1, iz+2, :)    - vg(ix+2, iy, iz)          + vg(ix, iy+1, iz));
            % d2f/dydz
                0.25 * (vg(ix, iy+1, iz+1, :)    - vg(ix, iy-1, iz+1, :)    - vg(ix, iy+1, iz-1)        + vg(ix, iy-1, iz-1));
                0.25 * (vg(ix, iy+1, iz+2, :)    - vg(ix, iy-1, iz+2, :)    - vg(ix, iy+1, iz)     + vg(ix, iy-1, iz));
                0.25 * (vg(ix, iy+2, iz+1, :)    - vg(ix, iy, iz+1, :)      - vg(ix, iy+2, iz-1)        + vg(ix, iy, iz-1));
                0.25 * (vg(ix, iy+2, iz+2, :)    - vg(ix, iy, iz+2, :)      - vg(ix, iy+2, iz)          + vg(ix, iy, iz));
                0.25 * (vg(ix+1, iy+1, iz+1, :)  - vg(ix+1, iy-1, iz+1, :)  - vg(ix+1, iy+1, iz-1)      + vg(ix+1, iy-1, iz-1));
                0.25 * (vg(ix+1, iy+1, iz+2, :)  - vg(ix+1, iy-1, iz+2, :)  - vg(ix+1, iy+1, iz)        + vg(ix+1, iy-1, iz));
                0.25 * (vg(ix+1, iy+2, iz+1, :)  - vg(ix+1, iy, iz+1, :)    - vg(ix+1, iy+2, iz-1)      + vg(ix+1, iy, iz-1));
                0.25 * (vg(ix+1, iy+2, iz+2, :)  - vg(ix+1, iy, iz+2, :)    - vg(ix+1, iy+2, iz)        + vg(ix+1, iy, iz));
            % d3f/dxdydz
                0.125 * (vg(ix+1, iy+1, iz+1, :) - vg(ix-1, iy+1, iz+1, :)  - vg(ix+1, iy-1, iz+1)      + vg(ix-1, iy-1, iz+1, :)   - vg(ix+1, iy+1, iz-1)  + vg(ix-1, iy+1, iz-1) + vg(ix+1, iy-1, iz-1)   -  vg(ix-1, iy-1, iz-1));
                0.125 * (vg(ix+1, iy+1, iz+2, :) - vg(ix-1, iy+1, iz+2, :)  - vg(ix+1, iy-1, iz+2)      + vg(ix-1, iy-1, iz+2, :)   - vg(ix+1, iy+1, iz)    + vg(ix-1, iy+1, iz) + vg(ix+1, iy-1, iz)       -  vg(ix-1, iy-1, iz));
                0.125 * (vg(ix+1, iy+2, iz+1, :) - vg(ix-1, iy+2, iz+1, :)  - vg(ix+1, iy, iz+1)        + vg(ix-1, iy, iz+1, :)     - vg(ix+1, iy+2, iz-1)  + vg(ix-1, iy+2, iz-1) + vg(ix+1, iy, iz-1)     -  vg(ix-1, iy, iz-1));
                0.125 * (vg(ix+1, iy+2, iz+2, :) - vg(ix-1, iy+2, iz+2, :)  - vg(ix+1, iy, iz+2)        + vg(ix-1, iy, iz+2, :)     - vg(ix+1, iy+2, iz)    + vg(ix-1, iy+2, iz) + vg(ix+1, iy, iz)         -  vg(ix-1, iy, iz));
                0.125 * (vg(ix+2, iy+1, iz+1, :) - vg(ix,   iy+1, iz+1, :)  - vg(ix+2, iy-1, iz+1)      + vg(ix, iy-1, iz+1, :)     - vg(ix+2, iy+1, iz-1)  + vg(ix,   iy+1, iz-1) + vg(ix+2, iy-1, iz-1)   -  vg(ix, iy-1, iz-1));
                0.125 * (vg(ix+2, iy+1, iz+2, :) - vg(ix,   iy+1, iz+2, :)  - vg(ix+2, iy-1, iz+2)      + vg(ix, iy-1, iz+2, :)     - vg(ix+2, iy+1, iz)    + vg(ix,   iy+1, iz) + vg(ix+2, iy-1, iz)       -  vg(ix, iy-1, iz));
                0.125 * (vg(ix+2, iy+2, iz+1, :) - vg(ix,   iy+2, iz+1, :)  - vg(ix+2, iy, iz+1)        + vg(ix, iy, iz+1, :)       - vg(ix+2, iy+2, iz-1)  + vg(ix,   iy+2, iz-1) + vg(ix+2, iy, iz-1)     -  vg(ix, iy, iz-1));
                0.125 * (vg(ix+2, iy+2, iz+2, :) - vg(ix,   iy+2, iz+2, :)  - vg(ix+2, iy, iz+2)        + vg(ix, iy, iz+2, :)       - vg(ix+2, iy+2, iz)    + vg(ix,   iy+2, iz) + vg(ix+2, iy, iz)         -  vg(ix, iy, iz));
                ]);
             a_sol = obj.M \ D;
             field = [  obj.BFun(x, y, z, a_sol(:,1)'); ...
                    obj.BFun(x, y, z, a_sol(:,2)'); ...
                    obj.BFun(x, y, z, a_sol(:,3)')];
        end
        
        function normalized = getNormalizedPositions(obj, positions)
            pm = reshape(positions, [], 3);
            pmin = min(reshape(obj.NodePositions, [],3),[],1);
            pmax = max(reshape(obj.NodePositions, [],3),[],1);
            steps = (pmax - pmin) ./ (size(obj.NodePositions(:,:,:,1)) - 1);
            normalized = reshape((pm - pmin) ./ steps, size(positions));
        end
    end
    
end

