classdef TestGrad < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (ClassSetupParameter)
        grid = struct('size', {3, 4, 5}, 'eps', {10., 10., 30.});
        method = {'SimpleRBFInterpolator'};
    end
    
    properties (TestParameter)
        %eps = {10.};
    end
    
    methods (TestClassSetup, ParameterCombination='sequential')
        function ClassSetup(testCase, grid_size, eps)
            NODES_DATASET = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grid_size);
            EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

            % load currents file
            currents = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/currents_3787.h5', '/currents')';

            % load positions
            nodes_pos_fn = fullfile(NODES_DATASET,'/positions.h5');
            xg = h5read(nodes_pos_fn, '/xg');
            yg = h5read(nodes_pos_fn, '/yg');
            zg = h5read(nodes_pos_fn, '/zg');

            xg = permute(xg, [3, 2, 1]);
            yg = permute(yg, [3, 2, 1]);
            zg = permute(zg, [3, 2, 1]);

            xd = xg(:);
            yd = yg(:);
            zd = zg(:);

            nodes = cat(4, xg, yg, zg);

            % load positions
            eval_pos_fn = fullfile(EVAL_DATASET,'/positions.h5');

            xg_ev = h5read(eval_pos_fn, '/xg');
            yg_ev = h5read(eval_pos_fn, '/yg');
            zg_ev = h5read(eval_pos_fn, '/zg');

            xg_ev = permute(xg_ev, [3, 2, 1]);
            yg_ev = permute(yg_ev, [3, 2, 1]);
            zg_ev = permute(zg_ev, [3, 2, 1]);

            positions_ev = cat(4, xg_ev, yg_ev, zg_ev);

            NUM_CURRENTS = 50;
        end
    end
    
    methods (Test)
        function testMae(testCase)
            mae = zeros(NUM_CURRENTS, 3, 3);
            r2 = zeros(NUM_CURRENTS, 3, 3);
            for i=1:NUM_CURRENTS
                fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
                fields = permute(fields, [4, 3, 2, 1]);
                % adding random noise
                fields = fields + 200e-6 * randn(size(fields));

                gradients_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/gradients');
                gradients_ev = permute(gradients_ev, [5, 4, 3, 2, 1]);

                model = SimpleRBFInterpolator(nodes, fields, eps);

                gradients = model.getGradientsAtPositions(positions_ev);

                mae(i,:,:) = gradmae(gradients_ev, gradients);
                r2(i,:,:) = gradr2(gradients_ev, gradients);

                percent_done = 100 * i / NUM_CURRENTS;

                if mod(percent_done, 10) == 0
                    fprintf('%d %% done\n', percent_done);
                end

            end
        end
    end
    
end

