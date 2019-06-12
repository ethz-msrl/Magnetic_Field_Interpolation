function test_fields( model_name, grid, noise_std)
    
    results = struct('grid_size', {}, 'mae', {}, 'r2', {}, 'mean_div', {}, ...
        'mean_curl', {});
    for g=1:length(grid)

        grid_size = grid(g).size;
        
        fprintf('running grid size: %dx%d\n', grid_size, grid_size);

        nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grid_size);
        EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

        % load positions
        nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
        xg = h5read(nodes_pos_fn, '/xg');
        yg = h5read(nodes_pos_fn, '/yg');
        zg = h5read(nodes_pos_fn, '/zg');

        xg = permute(xg, [3, 2, 1]);
        yg = permute(yg, [3, 2, 1]);
        zg = permute(zg, [3, 2, 1]);

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

        NUM_CURRENTS = 20;

        mae = zeros(NUM_CURRENTS, 3);
        rmse = zeros(NUM_CURRENTS,3);
        r2 = zeros(NUM_CURRENTS, 3);

        for i=1:NUM_CURRENTS
            fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
            fields = permute(fields, [4, 3, 2, 1]);
            % adding random noise
            fields = fields + noise_std * randn(size(fields));
            
            fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
            fields_ev = permute(fields_ev, [4, 3, 2, 1]);

            if strcmp(model_name, 'RBF3D')
                eps = grid(g).eps;
                model = SimpleRBFInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'RBFD')
                eps = grid(g).eps;
                model = DivFreeRBFInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'SPL3D')
                degree = grid(g).degree;
                model = BSpline3DInterpolator(nodes, fields, degree);
            elseif strcmp(model_name, 'SPLLPL')
                degree = grid(g).degree;
                model = BSplineLaplacianInterpolator(nodes, fields, degree);
            elseif strcmp(model_name, 'TRI3D')
                load('tricubic_simple_M.mat');
                model = SimpleTricubicInterpolator(nodes, fields, M);
            elseif strcmp(model_name, 'TRILPL')
                load('tricubic_scalar_field_matrix.mat');
                model = TricubicScalarFieldInterpolator(nodes, fields, M);
            else
                error('invalid model name');
            end
            
            ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);

            mae(i,:) = ev.get_mae();
            r2(i,:) = ev.get_r2();

            percent_done = 100 * i / NUM_CURRENTS;

            if mod(percent_done, 10) == 0
                fprintf('%d %% done\n', percent_done);
            end

        end
        
        results(g).grid_size = grid_size;
        results(g).mae = reshape(squeeze(mean(mae,1)), [], 1);
        results(g).r2 = reshape(squeeze(mean(r2,1)), [], 1);
        results(g).rmse = reshape(squeeze(mean(r2,1)), [], 1);

        disp('');
        fprintf('grid size: %d\n', grid_size);
        disp('MAE (mT):')
        disp(1000*squeeze(mean(mae, 1)));

        disp('R2:')
        disp(squeeze(mean(r2, 1)));

        fprintf('\n\n');

    end
    
    save_fn = sprintf('data/fields/%s_%d.mat', model_name, noise_std);
    save(save_fn, 'results');
end

