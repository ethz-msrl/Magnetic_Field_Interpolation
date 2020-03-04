function test_fields( model_name, grid, options)
% TEST_FIELDS Generates test data for fields and saves to a file
%   model_name: a string defining the interpolation model type
%   grid: a struct containing both the grid size and params that are needed
%   by the model
%   options: a struct containing the global options. See settings.m
%   
%   a file will be saved to data/fields/{model_name}_{noise_std}.mat
    
    % the results is a struct array with the calculated metrics for each
    % grid size
    results = struct('grid_size', {}, 'mae', {}, 'r2', {}, 'mean_div', {}, ...
        'mean_curl', {});
    for g=1:length(grid)

        grid_size = grid(g).size;
        
        fprintf('running grid size: %dx%d\n', grid_size, grid_size);

        nodes_dataset = sprintf(options.base_dataset, grid_size);

        % the are the nodes that are used by the interpolant
        nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
        xg = h5read(nodes_pos_fn, '/xg');
        yg = h5read(nodes_pos_fn, '/yg');
        zg = h5read(nodes_pos_fn, '/zg');
        
        % the positions are in the z,y,x format because of deep-fluids
        xg = permute(xg, [3, 2, 1]);
        yg = permute(yg, [3, 2, 1]);
        zg = permute(zg, [3, 2, 1]);

        nodes = cat(4, xg, yg, zg);

        % these are the positions at which the interpolation is evaluated
        % usually it's a fine grid with Ng=16
        eval_pos_fn = fullfile(options.eval_dataset,'/positions.h5');

        xg_ev = h5read(eval_pos_fn, '/xg');
        yg_ev = h5read(eval_pos_fn, '/yg');
        zg_ev = h5read(eval_pos_fn, '/zg');
        
        xg_ev = permute(xg_ev, [3, 2, 1]);
        yg_ev = permute(yg_ev, [3, 2, 1]);
        zg_ev = permute(zg_ev, [3, 2, 1]);
        
        % the positions are in the z,y,x format because of deep-fluids
        positions_ev = cat(4, xg_ev, yg_ev, zg_ev);
        
        % Mean absolute error
        mae = zeros(options.num_currents, 3);
        
        % Mean absolute error normalized (to range of mae)
        nmae = zeros(options.num_currents, 3);
        
        % root mean square error
        rmse = zeros(options.num_currents,3);
        
        % RMSE normalized to the range
        nrmse = zeros(options.num_currents,3);
        
        % coefficient of determination
        r2 = zeros(options.num_currents, 3);

        for i=1:options.num_currents
            fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
            fields = permute(fields, [4, 3, 2, 1]);
            % adding random noise
            % We use this even for experiments without noise where the
            % noise_std is 0
            fields = fields + options.noise_std * randn(size(fields));
            
            % these are the field values on the evaluation grid
            % note that they are obtained from synthetic data generated
            % with the MPEM so they are naturally curl and divergence free
            fields_ev = h5read(fullfile(options.eval_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
            fields_ev = permute(fields_ev, [4, 3, 2, 1]);

            if strcmp(model_name, 'RBF-G-3D')
                % the shape parameter of the RBF
                eps = grid(g).eps;
                model = RBF3DInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'RBF-MQ-3D')
                % the shape parameter of the RBF
                eps = grid(g).eps;
                model = RBF3DMultiquadricInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'RBF-G-DF')
                % the shape parameter of the RBF
                eps = grid(g).eps;
                model = RBFDivFreeInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'RBF-MQ-DF')
                % the shape parameter of the RBF
                eps = grid(g).eps;
                model = RBFDivFreeMultiquadricInterpolator(nodes, fields, eps);
            elseif strcmp(model_name, 'SPL-3D')
                 % this is the degree of the Bspline (d=4 for cubic)
                degree = grid(g).degree;
                model = BSpline3DInterpolator(nodes, fields, degree);
            elseif strcmp(model_name, 'SPL-LPL')
                 % this is the degree of the Bspline (d=4 for cubic)
                degree = grid(g).degree;
                model = BSplineLaplacianInterpolator(nodes, fields, degree);
            elseif strcmp(model_name, 'TRI-3D')
                load('tricubic_3D_M.mat');
                model = Tricubic3DInterpolator(nodes, fields, M);
            elseif strcmp(model_name, 'TRI-LPL')
                load('tricubic_scalar_field_M.mat');
                model = TricubicScalarFieldInterpolator(nodes, fields, M);
            elseif strcmp(model_name, 'TRI-DF')
                load('tricubic_divfree_M.mat');
                model = TricubicDivFreeInterpolator(nodes, fields, M);
            else
                error('invalid model name');
            end
            
            ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
            
            mae(i,:) = ev.get_mae();
            nmae(i,:) = ev.get_nmae();
            rmse(i,:) = ev.get_rmse();
            nrmse(i,:) = ev.get_nrmse();
            r2(i,:) = ev.get_r2();
            
            if options.print_progress
                percent_done = 100 * i / options.num_currents;
                
                if mod(percent_done, 10) == 0
                    fprintf('%d %% done\n', percent_done);
                end
            end

        end
        
        results(g).grid_size = grid_size;
        results(g).mae = reshape(squeeze(mean(mae,1)), [], 1);
        results(g).nmae = reshape(squeeze(mean(nmae,1)), [], 1);
        results(g).r2 = reshape(squeeze(mean(r2,1)), [], 1);
        results(g).rmse = reshape(squeeze(mean(rmse,1)), [], 1);
        results(g).nrmse = reshape(squeeze(mean(nrmse,1)), [], 1);
        
        if options.print_progress
            disp('');
            fprintf('grid size: %d\n', grid_size);
            disp('MAE (%):')
            disp(100*squeeze(mean(nmae, 1)));
            
            disp('R2:')
            disp(squeeze(mean(r2, 1)));
            
            fprintf('\n\n');
        end

    end
    
    save_fn = sprintf('%s/fields/%s_%d.mat', options.data_base_path, ...
        model_name, options.noise_std);
    save(save_fn, 'results');
end

