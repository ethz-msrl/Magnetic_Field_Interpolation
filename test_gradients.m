function test_gradients( model_name, grid, options)
% TEST_GRADIENTS Generates test data for gradients and saves to a file
%   model_name: a string defining the interpolation model type
%   grid: a struct containing both the grid size and params that are needed
%   by the model
%   options: a struct containing the global options. See settings.m
%
%   a file will be saved to data/gradients/{model_name}_{noise_std}.mat
%   The results is a struct array with the calculated metrics for each
%   grid size and is saved in data/fields for use later
%
%   Copyright 2020, Samuel Charreyron
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
    mae = zeros(options.num_currents, 3, 3);
    
    % Mean absolute error normalized (to range of mae)
    nmae = zeros(options.num_currents, 3, 3);
    
    % root mean square error
    rmse = zeros(options.num_currents, 3, 3);
    
    % RMSE normalized to the range
    nrmse = zeros(options.num_currents, 3, 3);
    
    % coefficient of determination
    r2 = zeros(options.num_currents, 3, 3);
    
    % mean value of the absolute value of the divergence
    mean_div = zeros(options.num_currents, 1);
    
    % mean value of the curl vector norm
    mean_curl = zeros(options.num_currents, 1);
    
    for i=1:options.num_currents
        fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
        fields = permute(fields, [4, 3, 2, 1]);
        % adding random noise
        % We use this even for experiments without noise where the
        % noise_std is 0
        fields = fields + options.noise_std * randn(size(fields));
        
        % these are the gradient values on the evaluation grid
        % note that they are obtained from synthetic data generated
        % with the MPEM so they are naturally curl and divergence free
        gradients_ev = h5read(fullfile(options.eval_dataset, 'v', sprintf('%04d.h5', i)), '/gradients');
        gradients_ev = permute(gradients_ev, [5, 4, 3, 2, 1]);
        
        if strcmp(model_name, 'RBF-G-3D')
            % the shape parameter of the RBF
            eps = grid(g).eps;
            model = RBF3DGaussianInterpolator(nodes, fields, eps);
        elseif strcmp(model_name, 'RBF-MQ-3D')
            % the shape parameter of the RBF
            eps = grid(g).eps;
            model = RBF3DMultiquadricInterpolator(nodes, fields, eps);
        elseif strcmp(model_name, 'RBF-G-DF')
            % the shape parameter of the RBF
            eps = grid(g).eps;
            model = RBFDivFreeGaussianInterpolator(nodes, fields, eps);
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
        else
            error('invalid model name');
        end
        
        gradients = model.getGradientsAtPositions(positions_ev);
        grads = reshape(gradients,[],9);
        % we take the mean absolute divergence since positives and
        % negatives can cancel each other out in the mean
        divs = abs(sum(grads(:,1:3+1:9),2));
        % we take the magnitude of the curl vector
        curls = vecnorm([grads(:,6) - grads(:,8), grads(:,7) - grads(:,3), ...
            grads(:,2) - grads(:,4)], 2, 2);
        mean_curl(i) = mean(curls, 1);
        mean_div(i) = mean(divs(:));
        
        mae(i,:,:) = gradmae(gradients_ev, gradients);
        nmae(i,:,:) = mae(i,:,:) ./ (max(reshape(abs(gradients_ev), [], 3,3),[], 1) - ...
            min(reshape(abs(gradients_ev), [], 3,3),[], 1));
        rmse(i,:,:) = gradrmse(gradients_ev, gradients);
        nrmse(i,:,:) = rmse(i,:,:) ./ (max(reshape(abs(gradients_ev), [], 3, 3), [], 1) - ...
            min(reshape(abs(gradients_ev), [], 3, 3), [], 1));
        r2(i,:,:) = gradr2(gradients_ev, gradients);
        
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
    results(g).rmse = reshape(squeeze(mean(rmse,1)), [], 1);
    results(g).nrmse = reshape(squeeze(mean(nrmse,1)), [], 1);
    results(g).r2 = reshape(squeeze(mean(r2,1)), [], 1);
    
    results(g).mean_div = mean(mean_div,1);
    results(g).mean_curl = mean(mean_curl,1);
    
    if options.print_progress
        disp('');
        fprintf('grid size: %d\n', grid_size);
        disp('MAE (T/m):')
        disp(squeeze(mean(mae, 1)));
        
        disp('R2:')
        disp(squeeze(mean(r2, 1)));
        
        fprintf('Mean divergence abs: %2.3f (mT/m) \n', 1000*mean(mean_div));
        temp = 1000*mean(mean_curl,1);
        fprintf('Mean curl magnitude: %2.3f (mT/m) \n',  temp);
        
        fprintf('\n\n');
    end
    
end

save_fn = sprintf('%s/gradients/%s_%d.mat', options.data_base_path,...
    model_name, options.noise_std);
save(save_fn, 'results');
end

