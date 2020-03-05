function [eps_v, scores, cond_numbers] = test_eps_score(nodes_dataset, type, bounds, options)
% nodes_dataset: filename of the nodes file to use
% type 0: Gaussian RBF3D
% type 1: Multiquadric 3D
% type 2: Gaussian Div-free
% type 3: Multiquadric Div-free
% bounds: an array of two values containing the lower and upper bound of
% the eps param. It will calculate the score in a linspace array of 20
% elemnents between the bounds.

    % load positions
    nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
    xg = h5read(nodes_pos_fn, '/xg');
    yg = h5read(nodes_pos_fn, '/yg');
    zg = h5read(nodes_pos_fn, '/zg');
    
    % normalize positions'
    [xg, yg, zg, maxp, minp] = normalize_positions_minmax(xg, yg, zg); 

    xg = permute(xg, [3, 2, 1]);
    yg = permute(yg, [3, 2, 1]);
    zg = permute(zg, [3, 2, 1]);

    xd = xg(:);
    yd = yg(:);
    zd = zg(:);

    nodes = cat(4, xg, yg, zg);

    % load positions
    eval_pos_fn = fullfile(options.eval_dataset,'/positions.h5');

    xg_ev = h5read(eval_pos_fn, '/xg');
    yg_ev = h5read(eval_pos_fn, '/yg');
    zg_ev = h5read(eval_pos_fn, '/zg');
    
    % normalize positions
    [xg_ev, yg_ev, zg_ev] = normalize_positions_minmax(xg_ev, yg_ev, zg_ev, maxp, minp);

    xg_ev = permute(xg_ev, [3, 2, 1]);
    yg_ev = permute(yg_ev, [3, 2, 1]);
    zg_ev = permute(zg_ev, [3, 2, 1]);

    positions_ev = cat(4, xg_ev, yg_ev, zg_ev);
    
    % we use a smaller number of currents here or it would take forever
    NUM_CURRENTS = 10;

    eps_v = linspace(bounds(1), bounds(2), 20);

    scores = zeros(length(eps_v), 1);
    cond_numbers = zeros(length(eps_v), 1);
    for i = 1:length(eps_v)
        eps = eps_v(i);
        nrmse_scores = zeros(NUM_CURRENTS, 3);
        cond_numbers_temp = zeros(NUM_CURRENTS,1);
        for j=1:NUM_CURRENTS
            fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', j)), '/fields');
            fields = permute(fields, [4, 3, 2, 1]);
            
            % adding random noise
            fields = fields + options.noise_std * randn(size(fields));

            fields_ev = h5read(fullfile(options.eval_dataset, 'v', sprintf('%04d.h5', j)), '/fields');
            fields_ev = permute(fields_ev, [4, 3, 2, 1]);

            if type == 0
                model = RBF3DGaussianInterpolator(nodes, fields, eps);
            elseif type == 1
                model = RBF3DMultiquadricInterpolator(nodes, fields, eps);
            elseif type == 2
                model = RBFDivFreeGaussianInterpolator(nodes, fields, eps);
            else
                model = RBFDivFreeMultiquadricInterpolator(nodes, fields, eps);
            end
            cond_numbers_temp(j) = model.CondNumber;
            
            ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);

            nrmse_scores(j,:) = 100*ev.get_nrmse();
        end
        
        cond_numbers(i) = mean(cond_numbers_temp);
        scores(i,:) = mean(nrmse_scores(:));
    end

end
