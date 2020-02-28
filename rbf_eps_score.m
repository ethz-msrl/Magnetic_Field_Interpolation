function y = rbf_eps_score(nodes_dataset, eps, type)
% type 0: Gaussian RBF3D
% type 1: Multiquadric 3D
% type 2: Gaussian Div-free
% type 3: Multiquadric Div-free

    %nodes_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/';
    EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

    % load currents file
    currents = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/currents_3787.h5', '/currents')';

    % load positions
    nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
    xg = h5read(nodes_pos_fn, '/xg');
    yg = h5read(nodes_pos_fn, '/yg');
    zg = h5read(nodes_pos_fn, '/zg');
    
    [xg, yg, zg, maxp, minp] = normalize_positions_minmax(xg, yg, zg);

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
    
    % normalize positions
    [xg_ev, yg_ev, zg_ev] = normalize_positions_minmax(xg_ev, yg_ev, zg_ev, maxp, minp);

    xg_ev = permute(xg_ev, [3, 2, 1]);
    yg_ev = permute(yg_ev, [3, 2, 1]);
    zg_ev = permute(zg_ev, [3, 2, 1]);

    positions_ev = cat(4, xg_ev, yg_ev, zg_ev);

    NUM_CURRENTS = 5;

    nrmse_scores = zeros(NUM_CURRENTS, 3);
    for j=1:NUM_CURRENTS
        fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', j)), '/fields');
        fields = permute(fields, [4, 3, 2, 1]);

        fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', j)), '/fields');
        fields_ev = permute(fields_ev, [4, 3, 2, 1]);
        
        if type == 0
            model = RBF3DInterpolator(nodes, fields, eps);
        elseif type == 1
            model = RBF3DMultiquadricInterpolator(nodes, fields, eps);
        elseif type == 2
            model = RBFDivFreeInterpolator(nodes, fields, eps);
        else
            model = RBFDivFreeMultiquadricInterpolator(nodes, fields, eps);
        end
        
        ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);

        nrmse_scores(j,:) = 100*ev.get_nrmse();
    end

    y = mean(nrmse_scores(:));
end
