function plot_eps_score(nodes_dataset)

    %nodes_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/';
    EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

    % load currents file
    currents = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/currents_3787.h5', '/currents')';

    % load positions
    nodes_pos_fn = fullfile(nodes_dataset,'/positions.h5');
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

    NUM_CURRENTS = 10;

    eps_v = linspace(5,30,20);

    scores = zeros(length(eps_v), 1);
    for i = 1:length(eps_v)
        eps = eps_v(i);
        nrmse_scores = zeros(NUM_CURRENTS, 3);
        for j=1:NUM_CURRENTS
            fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', j)), '/fields');
            fields = permute(fields, [4, 3, 2, 1]);

            fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', j)), '/fields');
            fields_ev = permute(fields_ev, [4, 3, 2, 1]);

            model = RBF3DInterpolator(nodes, fields, eps);
            ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);

            nrmse_scores(j,:) = 100*ev.get_nrmse();
        end

        scores(i,:) = mean(nrmse_scores(:));
    end

    figure();
    bar(eps_v, scores);
    xlabel('eps');
    ylabel('normalized RMSE (%)');

end
