clear all;

grid = struct('size', {3,4,5,6}, 'eps', {35, 35, 50, 50});

for g=1:length(grid)
    
    grid_size = grid(g).size;
    eps = grid(g).eps;
    fprintf('running grid size: %dx%d, eps: %2.f\n', grid_size, grid_size, eps);

    
    nodes_dataset = sprintf('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_%d_h5/', grid_size);
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

    NUM_CURRENTS = 50;

    mae = zeros(NUM_CURRENTS, 3, 3);
    r2 = zeros(NUM_CURRENTS, 3, 3);
    meandivs = zeros(NUM_CURRENTS, 1);
    
    for i=1:NUM_CURRENTS
        fields = h5read(fullfile(nodes_dataset, 'v', sprintf('%04d.h5', i)), '/fields');
        fields = permute(fields, [4, 3, 2, 1]);
        % adding random noise
        %fields = fields + 200e-6 * randn(size(fields));

        gradients_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/gradients');
        gradients_ev = permute(gradients_ev, [5, 4, 3, 2, 1]);

        model = DivFreeRBFInterpolator(nodes, fields, eps);

        gradients = model.getGradientsAtPositions(positions_ev);
        divs = sum(gradients(:,:,:,1:3+1:9));
        grads = reshape(gradients,[],9);
        curls = [grads(:,6) - grads(:,8), grads(:,7) - grads(:,3), grads(:,2) - grads(:,4)];
        meancurls(i,:) = mean(curls, 1);
        meandivs(i) = mean(divs(:));

        mae(i,:,:) = gradmae(gradients_ev, gradients);
        r2(i,:,:) = gradr2(gradients_ev, gradients);

        percent_done = 100 * i / NUM_CURRENTS;

        if mod(percent_done, 10) == 0
            fprintf('%d %% done\n', percent_done);
        end

    end
    
    save_fn = sprintf('data/divfree_rbf_%dx%d_grad_200mTnoise.mat', grid_size, grid_size);
    save(save_fn, 'mae', 'r2', 'meandivs', 'meancurls');

    disp('');
    fprintf('grid size: %d\n', grid_size);
    disp('MAE (T/m):')
    disp(squeeze(mean(mae, 1)));

    disp('R2:')
    disp(squeeze(mean(r2, 1)));
    
    fprintf('Mean divergence: %2.3f (mT/m) \n', 1000*mean(meandivs));
    temp = 1000*mean(meancurls,1);
    fprintf('Mean curl: [%2.3f, %2.3f, %2.3f] (mT/m) \n',  temp(1), temp(2), temp(3));

    fprintf('\n\n');

end
