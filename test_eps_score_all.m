% calculates the N-RMSE fit for each RBF method over a
% linear range of shape parameter values and saves the data for use later
%
% The bounds are set to be between 0.2 and 2 for all methods, because we
% are using positions that are normalized to be between 0 and 1.
%
%   Copyright 2020, ETH Zurich, Multi Scale Robotics Lab, Samuel Charreyron

clear all;
load_settings;

nodes_dataset = sprintf(options.base_dataset, options.rbf_eps_grid);

% Grabbing grid number from nodes_dataset

if options.recompute ~= 0
    [eps_v, scores, cond_numbers] = test_eps_score(nodes_dataset, 0, [0.1, 2], options);
    save(sprintf('%s/eps_score/RBF-G-3D_%d.mat', options.data_base_path, options.rbf_eps_grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = test_eps_score(nodes_dataset, 1, [0.1, 2], options);
    save(sprintf('%s/eps_score/RBF-MQ-3D_%d.mat', options.data_base_path, options.rbf_eps_grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = test_eps_score(nodes_dataset, 2, [0.1, 2], options);
    save(sprintf('%s/eps_score/RBF-G-DF_%d.mat', options.data_base_path, options.rbf_eps_grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = test_eps_score(nodes_dataset, 3, [0.1, 2], options);
    save(sprintf('%s/eps_score/RBF-MQ-DF_%d.mat', options.data_base_path, options.rbf_eps_grid), 'eps_v', 'scores', 'cond_numbers');
end

