clear all;
load_settings;

% MAKE SURE THERE IS A TRAILING SLASH
% we use grid size of 5

nodes_dataset = sprintf(options.base_dataset, options.rbf_eps_grid);

% Grabbing grid number from nodes_dataset
% temp = split(nodes_dataset, '/');
% temp = regexp(temp{end-1},'_(\d)_', 'tokens');
% grid = str2num(temp{1}{1});

if options.recompute ~= 0
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 0, [0.1, 2], options);
    save(sprintf('data/eps_score/RBF-G-3D_%d.mat', grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 1, [0.1, 2], options);
    save(sprintf('data/eps_score/RBF-MQ-3D_%d.mat', grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 2, [0.1, 2], options);
    save(sprintf('data/eps_score/RBF-G-DF_%d.mat', grid), 'eps_v', 'scores', 'cond_numbers');
    
    [eps_v, scores, cond_numbers] = get_eps_score(nodes_dataset, 3, [0.1, 2], options);
    save(sprintf('data/eps_score/RBF-MQ-DF_%d.mat', grid), 'eps_v', 'scores', 'cond_numbers');
end

