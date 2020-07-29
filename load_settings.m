options = struct();

% the number of currents that are picked randomly from the dataset
% to form the evaluation set
options.num_currents = 100;

% the standard deviation of some normal noise to add to the true references
% leave at 0 for no noise added
options.noise_std = 0;

% path to the dataset that we use for evaluation
options.eval_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_synthetic_mpem_v1/mpem_synthetic_16_h5/';

% path to the base dataset that we use for the interpolant
% this string will be formatted later to replace in the right grid size
options.base_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_synthetic_mpem_v1/mpem_synthetic_%d_h5/';

% path to the currents
options.currents_dataset = '/Volumes/msrl/users/samuelch/datasets/cmag_synthetic_mpem_v1/currents_3787.h5';

% path at which to store all the results
options.data_base_path = 'data';

% if 0, do not print anything out
options.print_progress = 1;

% if this is set to 1, it will recompute all the values
% this takes quite a while
options.recompute = 1;

% these are the grid sizes used in the field and gradient experiments
options.grid_sizes = {3,4,5,6};

% needs to be same size as grid_sizes
options.bspline_degrees = {3, 4, 5, 6};

% can be either ieee or thesis
% IEEE is for the T-RO submission
options.plot_mode = 'ieee';

options.ieee_font_size = 8;
options.thesis_font_size = 11;

% colors are stored in 
% data/colors.mat

% the grid size to use for the plots showing the relationship of RBF
% performance on the shape param
options.rbf_eps_grid = 5;
options.rbf_eps_method = 'RBF-MQ-3D';
