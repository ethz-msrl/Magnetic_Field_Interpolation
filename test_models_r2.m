clear all;

NODES_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_4_h5/';
EVAL_DATASET = '/Volumes/msrl/users/samuelch/datasets/cmag_calibration/mpem_synthetic_16_h5/';

% load currents file
currents = h5read('/Volumes/msrl/users/samuelch/datasets/cmag_calibration/currents_3787.h5', '/currents')';

% load positions
nodes_pos_fn = fullfile(NODES_DATASET,'/positions.h5');
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
EPS = 10.;

r2_scores = zeros(NUM_CURRENTS, 3);
% 
% for i=1:NUM_CURRENTS
%     fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields = permute(fields, [4, 3, 2, 1]);
%     
%     fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields_ev = permute(fields_ev, [4, 3, 2, 1]);
%     
%     model = SimpleRBFInterpolator(nodes, fields, EPS);
%     ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
%     
%     r2_scores(i,:) = ev.get_r2();
% end
% 
% disp('Simple RBF mean r2 scores:')
% disp(mean(r2_scores, 1));
% 
% EPS = 21.;
% 
% r2_scores = zeros(NUM_CURRENTS, 3);
% 
% for i=1:NUM_CURRENTS
%     fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields = permute(fields, [4, 3, 2, 1]);
%     
%     fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields_ev = permute(fields_ev, [4, 3, 2, 1]);
%     
%     model = DivFreeRBFInterpolator(nodes, fields, EPS);
%     ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
%     
%     r2_scores(i,:) = ev.get_r2();
% end
% 
% disp('Div-Free RBF mean r2 scores:')
% disp(mean(r2_scores, 1));

% for i=1:NUM_CURRENTS
%     fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields = permute(fields, [4, 3, 2, 1]);
%     
%     fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields_ev = permute(fields_ev, [4, 3, 2, 1]);
%     
%     model = BSpline3DInterpolator(nodes, fields, 4);
%     ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
%     
%     r2_scores(i,:) = ev.get_r2();
% end
% 
% disp('BSpline3D mean r2 scores');
% disp(mean(r2_scores, 1));

% for i=1:NUM_CURRENTS
%     fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields = permute(fields, [4, 3, 2, 1]);
%     
%     fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields_ev = permute(fields_ev, [4, 3, 2, 1]);
%     
%     model = BSplineLaplacianInterpolator(nodes, fields, 4);
%     ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
%     
%     r2_scores(i,:) = ev.get_r2();
% end
% 
% disp('BSplineLaplacian mean r2 scores');
% disp(mean(r2_scores, 1));

%[M] = get_tricubic_3d_matrix();
load('tricubic_simple_M.mat', 'M');

for i=1:NUM_CURRENTS
    disp(i);
    fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
    fields = permute(fields, [4, 3, 2, 1]);
    
    fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
    fields_ev = permute(fields_ev, [4, 3, 2, 1]);
    
    model = SimpleTricubicInterpolator(nodes, fields, M);
    ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
    
    r2 = ev.get_r2();
    disp(r2);
    r2_scores(i,:) = r2;
end

disp('Simple tricubic scores:')
disp(mean(r2_scores, 1));

% M = get_tricubic_scalar_field_matrix();
% 
% for i=1:NUM_CURRENTS
%     disp(i);
%     fields = h5read(fullfile(NODES_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields = permute(fields, [4, 3, 2, 1]);
%     
%     fields_ev = h5read(fullfile(EVAL_DATASET, 'v', sprintf('%04d.h5', i)), '/fields');
%     fields_ev = permute(fields_ev, [4, 3, 2, 1]);
%     
%     model = TricubicScalarFieldInterpolator(nodes, fields, M);
%     ev = FieldInterpolatorEvaluator(model, positions_ev, fields_ev);
%     
%     r2 = ev.get_r2();
%     disp(r2);
%     r2_scores(i,:) = r2;
% end
% 
% disp('Tricubic scalar field scores:')
% disp(mean(r2_scores, 1));