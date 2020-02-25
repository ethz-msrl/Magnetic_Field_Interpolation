output_files = dir('data/fields/*.mat');

load('data/colors');
load('data/idx');
load('data/grid_sizes');

Nf = length(output_files);
cmap = cbrewer('qual', 'Set1', Nf);

mae = zeros(Nf, length(grid_sizes));
nmae = zeros(Nf, length(grid_sizes));
rmse = zeros(Nf, length(grid_sizes));
nrmse = zeros(Nf, length(grid_sizes));
r2 = zeros(Nf, length(grid_sizes));
mean_div = zeros(Nf, length(grid_sizes));
mean_curl = zeros(Nf, length(grid_sizes));
model_names = {};
for i=1:Nf
    filename = fullfile(output_files(i).folder, output_files(i).name);
    
    temp = strsplit(output_files(i).name, '_');
    model_name = temp{1};
    model_names{i} = model_name;
    noise_std = temp{2} * 1e-6;
    
    results = importdata(filename);
    mae(i,:) = mean([results.mae]);
    nmae(i,:) = mean([results.nmae]);
    rmse(i,:) = mean([results.rmse]);
    nrmse(i,:) = mean([results.nrmse]);
    r2(i,:) = mean([results.r2]);
end
