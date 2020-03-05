% This makes a figure with side by side bar plots comparing the interpolation
% performance of all interpolation methods of the magnetic field and
% magnetic field gradient. This plot is intendended for a double column
% IEEE figure. 
%
% The plots show the mean absolute error normalized by its range so that
% error is unitless. For readability, we show it as a percentage.
% The error metric is averaged over all positions, currents, and components. 
% For the gradients we also average over all components of the gradient.
%
%   Copyright 2020, Samuel Charreyron

clear variables;
close all;
load_settings;

%% Fields

load_field_data;

fh_nmae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', [0, 0, 7, 2.2], 'color', 'w', 'DefaultAxesFontSize', 8);
ax = subplot(1,2,1);
colormap(cmap);

% the idx variable controls the order in which the bars appear.
% We actually set this and saved it in the data folder so that the bars
% roughly appear in increasing order from left to right at higher grid
% sizes
load(sprintf('%s/idx.mat', options.data_base_path));
%[~, idx] = sort(nmae(:,3), 1);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = colors(model_names{idx(i)});
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = opt.ieee_font_size;

xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
title('Magnetic Field');
legend(model_names(idx));


%% Gradients
load_gradient_data;

ax = subplot(1,2,2);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');

% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = colors(model_names{idx(i)});
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = options.ieee_font_size;

xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
title('Magnetic Field Gradient');

% if strcmp(options.plot_mode, 'ieee') 
%     set(fh_nmae, 'PaperUnits', 'inches');
%     set(fh_nmae, 'PaperSize', [3.45, 2.1]);
% end
%
set(fh_nmae, 'PaperUnits', 'inches');
set(fh_nmae, 'PaperSize', [7, 2.1]);

export_fig(fh_nmae, sprintf('figures/interp_field_gradient_nmae_%s.pdf', options.plot_mode));