% This makes a figure with side by side bar plots comparing the mean value 
% of the absolute value of the divergence (so that negative values don't
% compensate for positive ones when taking the mean) for all methods. It
% also shows the mean value of the magnitude of the curl vector for all
% methods.
%
% Note that some methods are naturally curl or divergence-free. Those
% methods are ommitted from the plot.
% 
% This plot is intendended for a double column
% IEEE figure. 
%
% The divergence and curl are averaged over all positions and currents.
%
% Note: because I couldn't find a way of putting the same legend as the
% field and gradient plots on this figure, I disabled the legend. I just
% copy it over manually onto the figure in Illustrator later.
%
%   Copyright 2020, Samuel Charreyron

clear variables;
close all;

load_gradient_data;

% filtering out values where the divergence is essentially 0 so those
% methods don't appear in the plot
idx_r = find(abs(mean_div(:,1)) > 1e-12);

fh = figure('Name', 'Mean Divergence and Curl', 'units', 'inch', ...
    'position', [0,0,7,2.1], 'color', 'w', 'DefaultAxesFontSize', options.ieee_font_size);

ax = subplot(1,2,1);

b = bar([results.grid_size], 1000*mean_div(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = options.ieee_font_size;

xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, '$ | \nabla \cdot \mathbf{b} |$ (mT/m)', 'Interpreter', 'latex');
title('Mean Absolute Divergence');

%% curl
% here we ignore all the ones that are zero
idx_r = find(abs(mean_curl(:,1)) > 1e-12);

ax = subplot(1,2,2);

%colormap(cmap);
b = bar([results.grid_size], 1000*mean_curl(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = options.ieee_font_size;

xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, '$\| \nabla \times \mathbf{b} \|$ (mT/m)', 'Interpreter', 'latex');
title('Mean Curl Magnitude');

export_fig(fh, 'Figures/interp_div_curl_ieee.pdf')