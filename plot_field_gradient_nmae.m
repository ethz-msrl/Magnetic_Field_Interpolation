clear variables;
close all;

MODE = 'ieee';

if strcmp(MODE, 'ieee')
    FONT_SIZE = 8;
else
    FONT_SIZE = 11;
end

load_field_data;

fh_nmae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', [0, 0, 7, 2.2], 'color', 'w', 'DefaultAxesFontSize', 8);
ax = subplot(1,2,1);
colormap(cmap);
[~, idx] = sort(nmae(:,1), 1);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = cmap(i,:);
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
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
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
title('Magnetic Field Gradient');

% if strcmp(MODE, 'ieee') 
%     set(fh_nmae, 'PaperUnits', 'inches');
%     set(fh_nmae, 'PaperSize', [3.45, 2.1]);
% end
%
set(fh_nmae, 'PaperUnits', 'inches');
set(fh_nmae, 'PaperSize', [7, 2.1]);

export_fig(fh_nmae, sprintf('figures/interp_field_gradient_nmae_%s.pdf', MODE));