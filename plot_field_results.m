load_field_data;

%% nmae
fh_mae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 11);
colormap(cmap);
[~, idx] = sort(nmae(:,1), 1);

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = colors(model_names{idx(i)});
end
ax = fh_mae.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';

xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
legend(model_names(idx));

%set(fh_mae, 'PaperUnits', 'inches');
%set(fh_mae, 'PaperSize', [3.45/2, 2.1]);

export_fig(fh_mae, 'figures/interp_field_nmae.pdf');

%% r2

fh_r2 = figure('Name', 'Mean R2', 'units', 'inch', ...
     'position', [0, 0, 4.6, 3], 'color', 'w', 'DefaultAxesFontSize', 8);
colormap(cmap);
% we sort the by the r2 in the lowest grid resolution
%[~, idx] = sort(r2(:,1), 1, 'descend');
b = bar([results.grid_size], r2(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = colors(model_names{idx(i)});
end
ax = fh_r2.CurrentAxes;
xticks(ax, cell2mat(options.grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, '$R^2$ ', 'Interpreter', 'latex');
ylim(ax, [min(r2(:))-0.1, 1])
legend(model_names(idx));

export_fig(fh_r2, 'figures/interp_field_r2.pdf');
