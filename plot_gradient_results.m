load_gradient_data;

MODE = 'ieee';

if strcmp(MODE, 'ieee')
    PAPER_SIZE = [0, 0, 4.6, 2.3];
    FONT_SIZE = 11;
else
    PAPER_SIZE = [0, 0, 3.5, 2.3];
    FONT_SIZE = 8;
end

fh_nmae = figure('Name', 'Mean NMAE', 'units', 'inch', ...
    'position', PAPER_SIZE, 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);

ax = gca;

b = bar([results.grid_size], 100*nmae(idx,:)', 'grouped', 'EdgeColor','none');

% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    b(i).FaceColor = colors(model_names{idx(i)});
end
ax = fh_nmae.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'N-MAE (\%)', 'Interpreter', 'latex');
legend(model_names(idx));

if strcmp(MODE, 'ieee') 
    set(fh_nmae, 'PaperUnits', 'inches');
    set(fh_nmae, 'PaperSize', [3.45, 2.1]);
end
    
export_fig(fh_nmae, sprintf('figures/interp_gradient_nmae_%s.pdf', MODE));

fh_r2 = figure('Name', 'Mean R2', 'units', 'inch', ...
     'position', PAPER_SIZE, 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);

% we sort the by the r2 in the lowest grid resolution
%[~, idx] = sort(r2(:,1), 1, 'descend');
% hold on;
% for i=1:Nf
%     b = bar([results.grid_size], r2(idx(i),:)');
%     b.FaceColor = cmap(idx(i),:);
% end

b = bar([results.grid_size], r2(idx,:)', 'grouped', 'EdgeColor','none');
% don't forget to also sort the colors so they match the other figure
for i=1:length(idx)
    %b(i).FaceColor = cmap(idx(i),:);
    b(i).FaceColor = colors(model_names{idx(i)});
end

ax = fh_r2.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, '$R^2$ ', 'Interpreter', 'latex');
ylim(ax, [min(r2(:))-0.1, 1])
legend(model_names(idx));

%export_fig(fh_r2, sprintf('figures/interp_gradient_r2_%s.pdf', MODE));

%% div
idx_r = find(abs(mean_div(:,1)) > 1e-12);

if strcmp(MODE, 'ieee')
    fh_md = figure('Name', 'Mean Divergence', 'units', 'inch', ...
        'position', [0,0,3.5,2.3], 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);
else
    fh_md = figure('Name', 'Mean Divergence', 'units', 'inch', ...
        'position', PAPER_SIZE, 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);
end
    %colormap(cmap);
b = bar([results.grid_size], 1000*mean_div(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax = fh_md.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');
ylabel(ax, 'Mean Divergence (mT/m)', 'Interpreter', 'latex');
legend(model_names(idx_r));

if strcmp(MODE, 'ieee')
    set(ax.Legend.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
end

export_fig(fh_md, sprintf('figures/interp_divergence_%s.pdf', MODE));

%% curl
idx_r = find(abs(mean_curl(:,1)) > 1e-12);

if strcmp(MODE, 'ieee')
    subplot(1,2,2);
    fh_mc = figure('Name', 'Mean Curl Magnitude', 'units', 'inch', ...
        'position', [0,0,3.5,2.3], 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);
else
    fh_mc = figure('Name', 'Mean Curl Magnitude', 'units', 'inch', ...
        'position', PAPER_SIZE, 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);
end

%colormap(cmap);
b = bar([results.grid_size], 1000*mean_curl(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax = fh_mc.CurrentAxes;
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $n_g$', 'Interpreter', 'latex');

ylabel(ax, 'Mean Curl Magnitude (mT/m)', 'Interpreter', 'latex');
legend(model_names(idx_r));


if strcmp(MODE, 'ieee')
    set(ax.Legend.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
end

%export_fig(fh_mc, sprintf('figures/interp_curl_%s.pdf', MODE));