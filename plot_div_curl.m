clear variables;
close all;

MODE = 'ieee';
FONT_SIZE = 8;

load_gradient_data;

idx_r = find(abs(mean_div(:,1)) > 1e-12);

fh = figure('Name', 'Mean Divergence and Curl', 'units', 'inch', ...
    'position', [0,0,7,2.1], 'color', 'w', 'DefaultAxesFontSize', FONT_SIZE);

ax = subplot(1,2,1);

    %colormap(cmap);
b = bar([results.grid_size], 1000*mean_div(idx_r,:)', 'grouped', 'EdgeColor','none');

for i=1:length(idx_r)
    b(i).FaceColor = colors(model_names{idx_r(i)});
end

ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, '$\nabla \cdot \mathbf{b}$ (mT/m)', 'Interpreter', 'latex');
title('Mean Divergence');

%legend(model_names(idx_r));

% if strcmp(MODE, 'ieee')
%     set(ax.Legend.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
% end

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
ax.FontSize = FONT_SIZE;

xticks(ax, cell2mat(grid_sizes));
xlabel(ax, 'Grid Size $N_g$', 'Interpreter', 'latex');
ylabel(ax, '$\| \nabla \times \mathbf{b} \|$ (mT/m)', 'Interpreter', 'latex');
title('Mean Curl Magnitude');
%legend(model_names(idx_r));

% if strcmp(MODE, 'ieee')
%     set(ax.Legend.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
% end

export_fig(fh, 'Figures/interp_div_curl_ieee.pdf')