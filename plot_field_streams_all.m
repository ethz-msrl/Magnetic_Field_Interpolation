grid_sizes = [3, 4, 6, 8];

clf;
f = figure('PaperUnits', 'inches', 'PaperSize', [6,6], 'color', 'white');

method = 'RBF-G-3D';
current = 5;

for i=1:length(grid_sizes)
    subtightplot(2,2,i, [0.1,0.05]);
    %subplot(2,2,i);
    plot_field_streamslice(method, grid_sizes(i), current, 0, 45);
    xlim([-120, 120]);
    ylim([-120, 120]);
end

export_fig(f, sprintf('figures/streams_%s_%d.pdf', method, current));
