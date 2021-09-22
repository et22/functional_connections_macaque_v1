function plot_cluster_axis(smoothed_ccgs, idx, num_clusters, cmap, sig_stars)
if ~exist('sig_stars', 'var')
    sig_stars = zeros(1,num_clusters);
end
pos = get(gca, 'position');
fig_pos = get(gcf, 'position');
w_to_h = fig_pos(3)/fig_pos(4);
width = pos(3)/(num_clusters*1.2);
height = width*1.2*w_to_h;
set(gca,'position', [pos(1),pos(2)+height, pos(3), max([0,pos(4)-height])]);
pos = get(gca, 'position');
for i = 1:num_clusters
    ccg = mean(smoothed_ccgs(idx==i,:),1);
    xstart = pos(1) + (2*i-1)*pos(3)/(num_clusters*2)-pos(3)/(num_clusters*2.4);
    ystart = pos(2)-height*1.1;

    axes('Position',[xstart, ystart, width, height]);
    plot(ccg, 'color', cmap(i,:), 'linewidth', 2);
    if sig_stars(i)
        ylims = ylim;
        xlims = xlim;
        text(gca, xlims(2)*.85, ylims(2)*.85, '*', 'fontsize', 12);
    end
    set(gca, 'linewidth',1, 'box', 'off', 'tickdir', 'out', 'xticklabel', [], 'yticklabel', [], 'xcolor', 'none', 'ycolor', 'none', "Color","none")
    xlabel('');
    ylabel('');
end

end