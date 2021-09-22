function h = cont_heatmap(cont, tick_labels, cmap, hand)
    if exist('hand', 'var')
        h = heatmap(cont, 'Colormap', cmap,'parent', hand);
    else
            h = heatmap(cont, 'Colormap', cmap);
    end
    h.XDisplayLabels = tick_labels;
    h.YDisplayLabels = tick_labels;
    xlabel("neuron k");
    ylabel("neuron j");
    set(gca, 'fontsize', 14);
end