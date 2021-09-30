%% load and format
load('data\ori1_data_session_3.mat')
[~, st_large_idx] = sort(sum(logical(data.spiketrain(:,:,data.Cluster_Included)), [2,3]));
condition_idx = st_large_idx(end-10); % get the condition with the 10th most spikes
spiketrain_subset = squeeze(logical(data.spiketrain(condition_idx,:,data.Cluster_Included)));
depths = data.Cluster_celldepth(data.Cluster_Included);
[~, i] = sort(depths);

sorted_depths = depths(i);
sorted_sts = spiketrain_subset(:,i);

%% plot
figure('position', [   360   198   314   500]);
for neu = 1:size(sorted_sts,2)
    for time = -100:999
        if sorted_sts(time+101,neu)
            plot([time,time],[neu,neu+1], 'linewidth', 1,'color', 'k');
            hold on;
        end
    end
end
xlim([400,1000]);
ylim([1,size(sorted_sts,2)]);
xline(0, '--k', 'linewidth', 1);
set_axis_defaults();

%xlabel("time rel. to stim. onset (ms)", 'fontsize', 12);
%ylabel("<- closest - neuron to probe tip - farthest ->", 'fontsize', 12);
set(gca, 'fontsize', 12);
set(gca, 'linewidth', 1);
save_close_figures("figures\matlab\spike_raster");