function mean_Hz = get_mean_rate(M, time_from, time_to, cell_index)
% function gives mean firing rate (Hz of specified population of cells in
% specified time frame
    spikelog = M.spikelog;
    spikelog_of_interest = spikelog(cell_index, time_from:time_to);
    mean_Hz = sum(spikelog_of_interest(:))./ (length(cell_index).*(time_to-time_from).*0.001);
end
