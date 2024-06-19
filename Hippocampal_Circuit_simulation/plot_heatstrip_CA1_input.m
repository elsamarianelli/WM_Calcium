function[sortIdx] = plot_heatstrip_CA1_input(input_1, input_2, input_label, pair_spikes)

% get colour map
colormap(flipud(gray)); % Use 'gray' colormap and flip it for dark to light mapping
cmap = colormap;

% sort input 
[input_2_sorted, sortIdx] = sort(input_2, 'descend');
input_1_sorted = input_1(sortIdx);
pair_spikes = pair_spikes(sortIdx);
% heat_subplot = figure; % Make the figure invisible

% Function to plot color strip with second odour as sorting priority
plot_color_strip = @(values, y_position) arrayfun(@(i) ...
    rectangle('Position', [i-1, y_position, 1, 1], 'FaceColor', cmap(round(values(i) * (size(cmap, 1) - 1)) + 1, :), 'EdgeColor', 'none'), ...
    1:length(values));

% Plot the color strips for A, B, and C
plot_color_strip(input_1_sorted/20, 1);
plot_color_strip(input_2_sorted/20, 2);
% Set x-axis limits and customize y-axis
set(gca, 'ytick', [0.5, 1.5, 2.5], 'yticklabel', {input_label(1), input_label(2)});
hold on
plot(1:length(pair_spikes), pair_spikes / max(pair_spikes) * 2 + 3, 'k'); % Adjust scaling as needed
ylim([1 5])

