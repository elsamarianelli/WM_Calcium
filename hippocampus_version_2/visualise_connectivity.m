function connectivity_plot = visualise_connectivity(C, CA3_populations, CA1_populations)
%% visualise connectivity between CA3 and CA1

    connectivity_plot = figure;
    connectivity_matrix = C;
    length_show = 20;
    connectivity_matrix = connectivity_matrix(1:length_show, 1:length_show);
    CA3_populations{1}(CA3_populations{1}>length_show) = [];
    CA3_populations{2}(CA3_populations{2}>length_show) = [];
    CA3_populations{3}(CA3_populations{3}>length_show) = [];
    CA1_populations{1}(CA1_populations{1}>length_show) = [];
    CA1_populations{2}(CA1_populations{2}>length_show) = [];
    CA1_populations{3}(CA1_populations{3}>length_show) = [];


    % Get the size of the matrix
    [num_inputs, num_outputs] = size(connectivity_matrix);
    
    % Define positions of input and output cells
    input_positions = 1:num_inputs;
    output_positions = (1:num_outputs);
    
    % Plot connections
    for i = 1:num_inputs
        for j = 1:num_outputs
            if connectivity_matrix(i, j) == 1
                plot([input_positions(i), output_positions(j)], [1, 2], 'k');
                hold on;
            end
        end
    end
    
    hold on;

    % Plot input cells
    for i = 1:length(input_positions)
        position = input_positions(i);
        if ismember(position, CA3_populations{1})
            scatter(position, 1, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'blue'); hold on;
        elseif ismember(position, CA3_populations{2})
            scatter(position, 1, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'red'); hold on;
        elseif ismember(position, CA3_populations{3})
            scatter(position, 1, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'magenta'); hold on;     
        else
            scatter(position, 1, 80, 'filled', 'Marker', 'o','MarkerFaceColor', [.5 .5 .5]); hold on;     
        end
    end
    
    % Plot output cells
    for i = 1:length(output_positions)
        position = output_positions(i);
        if ismember(position, CA1_populations{1})
            scatter(position, 2, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'blue'); hold on;
        elseif ismember(position, CA1_populations{2})
            scatter(position, 2, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'red'); hold on;
        elseif ismember(position, CA1_populations{3})
            scatter(position, 2, 80, 'filled', 'Marker', 'o','MarkerFaceColor', 'magenta'); hold on;     
        else
            scatter(position, 2, 80, 'filled', 'Marker', 'o','MarkerFaceColor', [.5 .5 .5]); hold on;     
        end
    end

    % Customize plot
    xlim([0, num_outputs + 1]);
    ylim([0, 3]);
    xlabel('Cells');
    ylabel('Layers');
    title('Connectivity Visualization');
    set(gca, 'YTick', [1, 2], 'YTickLabel', {'Input', 'Output'});
    
end