function[fig_handle ,means_perf, stds_perf] = plot_performance_across_variable_change(variable_range, main_folder, variable_type)
%% multiple runs and plot mean and SD of performance

numTests = 3;

% Initialize containers for performance metrics
% singlePerf = zeros(size(variable_range, 2), numTests);
multiPerf = zeros(size(variable_range, 2), numTests);
pleateau_iter_log =  zeros(size(variable_range, 2), numTests);

% Loop through each scenario
for idx = 1:size(variable_range, 2)
    
    variable = variable_range(idx);

    if contains(main_folder, 'ux_varied')
        folderName = variable;
        folderName = folderName{1};
    else
        folderName = num2str(variable);
    end


    % Full path for the new folder
    fullFolderPath = fullfile(main_folder, folderName);   


    % Load data
    [data, data_test] = loadData(fullFolderPath, main_folder);
    disp(variable)
    % Repeat tests
    for j = 1:numTests
        % % Train and test single layer perceptron
        % [~, ~, w] = run_perceptron_db(data);
        % singlePerf(idx, j) = test_perceptron_output(data_test, w);

        % Train and test multilayer perceptron
        [~, ~, w1, w2, plateau_iter] = run_multilayer_perceptron(data);
        pleateau_iter_log(idx, j)    = plateau_iter;
        multiPerf(idx, j) = test_multilayer_perceptron_output(data_test, w1, w2);
        disp(j)
    end
    
end

%% plot figure
figure;

% Colors
perf_color = [0 0 0]; % Dark for performance
learn_color = [.6 .6 .6]; % Lighter for learning time

% [1] plot STD and means for performance
% Calculate means and standard deviations for the performance
means_perf = mean(multiPerf, 2);
stds_perf = std(multiPerf, 0, 2);


labels = variable_range;

% Adding a line connecting the means
yyaxis left; % Left y-axis for performance
for i = 1:length(means_perf)-1
    % taking average color between two points
    plot([i, i+1], means_perf(i:i+1), '-', 'Color', perf_color, 'LineWidth', 2);
    hold on
end
ylim([0 1])

% plotting error bars and means
for i = 1:length(means_perf)
    errorbar(i, means_perf(i), stds_perf(i), 'o', 'Color', perf_color, 'MarkerSize', 10, 'MarkerFaceColor', perf_color);
    hold on
end

% Set left y-axis color
ax = gca;
ax.YColor = perf_color;

% [2] plot mean and STD for learning time
% Calculate means and standard deviations for the learning times
means_learn = mean(pleateau_iter_log, 2);
stds_learn = std(pleateau_iter_log, 0, 2);

% Adding a line connecting the means
yyaxis right; % Right y-axis for learning time
for i = 1:length(means_learn)-1
    % taking average color between two points
    plot([i, i+1], means_learn(i:i+1), '-', 'Color', learn_color, 'LineWidth', 2);
    hold on
end
ylim([0 3000]);

% plotting error bars and means
for i = 1:length(means_learn)
    errorbar(i, means_learn(i), stds_learn(i), 'o', 'Color', learn_color, 'MarkerSize', 10, 'MarkerFaceColor', learn_color);
    hold on
end

% Set right y-axis color and limits
ax.YColor = learn_color;
% ylim([0, 1500]); % Set y-axis limits for learning time

% [3] change plot appearance
set(gca, 'FontSize', 12);
xticks(1:length(labels));
xticklabels(labels);

% Labels and title
yyaxis left;
ylabel('Performance', 'FontSize', 12);

yyaxis right;
ylabel('Learning Time', 'FontSize', 12);

xlabel(variable_type, 'FontSize', 12);
title(variable_type);

xlim([0 length(variable_range)+1]);
% get figure handle 
fig_handle = gcf;

hold off

end