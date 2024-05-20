function[fig_handle] = plot_performance_across_variable_change(variable_range, main_folder, cmap, variable_type)


% Number of colors needed
ncols = length(variable_range);
% indices = linspace(1, 64, ncols);
% colors = cmap(ceil(indices), :); 
colors = zeros(ncols, 3);
colors_learning = 0.5*zeros(ncols, 3);
%% multiple runs and plot mean and SD of performance

numTests = 1;

% Initialize containers for performance metrics
% singlePerf = zeros(size(variable_range, 2), numTests);
multiPerf = zeros(size(variable_range, 2), numTests);
pleateau_iter_log =  zeros(size(variable_range, 2), numTests);

% Loop through each scenario
for idx = 1:size(variable_range, 2)
    
    variable = variable_range(idx);
    folderName = num2str(variable);
    % folderName = variable{1};

    % Full path for the new folder
    fullFolderPath = fullfile(main_folder, folderName);   

    % Load data
    [data, data_test] = loadData(fullFolderPath);
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

% Plot results
figure;

% [1] plot STD and means for performance

% Calculate means and standard deviations for the performance
means = mean(multiPerf, 2);
stds = std(multiPerf, 0, 2);

% labels for each delay time
labels = variable_range;

% Adding a line connecting the means
for i = 1:length(means)-1
    % taking average color between two  points
    interpColor = (colors(i, 1:3) + colors(i+1, 1:3)) / 2;
    plot([i, i+1], means(i:i+1), '-', 'Color', interpColor, 'LineWidth', 2);
    hold on
end

% plotting error bars and means
for i = 1:length(means)
    errorbar(i, means(i), stds(i), 'o', 'Color', colors(i, 1:3), 'MarkerSize', 10, 'MarkerFaceColor', colors(i, 1:3));
    hold on
end

% [2] plot mean and STD for learning time

% Calculate means and standard deviations for the learning times
means = mean(pleateau_iter_log, 2);
stds = std(pleateau_iter_log, 0, 2);

% Adding a line connecting the means
for i = 1:length(means)-1
    % taking average color between two  points
    interpColor = (colors_learning(i, 1:3) + colors_learning(i+1, 1:3)) / 2;
    plot([i, i+1], means(i:i+1), '-', 'Color', interpColor, 'LineWidth', 2);
    hold on
end

% plotting error bars and means
for i = 1:length(means)
    errorbar(i, means(i), stds(i), 'o', 'Color', colors_learning(i, 1:3), 'MarkerSize', 10, 'MarkerFaceColor', colors_learning(i, 1:3));
    hold on
end

% [3] change plot appearance
set(gca, 'FontSize', 12);
xticks(1:length(labels));
xticklabels(labels);
ylabel('Performance', 'FontSize', 12);
xlabel(variable_type, 'FontSize', 12)
title(['Multilayer Perceptron performance at different ', variable_type]);

% get figure handle 
fig_handle = gcf;

hold off

end