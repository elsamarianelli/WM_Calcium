function[fig_handle] = plot_performance_across_variable_change(variable_range, main_folder, cmap, variable_type)

% Number of colors needed
ncols = length(variable_range);
indices = linspace(1, 64, ncols);
colors = cmap(ceil(indices), :); 

%% multiple runs and plot mean and SD of performance

numTests = 15;

% Initialize containers for performance metrics
singlePerf = zeros(size(variable_range, 2), numTests);
multiPerf = zeros(size(variable_range, 2), numTests);

% Loop through each scenario
for idx = 1:size(variable_range, 2)

    variable = variable_range(idx);
    folderName = num2str(variable);
    
    % Full path for the new folder
    fullFolderPath = fullfile(main_folder, folderName);   

    % Load data
    [data, data_test] = loadData(fullFolderPath);

    % Repeat tests
    for j = 1:numTests
        % Train and test single layer perceptron
        [~, ~, w] = run_perceptron_db(data);
        singlePerf(idx, j) = test_perceptron_output(data_test, w);

        % Train and test multilayer perceptron
        [~, ~, w1, w2] = run_multilayer_perceptron(data);
        multiPerf(idx, j) = test_multilayer_perceptron_output(data_test, w1, w2);
    end
end

% Calculate means and standard deviations
means = mean(multiPerf, 2);
stds = std(multiPerf, 0, 2);

% labels for each delay time
labels = variable_range;

% Plot results
figure;

hold on;
for i = 1:length(means)
    errorbar(i, means(i), stds(i), 'o', 'Color', colors(i, 1:3), 'MarkerSize', 10, 'MarkerFaceColor', colors(i, 1:3));
end
set(gca, 'FontSize', 10);
xticks(1:length(labels));
xticklabels(labels);
xlim([0 length(variable_range)+1]);
ylim([0 1])
ylabel('Performance', 'FontSize', 10);
xlabel('Delay time', 'FontSize', 10)
title(['Multilayer Perceptron performance at different ', variable_type]);

% get figure handle 
fig_handle = gcf;

hold off;

end