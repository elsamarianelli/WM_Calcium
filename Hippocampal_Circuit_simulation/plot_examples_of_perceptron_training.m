function[fig_handle] = plot_examples_of_perceptron_training(variable_range, main_folder, cmap, variable_type)

% Number of colors needed
ncols = length(variable_range);
indices = linspace(1, 64, ncols); 
selectedColors = cmap(ceil(indices), :);  

figure;

% Loop through each scenario
for idx = 1%:4:size(variable_range, 2)

    variable = variable_range(idx);
    folderName = num2str(variable);
    
    % Full path for the new folder
    fullFolderPath = fullfile(main_folder, folderName);   
    color =selectedColors(idx, 1:3);

    % Load data
    [data, data_test] = loadData(fullFolderPath);

    % Train and test single layer perceptron
    [~, error, w] = run_perceptron_db(data);
    [performance_test_single] = test_perceptron_output(data_test, w);
    plotError(error, idx * 2 - 1, color, performance_test_single, ncols);

    % Train and test multilayer perceptron
    [~, error, w1, w2] = run_multilayer_perceptron(data);
    [performance_test_multi] = test_multilayer_perceptron_output(data_test, w1, w2);
    plotError(error, idx * 2, color, performance_test_multi,ncols);
    hold on
end

sgtitle(['perceptron performance over training for different ' variable_type])
fig_handle = gcf;

end