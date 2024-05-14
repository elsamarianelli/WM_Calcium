%% Code to train perceptron on odour discrimination task
%% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 2000;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;           % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 3000;
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

%% check dynamics over single trial 

% set extra paramaters for single trial
p.pattern_order         = 'AB';
input.simulation        = [p.start_time p.start_time+p.length_first];
input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
M                       = get_memory_hipp(p);

% take first and second odour to be presented 
stim                    = cell(2,1); 
first                   = double(upper(p.pattern_order(1))) - 64; 
stim{1}                 = ca3_ensembles{first}; clear first
second                  = double(upper(p.pattern_order(2))) - 64; 
stim{2}                 = ca3_ensembles{second}; clear second

% simulate dynamics
M                       = simulate_dynamics_hipp(p, C, J, input, M, stim);
output_plot             = get_output_plot(M,p.pattern_order, p, stim, C, ca3_ensembles, ca1_ensembles);


% get times to use to train perceptorn 
time1 = input.simulation(1);
time2 = input.simulation(2);
time3 = input.reactivation(1);
time4 = time3+100; %input.reactivation(2);

% % %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
% [spikeCounts,sequence]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
% [spikeCounts_test, ~]	= get_train_data_db(C, J, 10, p, ca3_ensembles, time3, time4);
% 
%% Save data and settings info to current working directory
folderName = 'fixIn_true_100_after_odour_2000_delay';

% Create the folder if it doesn't already exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Save the matrix in a .mat file
save(fullfile(folderName, 'spikeCounts.mat'), 'spikeCounts');
save(fullfile(folderName, 'spikeCounts_test.mat'), 'spikeCounts_test');
% Save the structure in a .mat file
save(fullfile(folderName, 'myStruct.mat'), 'p');

%% creating delay data using parpool

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;           % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
    

% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of delays
delays = 0:250:3000;

% Start a parallel pool if not already started
pool = gcp('nocreate'); 

if isempty(pool)
    parpool; 
end

% Parallel loop
parfor index = 2:3%length(delays)
    delay_time = delays(index);
    save_data_varying_delay_times(C, J, p, ca3_ensembles, delay_time, n_trials, main_folder);
end

%% plotting 

% Number of colors needed
ncols = 5;
indices = linspace(1, 64, ncols); % MATLAB colormaps by default contain 64 colors
cmap = winter(64);  % Get the full colormap`
selectedColors = cmap(ceil(indices), :);  % Use ceil to ensure indices are integers

% Define scenarios and their corresponding delays
data_loaded = {
    'fixIn_true_100_after_odour_300_delay', [3000, selectedColors(1, 1:3)]; %300 delay
    '300_trials_fixIn_true_100_after_odour_shorter_delay', [3000,selectedColors(2, 1:3)]; % short delay (700)
    '300_trials_fixIn_true_100_after_odour', [3000,selectedColors(3, 1:3)]; % 1500 delay
     'fixIn_true_100_after_odour_2000_delay', [3000, selectedColors(4, 1:3)]; % 2000 delay
    '300_trials_fixIn_true_100_after_odour_3000_delay', [3000,selectedColors(5, 1:3)]; % long delay (3000)
};

figure;

% Loop through each scenario
for idx = 1:size(data_loaded, 1)

    folderName = data_loaded{idx, 1};
    smoothVal = data_loaded{idx, 2}(1);
    color = data_loaded{idx, 2}(2:4);

    % Load data
    [data, data_test] = loadData(folderName);

    % Train and test single layer perceptron
    [performance_accuracy_single, error, w] = run_perceptron_db(data);
    [performance_test_single] = test_perceptron_output(data_test, w);
    plotError(error, idx * 2 - 1, color, performance_test_single, ncols);

    % Train and test multilayer perceptron
    [~, error, w1, w2] = run_multilayer_perceptron(data);
    [performance_test_multi] = test_multilayer_perceptron_output(data_test, w1, w2);
    plotError(error, idx * 2, color, performance_test_multi,ncols);

end

%% multiple runs and plot mean and SD of performance

numTests = 20;

% Initialize containers for performance metrics
singlePerf = zeros(size(data_loaded, 1), numTests);
multiPerf = zeros(size(data_loaded, 1), numTests);

% Loop through each scenario
for idx = 1:size(data_loaded, 1)

    folderName = data_loaded{idx, 1};

    % Load data
    [data, data_test] = loadData(folderName);

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
meanSingle = mean(singlePerf, 2);
stdSingle = std(singlePerf, 0, 2);
meanMulti = mean(multiPerf, 2);
stdMulti = std(multiPerf, 0, 2);
means = [meanSingle; meanMulti];
stds = [stdSingle; stdMulti];
colors = [selectedColors;selectedColors];    
% labels for each delay time
labels = {'300'
          '700'
          '1500'
          '2000'
          '3000'};
labels = [labels;labels];
% Plot results
figure;
hold on;
for i = 1:length(means)
    errorbar(i, means(i), stds(i), 'o', 'Color', colors(i, 1:3), 'MarkerSize', 10, 'MarkerFaceColor', colors(i, 1:3));
end
set(gca, 'FontSize', 10);
xticks(1:length(labels));
xticklabels(labels);
ylabel('Performance', 'FontSize', 10);
title('Performance Comparison Across Scenarios');
hold off;

