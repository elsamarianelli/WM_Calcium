%% looking into MLP behaviour

% check dynamics over a single trial
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;            % Length of time for which the first odour is presented (ms)
p.delay_time            = 500;            % Delay between odour presentations (ms)
p.length_second         = 250;            % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;          % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 2000;
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

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
output_plot             = get_output_plot(M,p,ca3_ensembles, C);

% looking at performance for a single setting
time3 = input.reactivation(1);
time4 = time3+100; 

%% simulate data
n_trials = 50;
[spikeCounts,sequenceID]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
[spikeCounts_test, sequenceID_test]	= get_train_data_db(C, J, 10, p, ca3_ensembles, time3, time4);
parent_dir = fileparts(pwd);
    
% Full path for the new folder
fullFolderPath = fullfile('Simulated_data', 'false_CM_data');

% Check if the directory exists; if not, create it
if ~exist(fullFolderPath, 'dir')
    mkdir(fullFolderPath);
end

% Save in .mat file
save(fullfile(fullFolderPath, 'spikeCounts.mat'), 'spikeCounts');
save(fullfile(fullFolderPath, 'spikeCounts_test.mat'), 'spikeCounts_test');
save(fullfile(fullFolderPath, 'sequenceID.mat'), 'sequenceID');
save(fullfile(fullFolderPath, 'sequenceID_test.mat'), 'sequenceID_test');

%%  how does CA1 activity during second odour correlate with CA1 input
%%  from selective CA3 populaions for each odour combination?

spikes = spikeCounts_test;
IDs = sequenceID_test;

uniqueIDs = unique(sequenceID);
spike_mean_log = zeros(length(uniqueIDs), 200);

% Separate spiking into  blocks based on sequence IDs
for i = 1:length(uniqueIDs)
    currentID = uniqueIDs(i);
    spikes_for_trial = spikes(find(strcmp(sequenceID_test,currentID)), 1:200);
    spike_mean_log(i, 1:200) = mean(spikes_for_trial);
end

% get CA1 populaitons which receive input from each CA3 group
A = C(ca3_ensembles{1}, :);
CA1_A = sum(A,1);
B = C(ca3_ensembles{2}, :);
CA1_B = sum(B,1);
C_ = C(ca3_ensembles{3}, :);
CA1_C = sum(C_,1);

% stringList     = {'CB'; 'BA'; 'AC'; 'AB'; 'CA'; 'BC'};
stringList = {'AB'; 'AC'; 'BA'; 'BC'; 'CA'; 'CB'};
% Create a figure
figure;
% Loop through each string in the list
hold on
for i = 1:length(stringList)
    str = stringList{i};
    
    % Generate the inputs dynamically
    first_input = eval(['CA1_', str(1)]);
    second_input = eval(['CA1_', str(2)]);
    input_label = str;
    pair_spikes = spike_mean_log(i, :);
    % plot
    subplot(3, 2, i);
    [sortIdx] = plot_heatstrip_CA1_input(first_input, second_input, input_label, pair_spikes);
    hold on;
end



%% 
% Example data
data = sorted_spikes;
% Define vertical offset
offset = 2;  % Adjust this value to control the spacing between the lines

% Create a figure
figure;

% Hold on to plot multiple lines on the same axes
hold on;

% Plot each row of the matrix with vertical offset
for i = 1:size(data, 1)
    plot(data(i, :) + (i-1) * offset, 'k');  % 'k' specifies black color for the lines
end

% Customize y-axis to show row labels
yticks((0:size(data, 1)-1) * offset);
yticklabels(arrayfun(@(x) sprintf('Row %d', x), 1:size(data, 1), 'UniformOutput', false));

% Add labels and title
xlabel('Columns');
ylabel('Rows');
title('Matrix Rows Plotted as Separate Lines with Offset');

% Adjust y-axis limits to fit all lines
ylim([-offset, size(data, 1) * offset]);

% Hold off to stop adding to the current plot
hold off;

%% looking at CA1 populations 
% no overlap populations

A = C(ca3_ensembles{1}, :);
CA1_A = find(sum(A,1)>0);
B = C(ca3_ensembles{2}, :);
CA1_B = find(sum(B,1)>0);
C_= C(ca3_ensembles{3}, :);
CA1_C = find(sum(C_,1)>0);

%overlap populations
A_B = find(sum(C(intersect(ca3_ensembles{1}, ca3_ensembles{2}), :))>0);
B_C = find(sum(C(intersect(ca3_ensembles{2}, ca3_ensembles{3}), :))>0);
A_C = find(sum(C(intersect(ca3_ensembles{1}, ca3_ensembles{3}), :))>0);

CA1_AB_A = intersect(A_B, CA1_A);
CA1_AC_A = intersect(A_C, CA1_A);
CA1_BC_B = intersect(B_C, CA1_B);
CA1_AB_B = intersect(A_B, CA1_B);
CA1_BC_C = intersect(B_C, CA1_C);
CA1_AC_C = intersect(A_C, CA1_C);

CA1_A_only = setdiff(CA1_A, CA1_AB_A, CA1_AC_A);
CA1_B_only = setdiff(CA1_B, CA1_BC_B, CA1_AB_B);
CA1_C_only = setdiff(CA1_C, CA1_BC_C, CA1_AC_C);
%%
[~, error, w1, w2, plateau_iter] = run_multilayer_perceptron(spikeCounts, 6);
[performance_test_multi, hidden_output_log, hidden_input_log] = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);

disp(performance_test_multi)

% Get unique sequence IDs
uniqueIDs = unique(sequenceID);

% Initialize a cell array to store blocks of activation numbers
hidden_output_mean_log = zeros(length(uniqueIDs), 6);
output_input_log = zeros(length(uniqueIDs), 6);
hidden_input_mean_log = zeros(length(uniqueIDs), 200);
input_ouput_log = zeros(length(uniqueIDs), 200);
activationNumbers = hidden_output_log;

% Separate activation numbers into blocks based on sequence IDs
for i = 1:length(uniqueIDs)
    currentID = uniqueIDs(i);

    % input layer inputs 
    hidden_input = hidden_input_log(find(strcmp(sequenceID_test,currentID)), 1:200);
    hidden_input_mean_log(i, 1:200) = mean(hidden_input);


    % hidden layer inputs
    activation_values = activationNumbers(find(strcmp(sequenceID_test,currentID)), 1:6);
    hidden_output_mean_log(i, 1:6) = mean(activation_values);

    % output layer inputs 
    output_input = (mean(activation_values)).* w2'; 
    output_input_log(i, 1:6) = output_input;
    
end

% Define colors
teal = [0, 0.5, 0.5];
grey = [0.5, 0.5, 0.5];

figure;
for i = 1:length(uniqueIDs)
    subplot(2, 3, i);
    if ismember(i, [2, 3, 6])
        bar(hidden_output_mean_log(i, :), 'FaceColor', teal);
    else
        bar(hidden_output_mean_log(i, :), 'FaceColor', grey);
    end
    title( uniqueIDs{i});
    xlabel('Cell');
    ylabel('Average Activation');
    ylim([-1 1]);
end

% Plotting output layer input mean 
figure;
for i = 1:length(uniqueIDs)
    subplot(2, 3, i);
    if ismember(i, [2, 3, 6])
        bar(output_input_log(i, :), 'FaceColor', teal);
    else
        bar(output_input_log(i, :), 'FaceColor', grey);
    end
    title(uniqueIDs{i});
    xlabel('Cell');
    ylabel('Average Output');
    ylim([-1 1]);
end

% input layer activaiton
figure;
for i = 1:length(uniqueIDs)
    subplot(2, 3, i);
    if ismember(i, [2, 3, 6])
        bar(hidden_input_mean_log(i, :), 'FaceColor', teal);
    else
        bar(hidden_input_mean_log(i, :), 'FaceColor', grey);
    end
    title( uniqueIDs{i});
    xlabel('Cell');
    ylabel('Average Activation');
    ylim([-1 1]);
end

%% how does number of hidden layer cells effect learning
% Train and test multilayer perceptron

[~, error, w1, w2, plateau_iter] = run_multilayer_perceptron(spikeCounts, n_hidden);
[performance_test_multi, ~] = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);
disp(performance_test_multi)

% Initialize variables
n_hidden_range = 1:8;
n_iterations = 50;  % Number of iterations to compute mean and standard deviation

mean_performances = zeros(size(n_hidden_range));
std_performances = zeros(size(n_hidden_range));

for i = 1:length(n_hidden_range)
    n_hidden = n_hidden_range(i);
    performances = zeros(1, n_iterations);
    
    for j = 1:n_iterations
        % Run the multilayer perceptron
        [~, ~, w1, w2, ~] = run_multilayer_perceptron(spikeCounts, n_hidden);
        
        % Test the multilayer perceptron
        performance_test_multi = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);
        
        % Store the performance
        performances(j) = performance_test_multi;
    end
    
    % Compute mean and standard deviation
    mean_performances(i) = mean(performances);
    std_performances(i) = std(performances);
end

% Plot the results
figure;
errorbar(n_hidden_range, mean_performances, std_performances, '-o', 'Color', 'k', 'MarkerFaceColor', 'k', 'CapSize', 5);
ylabel('Performance (Mean ± SD)');
title('Performance of Multilayer Perceptron');
xlim([min(n_hidden_range)-1, max(n_hidden_range)+1]);xlabel('Number of Hidden Layer Neurons');
xticks(n_hidden_range);
ylabel('Performance (Mean ± SD)');
title('Performance of Multilayer Perceptron');
