%% Code to train perceptron on odour discrimination task
%% Set parameters for the simulation
n_trials                = 200;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 700;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;           % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 1500;
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

%%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
% [spikeCounts,sequence]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
% [spikeCounts_test, ~]	= get_train_data_db(C, J, 10, p, ca3_ensembles, time3, time4);

%% Save data and settings info to current working directory
% folderName = '300_trials_fixIn_true_100_after_odour_shorter_delay';
% 
% % Create the folder if it doesn't already exist
% if ~exist(folderName, 'dir')
%     mkdir(folderName);
% end
% 
% % Save the matrix in a .mat file
% save(fullfile(folderName, 'spikeCounts.mat'), 'spikeCounts');
% save(fullfile(folderName, 'spikeCounts_test.mat'), 'spikeCounts_test');
% % Save the structure in a .mat file
% save(fullfile(folderName, 'myStruct.mat'), 'p');

%% load data to use 
folderName = '300_trials_fixIn_true_100_after_odour_shorter_delay';

% Load data from .mat files
load(fullfile(folderName, 'spikeCounts.mat'));  % Loads 'spikeCounts'
load(fullfile(folderName, 'spikeCounts_test.mat'));  % Loads 'spikeCounts_test'
load(fullfile(folderName, 'myStruct.mat'));  % Loads 'p'
%% use to train and test single layer perceptron 
[performance_accuracy_single, error, w]    = run_perceptron_db(spikeCounts(:, :));
[performance_test_single] = test_perceptron_output(spikeCounts_test, w);

% Debug plot (requires fastsmooth function)
figure;
plot(fastsmooth(abs(error),6000)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

%% use to train and test multilayer perceptron 
[~, error, w1, w2] = run_multilayer_perceptron(spikeCounts(:, :));
[performance_test_multi]  = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);

%  Debug plot (requires fastsmooth function)
figure;
plot(fastsmooth(abs(error),6000)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

%% fake data test
% assume binary input for subpopulations A B C only, and overlap
% populations, to show problem is solvable using a multilayer perceptron,
% but not a single one
target_output =  [1 1 1 0 0 0]; 
input =   [0 0 1 0 0 1
          0 1 0 0 1 0 
          1 0 0 1 0 0 
          0 1 0 1 0 0 
          0 0 1 0 1 0 
          1 0 0 0 0 1];

data = [input, target_output'];

% with single layer perceptron 
[performance_accuracy_single, error, w]    = run_perceptron_db(spikeCounts);
[performance_test_single] = test_perceptron_output(spikeCounts_test, w);

figure;
subplot(2, 1, 1)
plot(fastsmooth(abs(error),6000)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)
title('single layer performance =', performance_test_single)

hold on 

% with multilayer perceptron 
[output, error, w1, w2] = run_multilayer_perceptron(data);
[performance_test_multi] = test_multilayer_perceptron_output(data, w1, w2);

subplot(2, 1, 1)
plot(fastsmooth(abs(error),100)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)
title('multilayer performance =', performance_test_multi)