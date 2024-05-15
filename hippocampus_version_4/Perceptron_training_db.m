%% Code to train perceptron on odour discrimination task

%% check dynamics over a single trial
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 700;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.84;           % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 1500;
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

%% looking at performance for a single setting
time3 = input.reactivation(1);
time4 = time3+100; 

% %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
[spikeCounts,~]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
[spikeCounts_test, ~]	= get_train_data_db(C, J, 20, p, ca3_ensembles, time3, time4);
% Train and test multilayer perceptron
[~, error, w1, w2] = run_multilayer_perceptron(spikeCounts);
[performance_test_multi] = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);
disp(performance_test_multi)

%% Save data and settings info to current working directory
folderName = 'fixIn_true_100_after_odour_2000_delay';

% Create the folder if it doesn't already exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Save ouput to file
save(fullfile(folderName, 'spikeCounts.mat'), 'spikeCounts');
save(fullfile(folderName, 'spikeCounts_test.mat'), 'spikeCounts_test');
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
    
% % Create the folder if it doesn't already exist
% main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs';
% 
% if ~exist(main_folder, 'dir')
%     mkdir(main_folder);
% end
% 
% % Define the range of delays
% delays = 0:250:3000;
% parfor index = 2:length(delays)
%     delay_time = delays(index);
%     disp(delay_time)
%     save_data_varying_delay_times(C, J, p, ca3_ensembles, delay_time, n_trials, main_folder);
% end

%% plotting perceptron training and mean performance with varying delay times 

delays = [250:250:2500 3000];
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs';
cmap = winter(64); 

delay_performance = plot_performance_across_variable_change(delays, main_folder, cmap, 'delay');
% delay_perceptron_training = plot_examples_of_perceptron_training(delays, main_folder, cmap, 'delay');

%% plotting perceptron training and mean performance with varying scale factors  
scaleFs = (0.84:0.002:0.858);

main_folder = 'FixIn_true_delay_2000_CA3overlap_0.2_trials_100_1st_100_secs_CF_varied';
cmap = spring(64); 

cf_performance = plot_performance_across_variable_change(scaleFs, main_folder, cmap, 'contrast factor');
% cf_perceptron_training = plot_examples_of_perceptron_training(scaleFs, main_folder, cmap, 'contrast factor');

