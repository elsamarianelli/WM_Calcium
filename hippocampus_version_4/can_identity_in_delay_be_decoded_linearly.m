%% Code to train perceptron on odour discrimination task
% Finding a level of background noise that gives 'silent' delay period 
% (no linear decoding of previous odour identity) and 'active' delay period
% (significantly elevated CA1 rates, linear decoding is possible)

%% Set parameters for the simulation
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 700;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;           % Constant by which to scale random currents (to modulate baseline activity levels)
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
output_plot             = get_output_plot(M,p);

%%  Simulate hippocampal dynamics  over many trials with simplest version of task (only 2 odours)
time_1              = input.simulation(2);
time_2              = input.reactivation(1);
  
% simulate train data
n_trials            = 100;
[spikeCounts, ~]	= get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, time_1, time_2, input);
% simulate test data 
n_trials            = 10;
[spikeCounts_test, ~] = get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, time_1, time_2, input);

% run perceptron 
[performance_accuracy, error, w]    = run_perceptron_db(spikeCounts);
[performance_test]                  = test_perceptron_output(spikeCounts_test, w);
disp(performance_test)

%  Debug plot (requires fastsmooth function)
figure;
plot(fastsmooth(abs(error),100)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

%% run for different levels of background input to see how performance changes

% scaleF_list = (0.82:0.002:0.87);
scaleF_list = (0.85:0.002:0.86);
performance_log = zeros(length(scaleF_list), 1:20);

% Create the folder for data storage it doesn't already exist
folderName = 'Linear_decoding_delay_period_focus_more_trials';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

for i = 1:length(scaleF_list)

    scaleF = scaleF_list(i);
    p.scaleF = scaleF;
    
    %%  Simulate hippocampal dynamics  over many trials with simplest version of task (only 2 odours)
    time_1              = input.simulation(2);
    time_2              = input.reactivation(1);
    
    % simulate train data
    n_trials            = 100;
    [spikeCounts, ~]	= get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, time_1, time_2, input);
    % simulate test data 
    n_trials            = 10;
    [spikeCounts_test, ~] = get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, time_1, time_2, input);

  
    % log different performance runs (can average later)
    for performance_run = 1:20
          %% run perceptron and log performance
        [performance_accuracy, error, w]    = run_perceptron_db(spikeCounts);
        [performance_test]                  = test_perceptron_output(spikeCounts_test, w);
        disp(performance_test)
        performance_log(i, performance_run) = performance_test;
    end

    % Save the trin and test data for each contrast factor 
    sub_folderName = ['Contrast_factor_', num2str(scaleF), '.mat'];

    if ~exist(fullfile(folderName, sub_folderName), 'dir')
        mkdir(fullfile(folderName, sub_folderName));
    end
    
    performance = performance_log(i, 1:20);
    performance_mean = mean(performance);

    save(fullfile(folderName, sub_folderName, 'train_data.mat'), 'spikeCounts');
    save(fullfile(folderName, sub_folderName, 'test_data.mat'), 'spikeCounts_test');
    save(fullfile(folderName, sub_folderName, 'performance.mat'), 'performance');
    save(fullfile(folderName, sub_folderName, 'performance_mean.mat'), 'performance_mean');
end


%% plot perceptron performance across delay times

folderName = 'Linear_decoding_delay_period';
scaleF_list = (0.82:0.002:0.87);  % Ensure this matches the scaleF_list used during saving
performance_log_reloaded = zeros(1, length(scaleF_list));

for i = 1:length(scaleF_list)
    scaleF = scaleF_list(i);
    sub_folderName = ['Contrast_factor_', num2str(scaleF), '.mat'];
    
    % Construct the full path to the performance.mat file
    performance_file = fullfile(folderName, sub_folderName, 'performance.mat');
    
    % Load the performance data from the file
    data = load(performance_file, 'performance_test');
    
    % Store the loaded performance data into an array
    performance_log_reloaded(i) = data.performance_test;
end

% Display the loaded performance data
disp('Reloaded performance values:');
disp(performance_log_reloaded);

plot(scaleF_list, performance_log_reloaded); 
ylim([0 1])