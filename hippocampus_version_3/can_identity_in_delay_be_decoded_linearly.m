%% Code to train perceptron on odour discrimination task

%% Set parameters for the simulation
p.degree_overlap_CA3    = 0.25;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 1000;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;           % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 2500;
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
output_plot             = get_output_plot(M,p.pattern_order, p, stim, C);

%%  Simulate hippocampal dynamics  over many trials with simplest version of task (only 2 odours)
time_1              = input.simulation(2);
time_2              = input.reactivation(1);
% simulate train data
n_trials            = 100;
[spikeCounts, ~]	= get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, input.simulation(2), input.reactivation(1));
% simulate test data 
n_trials            = 20;
[spikeCounts_test, ~] = get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, input.simulation(2), input.reactivation(1));

%% run perceptron 
[performance_accuracy, error, w]    = run_perceptron_db(spikeCounts);
[performance_test]                  = test_perceptron_output(spikeCounts_test, w);

%%  Debug plot (requires fastsmooth function)
figure;
plot(fastsmooth(abs(error),50)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

%% run for different levels of background input to see how performance changes
scaleF_list = (0.84:0.002:0.85);
performance_log = zeros(1, length(scaleF_list));
for i = 1:length(scaleF_list)
    scaleF = scaleF_list(i);
    p.scaleF = scaleF;
    %%  Simulate hippocampal dynamics  over many trials with simplest version of task (only 2 odours)
    time_1              = input.simulation(2);
    time_2              = input.reactivation(1);
    % simulate train data
    n_trials            = 100;
    [spikeCounts, ~]	= get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, input.simulation(2), input.reactivation(1));
    % simulate test data 
    n_trials            = 20;
    [spikeCounts_test, ~] = get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, input.simulation(2), input.reactivation(1));
    
    %% run perceptron 
    [performance_accuracy, error, w]    = run_perceptron_db(spikeCounts);
    [performance_test]                  = test_perceptron_output(spikeCounts_test, w);
   
    % log performance 
    performance_log(scaleF) = performance_test;
end
        