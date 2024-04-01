%% code to get perceptron accuracy for different delay times 
% (and plot onvergence towards 0 error for training period)

%% perceptron training on segregated CA1 poulation ouptut
%% Set parameters for the simulation
p.degree_overlap_CA3    = 0.2;              % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.2;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 40;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 600;            % Delay between odour presentations (ms)
p.length_second         = 40;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;           % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
% see_connectivity        = visualise_connectivity(C, ca3_ensembles, ca1_ensembles);

%  Specify times that each odour is presented, assign memory for the output
input.simulation        = [p.start_time p.start_time+p.length_first];
input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
M                       = get_memory_hipp(p);

%%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
n_trials = 6.*100;
data = get_train_data(C, J, input, n_trials, p, ca3_ensembles);
performance_accuracy = run_perceptron(data, n_trials, p);
