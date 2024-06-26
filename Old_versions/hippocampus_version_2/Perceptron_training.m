%% code to get perceptron accuracy for different delay times 
% (and plot onvergence towards 0 error for training period)

%% perceptron training on segregated CA1 poulation ouptut
%% Set parameters for the simulation
p.degree_overlap_CA3  = 0.2;          % Overlap between neural representations of each odour
p.degree_overlap_CA1  = 0.0;
p.start_time          = 100;          % Time at which the first odour is presented (ms)
p.length_first        = 40;           % Length of time for which the first odour is presented (ms)
p.delay_time          = 700;          % Delay between odour presentations (ms)
p.length_second       = 40;           % Length of time for which the second odour is presented (ms)
p.scaleF              = 0.85;         % Constant by which to scale random currents (to modulate baseline activity levels)
p                     = get_params_hipp(p);

% do we want to constrain the overlap of CA1 populations which receive
% input from corresponding CA3 populaitons? or have projections be
% completely random

%  Assign CA3 cells to each odour representation
CA3_populations = get_odours_hipp(p, "OFF");
CA1_populations = get_odours_hipp(p, "ON");

%  Generate connectivity and synaptic efficacy matrix
p.c = 0.2;
overlap_control     = "ON";
[C, J]              = connectivity_matrix_hipp(p, overlap_control, CA3_populations, CA1_populations);

%  Specify times that each odour is presented, assign memory for the output
input.simulation    = [p.start_time p.start_time+p.length_first];
input.reactivation  = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
M                   = get_memory_hipp(p);
M = simulate_dynamics_hipp(p, C, J, input, M, mems);

% see U and X over time in synapses whihc receive overlapping input and
% those which only fire with odour 2, also firing across CA3 and CA1
% populations
output_plot         = get_output_plot(M,p.pattern_order, p, mems, C);
% visualise connectivity in part of the network, with different odour
% corresponding populations coloured differently, and connections for that
% portion of the network shown
see_connectivity    = visualise_connectivity(C, CA3_populations, CA1_populations);

n_trials = 6.*200;
data = get_train_data(C, J, input, n_trials, p, CA3_populations);

performance_accuracy = run_perceptron(data, n_trials, p);




%% generate training data for incremental increases in connectivity level between CA1 and CA3 --> needs updating still from v1
n_trials = 6.*200;
% Define the size of the cell array
num_iterations = 15;
connectivity_trial_data = cell(2, num_iterations);

for i = 1:5:num_iterations
    connectivity_level = 0.02.*i;
    p.c = connectivity_level;
    % simulate training + test data for 0.2 overlap
    data = get_train_data(C, J, input, n_trials, 0.2, p);
    connectivity_trial_data{1, i} = data;
    connectivity_trial_data{2, i} = connectivity_level;
    disp(i)
end


data = connectivity_trial_data{:, 1};
% save for later
save("Perceptron_performance_data_connectivity.mat", "delay_trial_data")

%% generate training data for incremental delay times
n_trials = 6.*200;
% Define the size of the cell array
num_iterations = 20;
delay_trial_data = cell(2, num_iterations);

for i = 1:num_iterations

    % Variable Times that memory is 'on' (ms) in each loop
    delay_time = 100.*i;
    input.simulation = [start_time (start_time+length_first)];
    input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];

    % simulate training + test data for 0.2 overlap
    data = get_train_data(C, J, input, n_trials, 0.2, p);
    delay_trial_data{1, i} = data;
    delay_trial_data{2, i} = delay_time;
    disp(i)

end

% save for later
% save("Perceptron_performance_data.mat", "delay_trial_data")

% train and test perceptron for each round of simulated data
for i = 1:5:15

    delay_trial = delay_trial_data(:, 1:end);
    performance_accuracy_all = zeros(size(delay_trial, 2));

for i=1:size(delay_trial, 2)
    data = delay_trial{1, i};
    delay = delay_trial{2, i};
    performance_accuracy = run_perceptron(data, n_trials, p);
    delay_trial{3, i} = performance_accuracy;
end

% plot variable delays
plot(1:22, cell2mat(delay_trial(3, :)))
custom_ticks = 1:22;
custom_labels = (1:22).*100; 
xticks(custom_ticks);
xticklabels(custom_labels);
xlabel('delay time between first and second odour')
ylabel('performance of perceptron on test data after training')
box off

