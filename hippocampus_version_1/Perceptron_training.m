%% code to run perceptron with an online learning algorithm 
%  MATLAB function training was being done using batch learning (aka updates
%  weights with cumulative errors after each epoch) 
%  This is to do online learning (aka updating weights after each trial
%  presentation)
%  and just to have a better idea of what's ouput in general.

%% get input
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_first = 30;
lenght_second = 30;
delay_time = 100;
start_time = 200;
% Get non-programmable paramaters
p = get_params_hipp(0.85);
% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);
% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_first)];
input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
% generate memory
M = get_memory_hipp(p);

%% generate training data for incremental delay times
n_trials = 6.*200;
% Define the size of the cell array
num_iterations = 20;
delay_trial_data = cell(2, num_iterations);

for i = 22:num_iterations
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
delay_trial = delay_trial_data(:, 1:22);
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

