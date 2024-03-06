% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers

%% initial set up 
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_first = 40;
lenght_second = 40;
delay_time = 800;
start_time = 200;

% Get non-programmable paramaters
p = get_params_hipp(0.85);

% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);

% get cells to be activated by each odour in layer CA3
mems = get_odours_hipp(p, degree_overlap, pattern_order);

% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_first)];
input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
M = get_memory_hipp(p);

% simulate hippocampal dynamics 
M = simulate_dynapics_hipp(p, C, J, input, M, mems);
% plot output for a single trial 
output_plot = get_output_plot(M,pattern_order, p, mems);

%plotting mean spiking in overlapping cells vs non overlapping cells
%during second odour, filtering for cells which recived increasing number of inputs
mean_firing_second_odour = get_mean_firing_second_odour(p, C, J, input, M, mems, length_second);

% plot SVM loss function for increasing number of trials 
SVM_plot = run_classifier(p, C, J, input, M, mems);