% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers
 
% Elsa Marianelli, UCL (2024)

%% Set parameters for the simulation
p.degree_overlap_CA3    = 0;              % Overlap between neural representations of each odour
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

% take first and second odour to be presented 
stim                    = cell(2,1); 
first                   = double(upper(p.pattern_order(1))) - 64; 
stim{1}                 = ca3_ensembles{first}; clear first
second                  = double(upper(p.pattern_order(2))) - 64; 
stim{2}                 = ca3_ensembles{second}; clear second

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
see_connectivity        = visualise_connectivity(C, ca3_ensembles, ca1_ensembles);

%  Specify times that each odour is presented, assign memory for the output
input.simulation        = [p.start_time p.start_time+p.length_first];
input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
M                       = get_memory_hipp(p);

%%  Simulate hippocampal dynamics 
M                       = simulate_dynamics_hipp(p, C, J, input, M, stim);

%  Plot output for a single trial 
output_plot             = get_output_plot(M,p.pattern_order, p, stim, C);

%% 
% % plotting mean spiking in overlapping cells vs non overlapping cells
% % during second odour, filtering for cells which recived increasing number of inputs
% mean_firing_second_odour = get_mean_firing_second_odour(p, C, J, input, M, mems, p.length_second);

%% keep for later but this isn't being used now
% plot SVM loss function for increasing number of trials 
% SVM_plot = run_classifier(p, C, J, input, M, mems);
% %% training perceptron with matlab function 
% % create perceptron with CA1 mean firing during odour 2 for each cell as 
% % input, and a single "lick" output neuron 
% n_trials = 6.*50;
% train_data_overlap = get_train_data(C, J, input, n_trials, degree_overlap, p);
% 
% train_data = train_data_no_overlap;
% P = train_data(1:n_trials, 1:p.out)'; 
% T = train_data(1:n_trials, end)';
% shuffle_T = T(randperm(length(T)));
% 
% % single layer notshuffled
% net = perceptron;
% net.trainParam.epochs =100;
% [net, tr] = train(net,P, T);
% 
% % single layer shuffled
% net = perceptron;
% net.trainParam.epochs =100;
% [net, tr] = train(net,P, shuffle_T);
% 
% % %evaluate weights and biases
% % w = net.iw{1,1}; b = net.b{1};
% 
% % multilayer not shuffled 
% net = feedforwardnet(1);
% [net,tr] = train(net,P, T);
% y = net(P);
% perf = perform(net,y,t);
% 
% % multilayer shuffled
% net = feedforwardnet(1);
% [net,tr] = train(net,P, shuffle_T);
% y = net(P);
% perf = perform(net,y,t);
% 
% 
% % with pattern net 
% net = patternnet(1);
% net = train(net,P, T);
% view(net)