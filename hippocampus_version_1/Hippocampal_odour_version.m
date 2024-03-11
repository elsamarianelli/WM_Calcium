% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers

%% initial set up 
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_first = 40;
lenght_second = 40;
delay_time = 600;
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
output_plot = get_output_plot(M,pattern_order, p, mems, C);

%plotting mean spiking in overlapping cells vs non overlapping cells
%during second odour, filtering for cells which recived increasing number of inputs
mean_firing_second_odour = get_mean_firing_second_odour(p, C, J, input, M, mems, length_second);

% plot SVM loss function for increasing number of trials 
SVM_plot = run_classifier(p, C, J, input, M, mems);


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

