%% Code to train perceptron on odour discrimination task

%% Set parameters for the simulation
n_trials                = 100;           % Number of training trials per odour pair
p.degree_overlap_CA3    = 0.2;          % Overlap between neural representations of each odour in CA3
p.degree_overlap_CA1    = 0.2;          % Overlap between neural representations of each odour in CA1
p.start_time            = 200;          % Time at which the first odour is presented (ms)
p.length_first          = 40;           % Length of time for which the first odour is presented (ms)
p.delay_time            = 600;          % Delay between odour presentations (ms)
p.length_second         = 40;           % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;         % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

%%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
[spikeCounts,sequence]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles);
[performance_accuracy, error, w]    = run_perceptron_db(spikeCounts);

%%  Debug plot (requires fastsmooth function)
figure
plot(fastsmooth(abs(error),20)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

%% test on test data generated with same connectivity matrix 
[spikeCounts_test, ~]	= get_train_data_db(C, J, 5, p, ca3_ensembles);
[performance_test] = test_perceptron_output(spikeCounts_test, w);

%% seeing if weighting matches what woudld be expected 
% % find position of highest weighted inputs
% [sorted_numbers, sorted_indices] = sort(w, 'descend');
% biggest = sorted_indices(1:10);
% smallest = sorted_indices(end-10:end);
% disp(w(biggest)); 
% disp(w(smallest));
% 
% % CA1 cells getting the most input CA3 cells which receive overlapping
% % input in REWARD pair cases 
% 
% % overlaps..
% AB = intersect(ca3_ensembles{1}, ca3_ensembles{2});
% BC = intersect(ca3_ensembles{2}, ca3_ensembles{3});
% AC = intersect(ca3_ensembles{1}, ca3_ensembles{3});
% 
% % CA1 overlap inputs...
% AB_CA1 = find(sum(C(AB, :))>3);
% BC_CA1 = find(sum(C(BC, :))>3);
% AC_CA1 = find(sum(C(AC, :))>3);
% 
% % use overlap receiving CA1 populaitons to index weights
% AB_w = w(AB_CA1);
% BC_w = w(BC_CA1);
% AC_w = w(AC_CA1);
% 
% % see if they intersect 
% figure
% scatter(AB_CA1, ones(length(AB_CA1)), "filled"); hold on; 
% scatter(BC_CA1, ones(length(BC_CA1)), "filled"); hold on; 
% scatter(AC_CA1, ones(length(AC_CA1)), "filled"); hold on; 
% scatter(biggest, 2.*ones(size(biggest)), "filled");  hold on;
% % scatter(smallest, 2.*ones(size(smallest)), "filled"); 
% ylim([0 3])
