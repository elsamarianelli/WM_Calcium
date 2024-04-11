% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers
 
% Elsa Marianelli, UCL (2024)

%% Set parameters for the simulation
n_trials                = 100;           % Number of training trials per odour pair
p.degree_overlap_CA3    = 0.2;          % Overlap between neural representations of each odour in CA3
p.degree_overlap_CA1    = 0.0;          % Overlap between neural representations of each odour in CA1
p.start_time            = 200;          % Time at which the first odour is presented (ms)
p.length_first          = 100;           % Length of time for which the first odour is presented (ms)
p.delay_time            = 1000;          % Delay between odour presentations (ms)
p.length_second         = 100;           % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;         % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength = 1500;                      % Length of simulation (ms)
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