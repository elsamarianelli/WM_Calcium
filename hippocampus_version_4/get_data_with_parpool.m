%% generating data using parpool

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.delay_time            = 1000;           % delay time
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;           % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
    
% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_delay_2000_CA3overlap_0.2_trials_100_1st_100_secs_CF_varied';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of contrast factors
scaleF_list = (0.84:0.002:0.86);

parfor index = 1:length(scaleF_list)
    scaleF = scaleF_list(index);
    disp(scaleF)                        
    save_data_varying_delay_times(C, J, p, ca3_ensembles, p.delay_time, scaleF, n_trials, main_folder);
end