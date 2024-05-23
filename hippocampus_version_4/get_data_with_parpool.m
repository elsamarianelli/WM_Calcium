%% generating data using parpool

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.delay_time            = 500;           % delay time
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
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_CF_varied';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of delays
CFs = 1 : 0.05 : 1.15;
for index = 1:length(CFs)
    cf = CFs(index);
    disp(cf)
    p.CF = cf;
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, 'contrast_factor', main_folder)
end

%% plotting perceptron training and mean performance with varying scale factors  

delay_performance = plot_performance_across_variable_change(cfs, main_folder, 'Contrast_factor');

% saving fig
fileName = 'CF_varied';
fileFormat = 'fig'; 
saveas(delay_performance, fileName, fileFormat);
