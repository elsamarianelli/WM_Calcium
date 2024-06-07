%% plotting MLP figures 

% get pathway to parent simulated data 
parent_dir = fileparts(pwd);
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% [1] Background excitation strength (p.CF varied)
file_name = 'CF_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.01 : 1.15;
CF_performance = plot_performance_across_variable_change(CFs, main_folder, 'Background excitation strength');

% [2] Odour strength (p.SF varied)
file_name = 'SF_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
SF_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [3] Delay time between odours (p.delay_time varied)
file_name = 'delay_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
delay_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [4] STSF vs STSD (p.tau_decay and p.tau_facil varied) 
file_name = 'ux_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
ux_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [5] Connectivity (p.c varied)
file_name = 'ux_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
connectivity_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [6] CA1 overlap (p.degree_overlap_CA1 varied) 
file_name = 'ux_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
overlap_CA1_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [7] CA3 overlap (p.degree_overlap_CA3 varied) 
file_name = 'ux_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.05 : 1.15;
overlap_CA3_performance = plot_performance_across_variable_change(CFs, main_folder, 'delay');

% [8] Just calcium decay longer (p.tau_facil, p.tau_decay fixed at 100ms)

% needs to be done a bit differently 
