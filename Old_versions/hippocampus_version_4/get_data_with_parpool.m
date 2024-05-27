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

%% generate data and plot contrast factor affect on performance 

% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_CF_varied';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of delays
variable_list = 1 : 0.01 : 1.15;

for index = 1:length(variable_list)

    variable = variable_list(index);
    disp(variable)

    % this is the line that needs to be changed to check different
    % parameters (and contents of variable list)
    p.CF = variable;
   
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder)
    
end

% plotting perceptron training and mean performance with varying scale factors  
plot_handle = plot_performance_across_variable_change(variable_list, main_folder, 'Contrast factor');

% saving fig
fileName = 'CF_varied';
fileFormat = 'fig'; 
saveas(plot_handle, fileName, fileFormat);

%% generate data and plot CA3 overlap affect on performance 
% Set parameters for the simulation
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_CA3_overlap_varied';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of delays
variable_list = 0 : 0.04 : 0.4;

for index = 1:length(variable_list)

    variable = variable_list(index);
    disp(variable)

    % this is the line that needs to be changed to check different
    % parameters (and contents of variable list)
    p.degree_overlap_CA1 = variable;
   
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder)
    
end

% plotting perceptron training and mean performance with varying scale factors  
plot_handle = plot_performance_across_variable_change(variable_list, main_folder, 'CA3 overlap');

% saving fig
fileName = 'CA3_overlap';
fileFormat = 'fig'; 
saveas(plot_handle, fileName, fileFormat);

%% generate data and plot CA3 overlap affect on performance 
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p                       = get_params_hipp(p);

% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_connectivity_varied';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of delays
variable_list = 0.05 : 0.01 : 0.15;

for index = 1:length(variable_list)

    variable = variable_list(index);
    disp(variable)

    % this is the line that needs to be changed to check different
    % parameters (and contents of variable list)
    p.c = variable;
   
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder)
    
end

% plotting perceptron training and mean performance with varying scale factors  
plot_handle = plot_performance_across_variable_change(variable_list, main_folder, 'connectivity');

% saving fig
fileName = 'connectivity';
fileFormat = 'fig'; 
saveas(plot_handle, fileName, fileFormat);

%% generate data and plot delay at different tau decays
% Set parameters for the simulation
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

variable_list_big = 500 : 500 : 3000;

for index = 1:length(variable_list_big)

    variable = variable_list_big(index);

    p.tau_facil = variable;
    % Create the folder if it doesn't already exist
    main_folder = ['FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_' num2str(variable) '_tau_facil'];
    
    if ~exist(main_folder, 'dir')
        mkdir(main_folder);
    end
    
    % Define the range of delays
    variable_list = 250 : 250 : 3000;
    
    for index = 1:length(variable_list)
    
        variable = variable_list(index);
        disp(variable)
    
        % this is the line that needs to be changed to check different
        % parameters (and contents of variable list)
        p.degree_overlap_CA1 = variable;
       
        save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder)
        
    end

end

% plotting perceptron training and mean performance with varying scale factors  
plot_handle = plot_performance_across_variable_change(variable_list, main_folder, 'CA3 overlap');

% saving fig
fileName = 'CA3_overlap';
fileFormat = 'fig'; 
saveas(plot_handle, fileName, fileFormat);