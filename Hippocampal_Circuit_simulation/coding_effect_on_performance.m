%% what is the effect of CA3 and CA1 odour overlap and connectivity on performance 
% (overlap of odour sensitive subpopulations in both layers) and
% connectivity between the two layers) on multilayer perceptron
% performance?

% Get the path to the parent directory of the current folder
parent_dir = fileparts(pwd);

%% Connectivity level (P.c aka probability of synaptic contact) effect on perceptron performance

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.delay_time            = 500;
p.scaleF                = 0.848;           % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

% Create the folder if it doesn't already exist
% Define the main folder name
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_connectivity_varied';

% Construct the path to the Simulated_data folder
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% Construct the path to the main folder within the Simulated_data folder
main_folder_path = fullfile(simulated_data_folder, main_folder);

% Check if the main folder exists, if not, create it
if ~exist(main_folder_path, 'dir')
    mkdir(main_folder_path);
end

% Define the range of connectivity levels

variable_list = 0:0.025:0.4;
for index = 1:length(variable_list)

    % change p to correct current variable 
    variable = variable_list(index);
    disp(variable)
    p.c = variable;

    %  Generate connectivity and synaptic efficacy matrix
    [C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
    
    %save data for training and testing of perceptron
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder_path)

end

% plot perceptron performance trained on this data
performance_variable_figure = plot_performance_across_variable_change(variable_list, main_folder_path, 'connectivity');

% Define the file name and format
fileName = 'connectivity_plot';
fileFormat = 'fig';

% Construct the path to the figure folder within the parent directory
figure_folder = fullfile(parent_dir, 'Figures');

% Construct the full path for the file to be saved
full_file_path = fullfile(figure_folder, [fileName, '.', fileFormat]);

% Save the figure to the specified folder
saveas(performance_variable_figure, full_file_path);

%% overlap of CA3 layer

%reset params
p                       = get_params_hipp(p);

% Create the folder if it doesn't already exist
% Define the main folder name
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_CA3_overlap_varied';

% Construct the path to the Simulated_data folder
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% Construct the path to the main folder within the Simulated_data folder
main_folder_path = fullfile(simulated_data_folder, main_folder);

% Check if the main folder exists, if not, create it
if ~exist(main_folder_path, 'dir')
    mkdir(main_folder_path);
end

% Define the range of connectivity levels
p.degree_overlap_CA1 = 0;

variable_list = 0:0.05:0.3;
for index = 1:length(variable_list)

    % change p to correct current variable 
    variable = variable_list(index);
    disp(variable)
    p.degree_overlap_CA3 = variable;

    %  Randomly assign CA3 and CA1 cells to each odour representation
    ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
    ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

    %  Generate connectivity and synaptic efficacy matrix
    [C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

    %save data for training and testing of perceptron
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder_path)

end

% plot perceptron performance trained on this data
performance_variable_figure = plot_performance_across_variable_change(variable_list, main_folder_path, 'CA3 overlap');

% Define the file name and format
fileName = 'CA3_overlap';
fileFormat = 'fig';

% Construct the path to the figure folder within the parent directory
figure_folder = fullfile(parent_dir, 'Figures');

% Construct the full path for the file to be saved
full_file_path = fullfile(figure_folder, [fileName, '.', fileFormat]);

% Save the figure to the specified folder
saveas(performance_variable_figure, full_file_path);

%% overlap of CA1 layer
p                       = get_params_hipp(p);

% Create the folder if it doesn't already exist
% Define the main folder name
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_CA1_overlap_varied';

% Construct the path to the Simulated_data folder
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% Construct the path to the main folder within the Simulated_data folder
main_folder_path = fullfile(simulated_data_folder, main_folder);

% Check if the main folder exists, if not, create it
if ~exist(main_folder_path, 'dir')
    mkdir(main_folder_path);
end

% Define the range of connectivity levels
p.degree_overlap_CA3 = 0;

variable_list = 0:0.05:0.3;
for index = 1:length(variable_list)

    % change p to correct current variable 
    variable = variable_list(index);
    disp(variable)
    p.degree_overlap_CA1 = variable;

    %  Randomly assign CA3 and CA1 cells to each odour representation
    ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
    ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

    %  Generate connectivity and synaptic efficacy matrix
    [C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

    %save data for training and testing of perceptron
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder_path)

end

% plot perceptron performance trained on this data
performance_variable_figure = plot_performance_across_variable_change(variable_list, main_folder_path, 'CA1 overlap');

% Define the file name and format
fileName = 'CA1_overlap';
fileFormat = 'fig';

% Construct the path to the figure folder within the parent directory
figure_folder = fullfile(parent_dir, 'Figures');

% Construct the full path for the file to be saved
full_file_path = fullfile(figure_folder, [fileName, '.', fileFormat]);

% Save the figure to the specified folder
saveas(performance_variable_figure, full_file_path);

