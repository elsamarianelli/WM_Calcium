%% Generating different facil and delay time data
% code to generate simulate data for performance means
% and SDs for increasing delay times, for a range if increasing
% facilitation time constants (to replicate some sort of 
parent_dir = fileparts(pwd);
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% Create the Simulated_data folder if it does not exist
if ~exist(simulated_data_folder, 'dir')
    mkdir(simulated_data_folder);
end

% Set parameters for the simulation
n_trials = 100;
p.degree_overlap_CA3 = 0.2;    % Overlap between neural representations of each odour
p.degree_overlap_CA1 = 0.0;
p.start_time = 200;            % Time at which the first odour is presented (ms)
p.length_first = 250;          % Length of time for which the first odour is presented (ms)
p.length_second = 250;         % Length of time for which the second odour is presented (ms)
p.delay_time = 500;
p.scaleF = 0.848;              % Constant by which to scale random currents (to modulate baseline activity levels)
p = get_params_hipp(p);

% Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles = get_odours_hipp(p.in + (1:p.out), p.f_o, p.degree_overlap_CA1);

[C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

file_name = 'varying_tau_facil';
facils = 500:1000:4500;
facils_names = arrayfun(@(v) sprintf('100_%d_tau_facil', v), facils, 'UniformOutput', false);

% Delays variable should be defined in your workspace
delays = 250:250:3000;

for i = 1:length(facils)
    facils_name = facils_names{i};
    main_folder = fullfile(simulated_data_folder, file_name, facils_name);

    % Create the main folder if it does not exist
    if ~exist(main_folder, 'dir')
        mkdir(main_folder);
    end

    p.facil = facils(i);
    disp(['Processing facil: ', num2str(facils(i))]);

    for d = 1:length(delays)
        variable = delays(d);
        p.delay_time = delays(d);
        save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder);
        disp(['Processing delay: ', num2str(delays(d))]);
    end
end

