%% Code to train perceptron on odour discrimination task

%% check dynamics over a single trial
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.2;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;            % Length of time for which the first odour is presented (ms)
p.delay_time            = 500;            % Delay between odour presentations (ms)
p.length_second         = 250;            % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;          % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 1500;
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);

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
output_plot             = get_output_plot(M,p,ca3_ensembles, C);

%% looking at performance for a single setting
time3 = input.reactivation(1);
time4 = time3+100; 
% 
% %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
% %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
n_trials = 100;
[spikeCounts,~]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
[spikeCounts_test, ~]	= get_train_data_db(C, J, 10, p, ca3_ensembles, time3, time4);

% Train and test multilayer perceptron
[~, error, w1, w2, plateau_iter] = run_multilayer_perceptron(spikeCounts);
[performance_test_multi] = test_multilayer_perceptron_output(spikeCounts_test, w1, w2);
disp(performance_test_multi)

%% Save data and settings info to current working directory
folderName = 'fixIn_true_100_after_odour_2000_delay';

% Create the folder if it doesn't already exist
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Save ouput to file
save(fullfile(folderName, 'spikeCounts.mat'), 'spikeCounts');
save(fullfile(folderName, 'spikeCounts_test.mat'), 'spikeCounts_test');
save(fullfile(folderName, 'myStruct.mat'), 'p');

%% creating delay data using parpool

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.delay_time            = 500;
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

%% plotting perceptron training and mean performance with varying delay times 

delays = [250:250:2500 3000];
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_delay_varied';

delay_performance = plot_performance_across_variable_change(delays, main_folder, 'delay');

% saving fig
fileName = 'delay_plot';
fileFormat = 'fig'; 
saveas(delay_performance, fileName, fileFormat);

%% plotting perceptron training and mean performance with varying delay times where trained perceptron is kept CONSTANT and only testig data changes delay length
parent_dir = fileparts(pwd);
main_folder = 'FixIn_true_CF_0.848_CA3overlap_0.2_trials_100_1st_100_secs_delay_varied';

% Construct the path to the Simulated_data folder
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');

% Construct the path to the main folder within the Simulated_data folder
main_folder_path = fullfile(simulated_data_folder, main_folder);

delays = [250:250:2500 3000];

delay_performance = plot_performance_across_variable_change_pretrained_perceptron(delays, main_folder_path, 'delay');
% Switch to the right y-axis
yyaxis right

% Set the limits of the second y-axis
ylim([-500 2000]); % Replace ymin and ymax with your desired limits
% saving fig
fileName = 'delay_plot_fixed_training';
fileFormat = 'fig'; 
saveas(delay_performance, fileName, fileFormat);

%% plotting perceptron training and mean performance with varying scale factors  
scaleFs = (0.82:0.005:0.875);

main_folder = 'FixIn_true_delay_2000_CA3overlap_0.2_trials_100_1st_100_secs_SF_varied';

sf_performance = plot_performance_across_variable_change(scaleFs, main_folder, 'scale factor');

% saving fig
fileName = 'scale_factor_performance_plot';
fileFormat = 'fig'; 
saveas(sf_performance, fileName, fileFormat);

%% plotting perceptron training and mean performance with varying synaptic variables
tau_facils          = (200:100:1500);
tau_decays          = flip(tau_facils);

a = arrayfun(@num2str, [tau_facils; tau_decays], 'UniformOutput', false);
variable_range = strcat(a(1, :), '_', a(2, :));

main_folder = 'FixIn_true_delay_500_CA3overlap_0.2_trials_100_1st_100_secs_ux_varied';

facil_decay_performance = plot_performance_across_variable_change(variable_range, main_folder, 'synaptic variables');

%% adding second axis in with decay and facil/decay line
figure(facil_decay_performance)

hold on;
% line to show where synapses cross from being facilitatory to depressing
xline(7.5, '--k'); 

% or alternatively shade in
% hold on;
% x_fill = [0, 7.5, 7.5, 0];
% y_limits = ylim;  
% y_fill = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
% fill(x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Get the current axes from the figure, assuming there is only one
ax1 = findobj(facil_decay_performance, 'Type', 'axes');

% Create the second axes that shares the x-axis with the original
ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ... % Make the axes transparent
           'YTick', [], ... % No y-ticks on the second axes
           'XColor', 'k', ... % Set color of x-axis (optional)
           'Box', 'off'); % Turn off the box to avoid double line on the right

% Link the x-axes of the two axes
linkaxes([ax1, ax2], 'x');

% rename axis to tau facil and tau decay
% Define `extraXTicks` as needed for your plot
% Set additional x-ticks for ax2 and label them
extraXTicks = 1:length(tau_decays)+1; 
ax2.XTick = extraXTicks; 
ax2.XTickLabel = tau_decays;
ax2.FontSize = 10; % Set font size of tick labels to 10

% Set additional x-ticks for ax1 and label them
extraXTicks = 1:length(tau_facils)+1; 
ax1.XTick = extraXTicks; 
ax1.XTickLabel = tau_facils;
ax1.FontSize = 10; % Set font size of tick labels to 10% relable axis 

xlabel(ax1, 'tau facil', FontSize=12);
xlabel(ax2, 'tau decay', FontSize=12);
title(ax1, '');

xlim([0 15])

% saving fig
fileName = 'ux_varied_plot';
fileFormat = 'fig'; 
saveas(facil_decay_performance, fileName, fileFormat);

