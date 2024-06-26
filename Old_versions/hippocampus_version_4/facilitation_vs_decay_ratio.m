%% generates data moving from the standard tau facil and tau decay times of 200 
%% and 1500 ms respectively, to 1500 and 200 ms
% aka moving from having facilitatory synapses to depressing synapses, this
% should generate a u curve for performance as even when synapses are
% depressing the difference in firing rate between depressed and non
% depressed synpases will hold?

%% Generate data

% Set parameters for the simulation
n_trials                = 100;
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.delay_time            = 3000;            % delay time
p.length_first          = 250;            % Length of time for which the first odour is presented (ms)
p.length_second         = 250;            % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.848;          % Constant by which to scale random currents (to modulate baseline activity levels)
p                       = get_params_hipp(p);

%  Randomly assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
    
% Create the folder if it doesn't already exist
main_folder = 'FixIn_true_delay_500_CA3overlap_0.2_trials_100_1st_100_secs_u_varied_x_constant';

if ~exist(main_folder, 'dir')
    mkdir(main_folder);
end

% Define the range of u and x
% tau_facils          = (200:100:1500);
tau_facils          = (600:200:2000);
% tau_decays          = flip(tau_facils);
tau_decays          = repmat(200, 1, length(tau_facils));           

for index = 1:length(tau_facils)
    p.tau_facil = tau_facils(index);
    p.tau_decay = tau_decays(index);
    disp([p.tau_facil p.tau_decay])                        
    save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, p.tau_facil, main_folder);
end

%% get plots 
% get correct variable names
a = arrayfun(@num2str, [tau_facils; tau_decays], 'UniformOutput', false);
variable_name = strcat(a(1, :), '_', a(2, :));

% run perceptron and plot
cmap = spring(64); 

facil_decay_performance = plot_performance_across_variable_change(variable_name, main_folder, cmap, 'synaptic variables');
facil_decay_perceptron_training = plot_examples_of_perceptron_training(variable_name, main_folder, cmap, 'synaptic variables');


% %% adding second axis in with decay and facil/decay line
% figure(facil_decay_performance)
% 
% hold on;
% % line to show where synapses cross from being facilitatory to depressing
% xline(7.5, '--k'); 
% 
% % or alternatively shade in
% % hold on;
% % x_fill = [0, 7.5, 7.5, 0];
% % y_limits = ylim;  
% % y_fill = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
% % fill(x_fill, y_fill, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
% 
% % Get the current axes from the figure, assuming there is only one
% ax1 = findobj(facil_decay_performance, 'Type', 'axes');
% 
% % Create the second axes that shares the x-axis with the original
% ax2 = axes('Position', ax1.Position, ...
%            'XAxisLocation', 'top', ...
%            'YAxisLocation', 'right', ...
%            'Color', 'none', ... % Make the axes transparent
%            'YTick', [], ... % No y-ticks on the second axes
%            'XColor', 'k', ... % Set color of x-axis (optional)
%            'Box', 'off'); % Turn off the box to avoid double line on the right
% 
% % Link the x-axes of the two axes
% linkaxes([ax1, ax2], 'x');
% 
% % rename axis to tau facil and tau decay
% % Define `extraXTicks` as needed for your plot
% % Set additional x-ticks for ax2 and label them
% extraXTicks = 1:length(tau_decays); 
% ax2.XTick = extraXTicks; 
% ax2.XTickLabel = tau_decays;
% ax2.FontSize = 10; % Set font size of tick labels to 10
% 
% % Set additional x-ticks for ax1 and label them
% extraXTicks = 1:length(tau_facils); 
% ax1.XTick = extraXTicks; 
% ax1.XTickLabel = tau_facils;
% ax1.FontSize = 10; % Set font size of tick labels to 10% relable axis 
% 
% xlabel(ax1, 'tau facil', FontSize=12);
% xlabel(ax2, 'tau decay', FontSize=12);
% title(ax1, '');
% 
% xlim([0 15])
