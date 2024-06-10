%% plotting MLP figures 

% get pathways to relevant folders
parent_dir = fileparts(pwd);
simulated_data_folder = fullfile(parent_dir, 'Simulated_data');
figure_folder = fullfile(parent_dir, 'figures_new');
fileFormat = 'fig'; 

% [1] Background excitation strength (p.CF varied)
file_name = 'CF_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CFs = 1 : 0.01 : 1.08;
CF_performance = plot_performance_across_variable_change(CFs, main_folder, 'Background excitation strength');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(CF_performance, full_file_path);

% [2] Odour strength (p.SF varied)
file_name = 'SF_varied';
main_folder = fullfile(simulated_data_folder, file_name);
SFs = 0.82:0.005:0.88;
SF_performance = plot_performance_across_variable_change(SFs, main_folder, 'Odour strength');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(SF_performance, full_file_path);

% [3] Delay time between odours (p.delay_time varied)
file_name = 'delay_varied';
main_folder = fullfile(simulated_data_folder, file_name);
delay_times = [250:250:2500 3000];
delay_performance = plot_performance_across_variable_change(delay_times, main_folder, ' Delay time');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(delay_performance, full_file_path);

% [4] STSF vs STSD (p.tau_decay and p.tau_facil varied) 
% do tau decay manually for now untill i manage to fix it)
file_name = 'ux_varied';
main_folder = fullfile(simulated_data_folder, file_name);

tau_facils = 200:100:1500;
tau_decays = flip(tau_facils);
a = arrayfun(@num2str, [tau_facils; tau_decays], 'UniformOutput', false);
uxs = strcat(a(1, :), '_', a(2, :));

facil_decay_performance = plot_performance_across_variable_change(variable_range, main_folder, 'synaptic variables');
% Adding second axis with decay and facil/decay line
figure(facil_decay_performance)
hold on;

% Line to show where synapses cross from being facilitatory to depressing
xline(7.5, '--r'); 

% Get the current axes from the figure
ax1 = findobj(facil_decay_performance, 'Type', 'axes');
ax1.XTickLabel = arrayfun(@num2str, tau_facils, 'UniformOutput', false);
saveas(facil_decay_performance, fileName, fileFormat);

% [5] Connectivity (p.c varied)
file_name = 'connectivity_varied';
main_folder = fullfile(simulated_data_folder, file_name);
cs = 0 : 0.025 : 0.175;
connectivity_performance = plot_performance_across_variable_change(cs, main_folder, 'connectivity');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(connectivity_performance, full_file_path);

% [6] CA1 overlap (p.degree_overlap_CA1 varied) 
file_name = 'CA1_overlap_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CA1s = 0 : 0.05: 0.25;
overlap_CA1_performance = plot_performance_across_variable_change(CA1s, main_folder, 'CA1 overlap');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(overlap_CA1_performance, full_file_path);

% [7] CA3 overlap (p.degree_overlap_CA3 varied) 
file_name = 'CA3_overlap_varied';
main_folder = fullfile(simulated_data_folder, file_name);
CA3s = 0 : 0.05: 0.25;
overlap_CA3_performance = plot_performance_across_variable_change(CA3s, main_folder, 'CA3 overlap');

full_file_path = fullfile(figure_folder, [file_name, '.', fileFormat]);
saveas(overlap_CA3_performance, full_file_path);

% % [8] Just calcium decay longer (p.tau_facil, p.tau_decay fixed at 100ms)
% 
% % plot the initial facl varied plot using function
% % need to run this data again 
% file_name = 'varying_tau_facil';
% 
% facils = 500:500:3000;
% facils = arrayfun(@(v) sprintf('100_%d_tau_facil', v), facils, 'UniformOutput', false);
% 
% delays = 250:250:3000;
% 
% facil = facils(1); facil = facil{1};
% main_folder = fullfile(simulated_data_folder, file_name, facil);
% first_plot = plot_performance_across_variable_change(delays, main_folder, 'delay time');
% hold on;
% 
% for i = 2:length(facils)
%     facil = facils(i); facil = facil{1};
%     main_folder = fullfile(simulated_data_folder, file_name, facil);
%     plot_performance_across_variable_change(delays, main_folder, 'delay time');
%     hold on; 
% end
% hold off
% % needs to be done a bit differently 

% Get the current figure and axes
fig = gcf;
ax = gca;

% Change the size of the figure (Width x Height in pixels)
fig.Position = [100, 100, 550, 450]; % Adjust these values as needed

% Change the left y-axis limits
yyaxis left;
ylim([0.4 1.05]);

% Change the right y-axis limits
yyaxis right;
ylim([-500 3000]);
ylabel('Iterations for learning'); % Right y-axis label

% Set the font size for x and y-axis labels
ax.XLabel.FontSize = 12;
ax.YLabel.FontSize = 12;
ax.YAxis(2).Label.FontSize = 12; % Right y-axis label font size

% Optionally, set the font size for the tick labels as well
ax.XAxis.FontSize = 12;
ax.YAxis(1).FontSize = 12; % Left y-axis tick labels
ax.YAxis(2).FontSize = 12; % Right y-axis tick labels
