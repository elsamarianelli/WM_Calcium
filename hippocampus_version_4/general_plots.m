%% general plots 
% plots firing rate during odour 1, delay period, and beginning of odour 2,
% for populations of interest
%% Set parameters for the simulation
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AB';
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 250;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 500;            % Delay between odour presentations (ms)
p.length_second         = 250;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.875;           % Constant by which to scale random currents (to modulate baseline activity levels)
p.SimLength             = 1500;
p                       = get_params_hipp(p);

%  Assign CA3 and CA1 cells to each odour representation
ca3_ensembles           = get_odours_hipp(1:p.in, p.f, p.degree_overlap_CA3);
ca1_ensembles           = get_odours_hipp(p.in+(1:p.out), p.f_o, p.degree_overlap_CA1);

% take first and second odour to be presented 
stim                    = cell(2,1); 
first                   = double(upper(p.pattern_order(1))) - 64; 
stim{1}                 = ca3_ensembles{first}; %clear first
second                  = double(upper(p.pattern_order(2))) - 64; 
stim{2}                 = ca3_ensembles{second}; %clear second

%  Generate connectivity and synaptic efficacy matrix
[C, J]                  = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles);
% see_connectivity        = visualise_connectivity(C, ca3_ensembles, ca1_ensembles);

%  Specify times that each odour is presented, assign memory for the output
input.simulation        = [p.start_time p.start_time+p.length_first];
input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
M                       = get_memory_hipp(p);
mems = stim;
%% 
%%  Simulate hippocampal dynamics 
M                       = simulate_dynamics_hipp(p, C, J, input, M, stim);

%  Plot output for a single trial 
output_plot             = get_output_plot(M, p, ca3_ensembles, C);
%% Simualte data: firing rate of CA1 cell populations during 1st and 2nd odours, and delay period
time1 = input.simulation(1);
time2 = input.simulation(2);
time3 = input.reactivation(1);
time4 = time3+100; %input.reactivation(2);

% get CA1 cells which are able to receive input from both A and B
CA3_overlap = intersect(ca3_ensembles{first}, ca3_ensembles{second});
ind = find((sum(C(CA3_overlap, :))>3)); 
CA1_overlap_first = [(intersect(ca1_ensembles{first}, ind+p.in))];
CA1_overlap_second =[(intersect(ca1_ensembles{second}, ind+p.in))];

not_ind = ~ismember(ca1_ensembles{second}, ind+p.in);
CA1_B = ca1_ensembles{second};
CA1_B = CA1_B(not_ind); 

not_ind = ~ismember(ca1_ensembles{first}, ind+p.in);
CA1_A = ca1_ensembles{first};
CA1_A = CA1_A(not_ind);

% Memory
number_of_trials = 10;

firing_rate_data = cell(2, 9);

firing_rate_data{1, 1} = 'A unique';
firing_rate_data{1, 2} = 'A overlap B';
firing_rate_data{1, 3} = 'B unique';

firing_rate_data{1, 4} = 'A unique';
firing_rate_data{1, 5} = 'A overlap B';
firing_rate_data{1, 6} = 'B unique';

firing_rate_data{1, 7} = 'A unique';
firing_rate_data{1, 8} = 'A overlap B';
firing_rate_data{1, 9} = 'B unique';

% for looking at firing rate of CA3 or CA1 populaitions
         % CA1
A       = CA1_A;
overlap_first = CA1_overlap_first;
overlap_second = CA1_overlap_second;
B       = CA1_B;
%         % CA3
% A       = ca3_ensembles{first};
% overlap = CA3_overlap;
% B       = ca3_ensembles{second};

for i = 1:number_of_trials
        
    % Simulate a trial
    M = get_memory_hipp(p);
    M = simulate_dynamics_hipp(p, C, J, input, M, stim);

    % [1] get firing rates during FIRST odour presentation

    % A unique
    Hz = get_mean_rate(M, time1, time2, A);
    firing_rate_data{2, 1} = [firing_rate_data{2, 1}, Hz]; clear Hz;
    % A overlap B 
    Hz = get_mean_rate(M, time1, time2, overlap_first);
    firing_rate_data{2, 2} = [firing_rate_data{2, 2}, Hz]; clear Hz;
    % B unique
    Hz = get_mean_rate(M, time1, time2, B);
    firing_rate_data{2, 3} = [firing_rate_data{2, 3}, Hz]; clear Hz;

    % [2] get firing rates during DELAY period

    % A unique
    Hz = get_mean_rate(M, time2, time3, A);
    firing_rate_data{2, 4} = [firing_rate_data{2, 4}, Hz]; clear Hz;
    % A overlap B 
    Hz = get_mean_rate(M, time2, time3, overlap_first);
    firing_rate_data{2, 5} = [firing_rate_data{2, 5}, Hz]; clear Hz;
    % B unique
    Hz = get_mean_rate(M, time2, time3, B);
    firing_rate_data{2, 6} = [firing_rate_data{2, 6}, Hz]; clear Hz;

    % [3] get firing rates during SECOND odour presentation
    % A unique
    Hz = get_mean_rate(M, time3, time4, A);
    firing_rate_data{2, 7} = [firing_rate_data{2, 7}, Hz]; clear Hz;
    % A overlap B 
    Hz = get_mean_rate(M, time3, time4, overlap_second);
    firing_rate_data{2, 8} = [firing_rate_data{2, 8}, Hz]; clear Hz;
    % B unique
    Hz = get_mean_rate(M, time3, time4, B);
    firing_rate_data{2, 9} = [firing_rate_data{2, 9}, Hz]; clear Hz;

    disp(i)

end


%% swarm chart with SD error bars
% Extract variable names and data
variable_names = firing_rate_data(1,:);
data = firing_rate_data(2,:);

% A / overlap / B colours
colors = {'b', 'm', 'r'}; 

% Group names 
group_names = {'first odour', 'delay', 'second odour'};

% Setup the figure
figure;

% Total number of subplots needed
totalSubplots = 3;

% Loop through each subplot
for subplotIdx = 1:totalSubplots
    subplot(1, totalSubplots, subplotIdx); % Arrange in 1 row, with totalSubplots columns
    
    % Loop through each category within a subplot
    for subpop = 1:3

        % Calculate index for data 
        dataIndex = (subplotIdx - 1) * 3 + subpop;
  
        if dataIndex > length(data)
            continue; % Skip if dataIndex exceeds the data length
        end
        
        % Generate jittering for axes to simulate the swarm effect
        xJitterValues = 0.1 * rand(1, numel(data{dataIndex})) - 0.05 + subpop;
        yJitterValues = data{dataIndex} + (0.05 * rand(1, numel(data{dataIndex})) - 0.025); % Adjust the magnitude of jittering as needed
        
        scatter(xJitterValues, yJitterValues, colors{subpop});
        hold on;
            
        % Calculate mean and SD (Standard Deviation)
        mean_value = mean(data{dataIndex});
        sd_value = std(data{dataIndex}); % SD calculation
        
        % mean bar
        plot([subpop-0.1, subpop+0.1], [mean_value, mean_value], '-k', 'LineWidth', 2);
        
        % error bar for SD
        errorbar(subpop, mean_value, sd_value, 'k', 'LineWidth', 2, 'CapSize', 10);
    end
    
    xticks(1:3);
    xticklabels({'A unique', 'A overlap B', 'B unique'});
    ylabel('Hz');
    title([group_names{subplotIdx}]);
    % ylim([0 20])
end