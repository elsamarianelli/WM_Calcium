%% general plots 

%% Set parameters for the simulation
p.degree_overlap_CA3    = 0.2;            % Overlap between neural representations of each odour
p.degree_overlap_CA1    = 0.0;
p.pattern_order         = 'AC';           % Order in which the odours should be presented
p.start_time            = 200;            % Time at which the first odour is presented (ms)
p.length_first          = 40;             % Length of time for which the first odour is presented (ms)
p.delay_time            = 600;            % Delay between odour presentations (ms)
p.length_second         = 40;             % Length of time for which the second odour is presented (ms)
p.scaleF                = 0.85;           % Constant by which to scale random currents (to modulate baseline activity levels)
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

%% 
%%  Simulate hippocampal dynamics 
M                       = simulate_dynamics_hipp(p, C, J, input, M, stim);

%  Plot output for a single trial 
output_plot             = get_output_plot(M,p.pattern_order, p, stim, C);

%% plot 1     M = simulate_dynamics_hipp(p, C, J, input, M, mems);
% [mean_firing_second_odour] = get_mean_firing_second_odour(p, C, J, input, M, mems, p.length_second);

% %% get rates for different populations
% number_of_trials = 30;
% rates_sorted_all = zeros(18, number_of_trials);
% sub_groups = [1, first+1, second+1];
% for i = 1:number_of_trials
%     M             = get_memory_hipp(p); 
%     M             = simulate_dynamics_hipp(p, C, J, input, M, stim);
%     rates_sorted  = sortRates(p,M,input,ca3_ensembles,ca1_ensembles);
% 
%     for n = 1:3
%         sub_group = sub_groups(n);
%         rates_sorted_all(n, i) = rates_sorted.ca1stim1Mn(sub_group);
%         rates_sorted_all(n+3, i) = rates_sorted.ca1delayMn(sub_group);
%         rates_sorted_all(n+6, i) = rates_sorted.ca1stim2Mn(sub_group);
% 
%         rates_sorted_all(n+9, i) = rates_sorted.ca3stim1Mn(sub_group);
%         rates_sorted_all(n+12, i) = rates_sorted.ca3delayMn(sub_group);
%         rates_sorted_all(n+15, i) = rates_sorted.ca3stim2Mn(sub_group);
% 
%         disp(sub_group)
%     end
% 
%     field_name = ['run', num2str(i)];
%     disp(field_name)
% end
% 
% %% plot swarm chart                            
% variable_names = {'ca1stim1Mn background'
%                  ['ca1stim1Mn ',p.pattern_order(1)]
%                  ['ca1stim1Mn ',p.pattern_order(2)]
%                  'ca1delayMn background'
%                  ['ca1delayMn ',p.pattern_order(1)]
%                  ['ca1delayMn ',p.pattern_order(2)]
%                  'ca1stim2Mn background'
%                  ['ca1stim2Mn ',p.pattern_order(1)]
%                  ['ca1stim2Mn ',p.pattern_order(2)]
%                  'ca3stim1Mn background'
%                  ['ca3stim1Mn ',p.pattern_order(1)]
%                  ['ca3stim1Mn ',p.pattern_order(2)]
%                  'ca3delayMn background'
%                  ['ca3delayMn ',p.pattern_order(1)]
%                  ['ca3delayMn ',p.pattern_order(2)]
%                  'ca3stim2Mn background'
%                  ['ca3stim2Mn ',p.pattern_order(1)]
%                  ['ca3stim2Mn ',p.pattern_order(2)]};
% 
% data = rates_sorted_all;
% 
% figure;
% for i = 1:numel(variable_names)
% 
%     % Swarm plot
%     swarmchart(i.*ones(1, number_of_trials), data(i, :));
%     hold on;
% 
%     % Add mean bar
%     mean_value = mean(data(i,:));
%     plot([i-0.1 i+0.1], [mean_value mean_value], '-k', 'LineWidth', 2);
% 
%     hold on;
% end  
% 
% % Set x-axis tick and label
% xticks(1:numel(variable_names)); 
% xticklabels(strrep(variable_names, '_', ' '));
% ylabel('Hz');
%% plot 2
number_of_trials = 50;

firing_rate_data = cell(2, 5);

firing_rate_data{1, 1} = 'background_CA3_Hz';
firing_rate_data{1, 2} = 'background_CA1_Hz';
firing_rate_data{1, 3} = 'CA3_stimulated_Hz';
% firing_rate_data{1, 3} = 'CA3_overlap_Hz';
firing_rate_data{1, 4} = 'CA1_stimulated_Hz';
firing_rate_data{1, 5} = 'CA1_overlap_Hz';


% get CA1 cells which are able to receive input from both A and B
overlap_CA3 = intersect(ca3_ensembles{first}, ca3_ensembles{second});
ind = find((sum(C(overlap_CA3, :))>0)); 
CA1_overlap = intersect(ca1_ensembles{second}, ind+p.in);
not_ind = ~ismember(ca1_ensembles{second}, ind+p.in);
CA1_non_overlap = ca1_ensembles{second};
CA1_non_overlap = CA1_non_overlap(not_ind);
% 
CA1_overlap = CA1_overlap(1:10);
CA1_non_overlap = CA1_non_overlap(1:10);

% get subpopulations to plot
% filter = 8;
% overlap_CA3 = intersect(ca3_ensembles{first}, ca3_ensembles{second});
% % overlap_CA1 = intersect(ca1_ensembles{first}, ca1_ensembles{second});
% ind1 = find((sum(C(mems{1}, :))>filter)); 
% ind2 = find((sum(C(mems{2}, :))>filter)); 
% coi = intersect(ind1,ind2);
% full = find((sum(C(mems{2}, :))>0)); 
% % Find elements in list1 that are not in list2
% logical_index = ~ismember(full, coi);
% % Use logical indexing to select elements from list1
% CA1_non_overlap = full(logical_index);
% overlap_CA1 = coi;
% % CA1_non_overlap = ca1_ensembles{second};
% % ind = ~ismember(CA1_non_overlap, overlap_CA1);
% CA1_non_overlap = CA1_non_overlap(ind);

time1 = input.reactivation(1);
time2 = input.reactivation(2);
% time1 = input.reactivation(2);
% time2 = input.reactivation(2)+200;

for i = 1:number_of_trials
        
    % simulate a trial
    M                   = get_memory_hipp(p);
    M = simulate_dynamics_hipp(p, C, J, input, M, stim);
    
    % get mean background firing rate CA3 (in)
    background_CA3_Hz = get_mean_rate(M, 1, p.SimLength, 1:p.in); 
    firing_rate_data{2, 1} = [firing_rate_data{2, 1}, background_CA3_Hz];

    % get mean background firing rate CA1 (out)
    background_CA1_Hz = get_mean_rate(M, 1, p.SimLength, p.in+1:p.full);
    firing_rate_data{2, 2} = [firing_rate_data{2, 2}, background_CA1_Hz];

    % get mean firing rate of presynaptic simulated cells during odour 2
    CA3_simulated_Hz = get_mean_rate(M, time1, time2, ca3_ensembles{second});
    firing_rate_data{2, 3} = [firing_rate_data{2, 3}, CA3_simulated_Hz];
    % 
    % % get mean firing rate of overlapping CA3 cells during odour 2
    % CA3_overlap_Hz = get_mean_rate(M, input.reactivation(1), input.reactivation(2), overlap_CA3);
    % firing_rate_data{2, 3} = [firing_rate_data{2, 3}, CA3_overlap_Hz];

    % get mean firing rate of postsynaptic simulated cells during odour 2
    CA1_simulated_Hz = get_mean_rate(M, time1, time2, CA1_non_overlap);
    firing_rate_data{2, 4} = [firing_rate_data{2, 4}, CA1_simulated_Hz];

    % get mean firing rate of overlaping CA1 cells during odour 2
    CA1_overlap_Hz = get_mean_rate(M, time1, time2, CA1_overlap);
    firing_rate_data{2, 5} = [firing_rate_data{2, 5}, CA1_overlap_Hz];

    disp(i)
end
        
%% plot swarm chart                            
variable_names = firing_rate_data(1,:);
data = firing_rate_data(2,:);

figure;
for i = 1:numel(variable_names)
    
    % Swarm plot
    swarmchart(i.*ones(1, number_of_trials), data{i});
    hold on;
    
    % Add mean bar
    mean_value = mean(data{i});
    plot([i-0.1 i+0.1], [mean_value mean_value], '-k', 'LineWidth', 2);
       
    hold on;
end  

% Set x-axis tick and label
xticks(1:numel(variable_names)); 
xticklabels(strrep(variable_names, '_', ' '));
ylabel('Hz');

% Extract variable names and data
variable_names = firing_rate_data(1,:);
data = firing_rate_data(2,:);
