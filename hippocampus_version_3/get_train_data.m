function  shuffled_spikes_x_trials = get_train_data(C, J, input, n_trials, p, mems_all)
% at present this code is being used to get input to perceptron, the output
% is a mean firing rate (Hz) for each CA1 cell during the second
% odour presentation. Code to run an SVM instead of using a perceptron, as
% well as code to run multiple delay times is at the bottom, commented out.

%% store n spikes during second odour 
% for all cells in CA1 with trial labelled as reward/no reward
spikes_x_trials = zeros(n_trials,p.out+1);

% what times is spiking being looked at from /to
start_time = input.reactivation(1); end_time = input.reactivation(2) + 100;
%% generate training data reward conditions
reward_patterns = {'CB'; 'BA'; 'AC'};

for pattern = 1:length(reward_patterns) 

    % get CA3 cell ensembles which odour 1 and 2 will activate
    mems_trial = cell(2,1); 
    pattern_order = reward_patterns{pattern}; disp(pattern_order)
    first = double(upper(pattern_order(1))) - 64; 
    mems_trial{1} = mems_all{first};
    second = double(upper(pattern_order(2))) - 64; 
    mems_trial{2} = mems_all{second};

    % simulate trials presenting odour A and then B, labelled no reward
    for i = 1:n_trials/6
        n = ((pattern.*n_trials/6)-n_trials/6) + i;
        M = get_memory_hipp(p);
        %  Simulate hippocampal dynamics 
        M = simulate_dynamics_hipp(p, C, J, input, M, mems_trial);
        spikes = M.spikelog(p.in+1:p.full, start_time:end_time);
        spikes_out = sum(spikes, 2);
        spikes_x_trials(n, 1:p.out) = spikes_out';
        % labelling as reward
        spikes_x_trials(n, p.out+1) = 1;
        disp(i); disp(p.pattern_order);

    end
end

%% generate training data no reward conditions
no_reward_patterns = {'AB'; 'CA'; 'BC'};

for pattern = 1:length(no_reward_patterns) 

    % get CA3 cell ensembles which odour 1 and 2 will activate
    mems_trial = cell(2,1); 
    pattern_order = no_reward_patterns{pattern}; disp(pattern_order)
    first = double(upper(pattern_order(1))) - 64; 
    mems_trial{1} = mems_all{first};
    second = double(upper(pattern_order(2))) - 64; 
    mems_trial{2} = mems_all{second};

    % simulate trials presenting odour A and then B, labelled no reward
    for i = 1:n_trials/6
        n = ((pattern.*n_trials/6)+n_trials/3) + i;
        M = get_memory_hipp(p);
        %  Simulate hippocampal dynamics 
        M = simulate_dynamics_hipp(p, C, J, input, M, mems_trial);
        spikes = M.spikelog(p.in+1:p.full, start_time:end_time);
        spikes_out = sum(spikes, 2);
        spikes_x_trials(n, 1:p.out) = spikes_out';
        % labelling as no-reward
        spikes_x_trials(n, p.out+1) = 0;
        disp(i); disp(pattern_order);

    end
end


%% shuffle trials randomly 
% (so that perceptron online training receives varied input instead of blocks of odour orders)
shuffled_spikes_x_trials = spikes_x_trials(randperm(size(spikes_x_trials,1)),:);
end