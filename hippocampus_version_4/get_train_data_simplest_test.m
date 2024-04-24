function[shuffled_spikes_x_trials,sequenceID] = get_train_data_simplest_test(C, J, n_trials, p, ca3_ensembles, time_1, time_2)

%%  Simulate hippocampal dynamics  over many trials with simplest version of task (only 2 odours)
odour_sequences     = {'CB'; 'BC'};
reward_outcome      = [1 0];
spikes_x_trials     = zeros(n_trials.*6,p.out+1);
sequenceID          = cell(n_trials.*6,1);

%% Generate training data for each odour sequence
%  Loop through each input pattern
for pattern         = 1 : length(odour_sequences) 
    
    % Identify input neurons for each odour in this sequence
    mems_trial      = cell(2,1); 
    p.pattern_order = odour_sequences{pattern}; 
    disp(['Starting ' p.pattern_order ' trials...']);
    first           = double(upper(p.pattern_order(1))) - 64; 
    mems_trial{1}   = ca3_ensembles{first}; clear first
    second          = double(upper(p.pattern_order(2))) - 64; 
    mems_trial{2}   = ca3_ensembles{second}; clear second
    stim = mems_trial;
    % simulate n_trials of these odour presentations
    for i           = 1 : n_trials
        
        % update the user
        disp(['Trial ' int2str(i) ' with odour pair ' p.pattern_order]);
        
        % assign memory and run the dynamics
        n           = (pattern-1)*n_trials + i;
        M           = get_memory_hipp(p);
        M                       = simulate_dynamics_hipp(p, C, J, input, M, stim);
        
        % log spiking activity during time of interest
        spikes      = M.spikelog(p.in+1:p.full,time_1:time_2);
        spikes      = sum(spikes, 2);
        spikes_x_trials(n, 1:p.out) = spikes'; clear spikes M
        sequenceID{n}               = p.pattern_order;
        
        % label as a rewarded or un-rewarded trial
        spikes_x_trials(n, p.out+1) = reward_outcome(pattern); clear n        
        
    end
    clear i mems_trial
    
end
clear C ca3_ensembles input J n_trials odour_sequences p pattern reward_outcome

% shuffle trials randomly
randOrd                     = randperm(size(spikes_x_trials,1));
shuffled_spikes_x_trials    = spikes_x_trials(randOrd,:); clear spikes_x_trials
sequenceID                  = sequenceID(randOrd); clear randOrd

end
