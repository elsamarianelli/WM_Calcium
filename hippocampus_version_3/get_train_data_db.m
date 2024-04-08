function[shuffled_spikes_x_trials,sequenceID] = get_train_data(C, J, n_trials, p, ca3_ensembles)
% at present this code is being used to get input to perceptron, the output
% is a mean firing rate (Hz) for each CA1 cell during the second
% odour presentation. Code to run an SVM instead of using a perceptron, as
% well as code to run multiple delay times is at the bottom, commented out.


%% Assign some parameters and memory for the output
odour_sequences     = {'CB'; 'BA'; 'AC'; 'AB'; 'CA'; 'BC'};
reward_outcome      = [1 1 1 0 0 0];
input.simulation    = [p.start_time p.start_time+p.length_first];
input.reactivation  = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];
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
    
    % simulate n_trials of these odour presentations
    for i           = 1 : n_trials
        
        % update the user
        disp(['Trial ' int2str(i) ' with odour pair ' p.pattern_order]);
        
        % assign memory and run the dynamics
        n           = (pattern-1)*n_trials + i;
        M           = get_memory_hipp(p);
        M           = simulate_dynamics_hipp(p, C, J, input, M, mems_trial);
        
        % log spiking activity during the second odour presentation
        spikes      = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
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