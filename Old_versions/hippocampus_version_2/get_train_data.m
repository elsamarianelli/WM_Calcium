function  shuffled_spikes_x_trials = get_train_data(C, J, input, n_trials, p, mems_all)
% at present this code is being used to get input to perceptron, the output
% is a mean firing rate (Hz) for each CA1 cell during the second
% odour presentation. Code to run an SVM instead of using a perceptron, as
% well as code to run multiple delay times is at the bottom, commented out.

%% looking at mean firing rate of selective population during second odour 
spikes_x_trials = zeros(n_trials,p.out+1);
reward_patterns = {'CB'; 'BA'; 'AC'};
% mems_all = get_odours_hipp(p, degree_overlap, "OFF");

% generate training data reward conditions
for pattern = 1:length(reward_patterns) 
    mems_trial = cell(2,1); 
    p.pattern_order = reward_patterns{pattern}; disp(p.pattern_order)
    first = double(upper(p.pattern_order(1))) - 64; 
    mems_trial{1} = mems_all{first};
    second = double(upper(p.pattern_order(2))) - 64; 
    mems_trial{2} = mems_all{second};
    % simulate 50 trials presenting odour A and then B, labelled no reward
    for i = 1:n_trials/6
        n = ((pattern.*n_trials/6)-n_trials/6) + i;
        M = get_memory_hipp(p);
        M = simulate_dynamics_hipp(p, C, J, input, M, mems_trial);
        spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
        spikes_out = sum(spikes, 2);
        spikes_x_trials(n, 1:p.out) = spikes_out';
        % labelling as reward
        spikes_x_trials(n, p.out+1) = 1;
        % disp(spikes_x_trials(n, p.out+1))
        disp(i); disp(p.pattern_order);
    end
end
no_reward_patterns = {'AB'; 'CA'; 'BC'};
% generate training data no reward conditions
for pattern = 1:length(no_reward_patterns) 
    mems_trial = cell(2,1); 
    pattern_order = no_reward_patterns{pattern}; disp(pattern_order)
    first = double(upper(pattern_order(1))) - 64; 
    mems_trial{1} = mems_all{first};
    second = double(upper(pattern_order(2))) - 64; 
    mems_trial{2} = mems_all{second};
    % simulate 50 trials presenting odour A and then B, labelled no reward
    for i = 1:n_trials/6
        n = ((pattern.*n_trials/6)+n_trials/3) + i;
        M = get_memory_hipp(p);
        M = simulate_dynamics_hipp(p, C, J, input, M, mems_trial);
        spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
        spikes_out = sum(spikes, 2);
        spikes_x_trials(n, 1:p.out) = spikes_out';
        % labelling as no-reward
        spikes_x_trials(n, p.out+1) = 0;
        % disp(spikes_x_trials(i, p.out+1))
        disp(i); disp(pattern_order);
    end
end


% shuffle trials randomly
shuffled_spikes_x_trials = spikes_x_trials(randperm(size(spikes_x_trials,1)),:);


% if run_SVM == "on"
%     % train on progressively more data to see if it improves
%     class_loss_log2 = [];
%     for i = 2:n_trials
%         SVMModel = fitcsvm(shuffled_spikes_x_trials(1:i, 1:p.out), shuffled_spikes_x_trials(1:i, end));
%         CVSVMModel = crossval(SVMModel);
%         classLoss = kfoldLoss(CVSVMModel);
%         disp(classLoss)
%         class_loss_log2 = [class_loss_log2; classLoss];
%     end
%     % plot class log for increasing number of trials provided to classifier
%     hold on;
%     SVM_plot = plot(2:n_trials, class_loss_log2);
%     set(0,'DefaultFigureVisible','off')
%     xlabel('number of observations');
%     ylabel('class loss');
% else 
%     SVM_plot = [];
% end
% % 
% % nothing
% % X = shuffled_spikes_x_trials(1:i, 1:p.out);
% % Y = shuffled_spikes_x_trials(1:i, end);
% % CVSVMModel = fitcsvm(X,Y,'Holdout',0.15);
% % CompactSVMModel = CVSVMModel.Trained{1}; % Extract the trained, compact classifier
% % testInds = test(CVSVMModel.Partition);   % Extract the test indices
% % XTest = X(testInds,:);
% % YTest = Y(testInds,:);
% % L = loss(CompactSVMModel,XTest,YTest);


% for different delay times 
% for delay_time = 100:100:1500
%     input.simulation = [start_time (start_time+length_first)];
%     input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
%     M = get_memory_hipp(p);
%     M_new = simulate_dynapics_hipp(p, C, J, input, M, mems);
% 
%     ind1 = find((sum(C(mems{1}, :))>filter)); 
%     ind2 = find((sum(C(mems{2}, :))>filter)); 
%     coi = intersect(ind1,ind2);
%     full = 1:p.out; full(coi) = [];
% 
%     n_trials = 10;
% 
%     overlapping_log = zeros(n_trials, 1);
%     non_overlapping_log = zeros(n_trials, 1);
% 
%     for i = 1:n_trials
% 
%         M = simulate_dynapics_hipp(p, C, J, input, M, mems);
%         output_mem2 = C(mems{2}, :);
%         overlapping = M.spikelog(p.in + coi, input.reactivation(1):input.reactivation(2));
%         non_overlapping = M.spikelog(p.in+full,  input.reactivation(1):input.reactivation(2));
% 
%         mean_overlapping = sum( overlapping, "all")/(size(overlapping, 1).*lenght_second.*0.001);
%         mean_non_overlapping = sum(non_overlapping, "all")/(size(non_overlapping, 1).*lenght_second.*0.001);
% 
%         overlapping_log(i) = mean_overlapping; non_overlapping_log(i) = mean_non_overlapping;
%         disp(i)
% 
%     end
% 
%     over = mean(overlapping_log); 
% %% for classifier section, WIP
end