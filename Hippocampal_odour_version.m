% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers

%% initial set up 
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_first = 50;
lenght_second = 20;
delay_time = 800;
start_time = 200;

% Get non-programmable paramaters
p = get_params_hipp(0.85);

% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);

% get cells to be activated by each odour in layer CA3
mems = get_odours_hipp(p, degree_overlap, pattern_order);

% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_first)];
input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
% M = get_memory_hipp(p);
M = simulate_dynapics_hipp(p, C, J, input, M, mems);

%% looking at mean firing rate of selective population during second odour 
% for different delay times 
for delay_time = 100:100:1500
    input.simulation = [start_time (start_time+length_first)];
    input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
    M = get_memory_hipp(p);
    M_new = simulate_dynapics_hipp(p, C, J, input, M, mems);
    
    ind1 = find((sum(C(mems{1}, :))>filter)); 
    ind2 = find((sum(C(mems{2}, :))>filter)); 
    coi = intersect(ind1,ind2);
    full = 1:p.out; full(coi) = [];
    
    n_trials = 10;
    
    overlapping_log = zeros(n_trials, 1);
    non_overlapping_log = zeros(n_trials, 1);
    
    for i = 1:n_trials

        M = simulate_dynapics_hipp(p, C, J, input, M, mems);
        output_mem2 = C(mems{2}, :);
        overlapping = M.spikelog(p.in + coi, input.reactivation(1):input.reactivation(2));
        non_overlapping = M.spikelog(p.in+full,  input.reactivation(1):input.reactivation(2));
    
        mean_overlapping = sum( overlapping, "all")/(size(overlapping, 1).*lenght_second.*0.001);
        mean_non_overlapping = sum(non_overlapping, "all")/(size(non_overlapping, 1).*lenght_second.*0.001);
    
        overlapping_log(i) = mean_overlapping; non_overlapping_log(i) = mean_non_overlapping;
        disp(i)

    end
    
    over = mean(overlapping_log); 
%% for classifier section, WIP
% % generate training data
n_trials = 50;
spikes_x_trials = zeros(n_trials,p.out+1);
% simulate 50 trials presenting odour A and then B, labelled no reward
for i = 1:n_trials/2
    M = simulate_dynapics_hipp(p, C, J, input, M, mems);
    spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
    spikes_out = sum(spikes, 2);
    spikes_x_trials(i, 1:p.out) = spikes_out';
    % labelling as no reward
    spikes_x_trials(i, p.out+1) = 0;
    disp(i)
end

% simulate another 50 trials but this time with C first, labelled reward
pattern_order = 'CB';
mems = get_odours_hipp(p, degree_overlap, pattern_order);
M = get_memory_hipp(p);

for i =  (n_trials/2)+1:n_trials
    M = simulate_dynapics_hipp(p, C, J, input, M, mems);
    spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
    spikes_out = sum(spikes, 2);
    spikes_x_trials(i, 1:p.out) = spikes_out';
    % labelling as reward
    spikes_x_trials(i, p.out+1) = 1;
    disp(i)
end

% shuffle trials randomly
shuffled_spikes_x_trials = spikes_x_trials(randperm(size(spikes_x_trials,1)),:);

% train on progressively more data to see if it improves
class_loss_log2 = [];
for i = 2:n_trials
    SVMModel = fitcsvm(shuffled_spikes_x_trials(1:i, 1:p.out), shuffled_spikes_x_trials(1:i, end));
    CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(CVSVMModel);
    disp(classLoss)
    class_loss_log2 = [class_loss_log2; classLoss];
end
% plot class log for increasing number of trials provided to classifier
hold on;
plot(2:n_trials, class_loss_log2)
xlabel('number of observations');
ylabel('class loss');

% nothing
% X = shuffled_spikes_x_trials(1:i, 1:p.out);
% Y = shuffled_spikes_x_trials(1:i, end);
% CVSVMModel = fitcsvm(X,Y,'Holdout',0.15);
% CompactSVMModel = CVSVMModel.Trained{1}; % Extract the trained, compact classifier
% testInds = test(CVSVMModel.Partition);   % Extract the test indices
% XTest = X(testInds,:);
% YTest = Y(testInds,:);
% L = loss(CompactSVMModel,XTest,YTest);

%% plotting 

% plotting vm in both patterns
fs = 10;
ns = 7; 
n = 0;
% % V mean for input and output layer
% subplot(ns, 1, 1)
% x = M.V_log_in;
% plot(1:p.SimLength, x)
% hold on;
% x = M.V_log_out;
% plot(1:p.SimLength, x)
% hold on;
% x_points = [input.simulation(1), input.simulation(1), input.simulation(2), input.simulation(2)];
% x2_points =[input.reactivation(1), input.reactivation(1), input.reactivation(2), input.reactivation(2)];
% y_points = [0, 1, 1, 0];
% color = [0, 0, 1];
% hold on;
% a = fill(x_points, y_points, color,'LineStyle','none');
% a.FaceAlpha = 0.1; 
% hold on;
% a = fill(x2_points, y_points, color,'LineStyle','none');
% a.FaceAlpha = 0.1; 
% hold on;
% ylabel('Vm in memory','FontSize',fs)
% legend('V out', 'V in')

% plotting synaptic parameters u an x for synapses potentiated by first
% odour which are reactivated during second odour presentation
subplot(ns,1,n+1)
av_u_memory = M.U_mem1_log;
av_x_memory = M.X_mem1_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(pattern_order(1),'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('overlapping')
%u and x for non overlapping 2nd odor cells
subplot(ns,1,n+2)
av_u_memory = M.U_mem2_log;
av_x_memory = M.X_mem2_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(pattern_order(2),'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('non-overlapping')

%plot  current logs
subplot(ns, 1, n+3)
av_e = mean(M.Iext_log(1:p.in, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('ext-->CA3','FontSize',fs)

subplot(ns, 1, n+4)
av_i = mean(M.Irec_log(1:p.out, :), 1); 
plot(1:p.SimLength, av_i)
ylabel('CA3-->CA1','FontSize',fs)

subplot(ns, 1, n+5)
av_e = mean(M.Iext_log(p.in+1:p.full, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('ext-->CA1','FontSize',fs)

%plotting spike raster for CA3 cells odour 1 and odour 2
color_ops = { 'b', 'r'};
ops = {0, 2};
subplot(ns, 1, n+6)
for m = 1:2
    tVec = 1:p.SimLength;
    spikeMat = M.spikelog(mems{m}, :);
    hold all; 
    for trialCount = 1:size(spikeMat,1)
        if sum(spikeMat(trialCount, :)) == 0
            continue
        else
            spikePos = tVec(find(spikeMat(trialCount, :)));
            for spikeCount = 1:length(spikePos)
                plot([spikePos(spikeCount)+ops{m} spikePos(spikeCount)+ops{m}], ...
                [trialCount-0.4 trialCount+0.4], color_ops{m});
                
            end
        end
    end
    hold on;
end

spikeMat = M.spikelog; spikeMat(mems{1},:) = []; spikeMat(mems{2},:) = [];
spikeMat = spikeMat(randperm(size(spikeMat,1)),:);
spikeMat = spikeMat(1:p.full/10, :);
for trialCount = 1:size(spikeMat,1)
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'Color', 'm');
        end
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('CA3','FontSize',fs)

L1 = plot(nan, nan, 'color', 'b');
L2 = plot(nan, nan, 'color', 'r');
L3 = plot(nan, nan, 'color', 'm');
legend([L1, L2, L3], {'odour 1' 'odour 2', '10% over non-odour cells'}, 'Location', 'southeast')

% plotting spike raster for CA1 ouput cell firing
ind1 = find((sum(C(mems{2}, :))>p.in.*0.01)); 
ind2 = find((sum(C(mems{1}, :))>p.in.*0.01)); 
coi = intersect(ind1,ind2);
full = 1:p.out; full(coi) = [];

% overlapping with more than 2 inputs
subplot(ns, 1, n+7)
spikeMat = M.spikelog(p.in+coi, :);
for trialCount = 1:size(spikeMat,1)
    hold all;
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount)+1 spikePos(spikeCount)+1], ...
            [trialCount-1 trialCount+1], 'Color', 'm');
        end
    end
end
hold on;
% others (non overlapping and those receiving <2 inputs)
spikeMat = M.spikelog(p.in+full, :);
for trialCount = 1:size(spikeMat,1)
    hold all;
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-1 trialCount+1], 'Color', 'b');
        end
    end
end
% ylim([0 p.in]);
xlim([0 p.SimLength])
ylabel('CA1','FontSize',fs)


L1 = plot(nan, nan, 'color', 'm');
L2 = plot(nan, nan, 'color', 'b');
legend([L1, L2], {'overlapping ' 'non-overlapping'}, 'Location', 'southeast')

%plotting mean spiking in overlapping cells vs non overlapping cells
%during second odour, filtering for cells which recived increasing number of inputs
for filter = 1:5
    
    ind1 = find((sum(C(mems{1}, :))>filter)); 
    ind2 = find((sum(C(mems{2}, :))>filter)); 
    coi = intersect(ind1,ind2);
    full = 1:p.out; full(coi) = [];
    
    n_trials = 10;
    
    overlapping_log = zeros(n_trials, 1);
    non_overlapping_log = zeros(n_trials, 1);
    
    for i = 1:n_trials

        M = simulate_dynapics_hipp(p, C, J, input, M, mems);
    
        overlapping = M.spikelog(p.in + coi, input.reactivation(1):input.reactivation(2));
        non_overlapping = M.spikelog(p.in+full,  input.reactivation(1):input.reactivation(2));
    
        mean_overlapping = sum( overlapping, "all")/(size(overlapping, 1).*lenght_second.*0.001);
        mean_non_overlapping = sum(non_overlapping, "all")/(size(non_overlapping, 1).*lenght_second.*0.001);
    
        overlapping_log(i) = mean_overlapping; non_overlapping_log(i) = mean_non_overlapping;
        disp(i)

    end
    
    over = mean(overlapping_log); 
    
    % plot non overlapping cells
    if filter == 1 
        non_over = mean(non_overlapping_log);
        scatter((1*ones(n_trials)), non_overlapping_log,'MarkerEdgeColor','r', 'MarkerFaceColor','r')
        hold on;
        plot([1-0.1 1.+0.1], [non_over non_over], 'Color', 'r', 'LineWidth',2);
        hold on;
    end

    scatter(filter+(1*ones(n_trials)),overlapping_log, 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
    hold on;
    plot([filter+1-0.1 filter+1+0.1], [over over],'Color', 'b', 'LineWidth',2);

end
xticks([1 2 3 4 5 6])
xticklabels({'non-overlap', '>1 input', '>2 input', '>3 input', '>4 input', '>5 input'})

xlabel('CA1 cells receiving non overlapping input and overlapping input from CA3')
ylabel('mean firing rate during second odour presentation (Hz)')

L1 = plot(nan, nan, 'color', 'r');
L2 = plot(nan, nan, 'color', 'b');
legend([L1, L2 ], {'non-overlapping' 'overlapping'}, 'Location', 'southeast')