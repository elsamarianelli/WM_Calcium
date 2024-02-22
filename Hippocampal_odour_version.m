% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers

%% initial set up 
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_stimulation = 50;
delay_time = 600;
start_time = 200;

% Get non-programmable paramaters
p = get_params_hipp(0.85);

% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);

% get cells to be activated by each odour in layer CA3
mems = get_odours_hipp(p, degree_overlap, pattern_order);

% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_stimulation)];
input.reactivation = [(start_time+length_stimulation+delay_time) (start_time+length_stimulation+length_stimulation+delay_time)];

% Generate memory to run simulation
M = get_memory_hipp(p);

% calculate degree that ouput layer cells receiving input overlap between 2
% memories 
overlap1 = mems{2};
ind1 = find((sum(C(overlap1, :)))); 
overlap2 = mems{1};
ind2 = find((sum(C(overlap2, :)))); 

degree_overlap_ouput = size(intersect(ind1,ind2), 2)/p.out;

M = simulate_dynapics_hipp(p, C, J, input, M, mems);
% %  %% generate training data
n_trials = 50;
spikes_x_trials = zeros(n_trials,p.out+1);
%simulate 50 trials presenting odour A and then B
for i = 1:n_trials/2
    M = simulate_dynapics_hipp(p, C, J, input, M, mems);
    spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
    spikes_out = sum(spikes, 2);
    spikes_x_trials(i, 1:p.out) = spikes_out';
    % labelled as no reward
    spikes_x_trials(i, p.out+1) = 0;
    disp(i)
end

% simulate another 50 trials but this time with C first, which is labelled
% as reward
pattern_order = 'CB';
mems = get_odours_hipp(p, degree_overlap, pattern_order);
M = get_memory_hipp(p);

for i =  (n_trials/2)+1:n_trials
     M = simulate_dynapics_hipp(p, C, J, input, M, mems);

    spikes = M.spikelog(p.in+1:p.full, input.reactivation(1):input.reactivation(2));
    spikes_out = sum(spikes, 2);
    spikes_x_trials(i, 1:p.out) = spikes_out';
    % labelled as reward
    spikes_x_trials(i, p.out+1) = 1;
    disp(i)
end
% shuffle trials randomly
shuffled_spikes_x_trials = spikes_x_trials(randperm(size(spikes_x_trials,1)),:);

% train on progressively more data to see if it improves
class_loss_log = [];
for i = 2:n_trials

    SVMModel = fitcsvm(shuffled_spikes_x_trials(1:i, 1:p.out), shuffled_spikes_x_trials(1:i, end));
    CVSVMModel = crossval(SVMModel);
    classLoss = kfoldLoss(CVSVMModel);
    disp(classLoss)
    class_loss_log = [class_loss_log; classLoss];

end
plot(2:n_trials, class_loss_log)
xlabel('number of observations');
ylabel('class loss');
% % %
% hi = shuffled_spikes_x_trials(1:n_trials, end);
% shuffled_spikes_x_trials(1:n_trials, end) = hi(randperm(size(hi,1)),:);

X = shuffled_spikes_x_trials(1:i, 1:p.out);
Y = shuffled_spikes_x_trials(1:i, end);
CVSVMModel = fitcsvm(X,Y,'Holdout',0.15);
CompactSVMModel = CVSVMModel.Trained{1}; % Extract the trained, compact classifier
testInds = test(CVSVMModel.Partition);   % Extract the test indices
XTest = X(testInds,:);
YTest = Y(testInds,:);
L = loss(CompactSVMModel,XTest,YTest);

%% plotting 
%% Plot the output, if required
% plotting vm in both patterns

fs = 10;
ns = 8; 

% V mean for input and output layer
subplot(ns, 1, 1)
x = M.V_log_in;
plot(1:p.SimLength, x)
hold on;
x = M.V_log_out;
plot(1:p.SimLength, x)
hold on;
x_points = [input.simulation(1), input.simulation(1), input.simulation(2), input.simulation(2)];
x2_points =[input.reactivation(1), input.reactivation(1), input.reactivation(2), input.reactivation(2)];
y_points = [0, 1, 1, 0];
color = [0, 0, 1];
hold on;
a = fill(x_points, y_points, color,'LineStyle','none');
a.FaceAlpha = 0.1; 
hold on;
a = fill(x2_points, y_points, color,'LineStyle','none');
a.FaceAlpha = 0.1; 
hold on;
ylabel('Vm in memory','FontSize',fs)
legend('V out', 'V in')

% plotting synaptic parameters u an x
subplot(ns,1,2)
av_u_memory = M.U_mem1_log;
av_x_memory = M.X_mem1_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(first_input,'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('odour 1')
%u and x for second memory 
subplot(ns,1,3)
av_u_memory = M.U_mem2_log;
av_x_memory = M.X_mem2_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(second_input,'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('odour 2')
%plot  current logs
subplot(ns, 1, 4)
av_e = mean(M.Iext_log(1:p.in, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('ext-->CA3','FontSize',fs)

subplot(ns, 1, 5)
av_i = mean(M.Irec_log(1:p.out, :), 1); 
plot(1:p.SimLength, av_i)
ylabel('CA3-->CA1','FontSize',fs)

subplot(ns, 1, 6)
av_e = mean(M.Iext_log(p.in+1:p.full, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('ext-->CA1','FontSize',fs)

%plotting spike raster
color_ops = { 'b', 'r'};
ops = {0, 2};
subplot(ns, 1, 7)
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
ylabel('CA3 spiking','FontSize',fs)

% overlap = intersect(mems{1}, mems{2});
overlap1 = mems{2};
ind1 = find((sum(C(overlap1, :))>3)); 
overlap2 = mems{1};
ind2 = find((sum(C(overlap2, :))>3)); 
coi = intersect(ind1,ind2);
% coi = p.in+1:p.full;

subplot(ns, 1, 8)
spikeMat = M.spikelog(coi, :);
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
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('odours CA1 overlap','FontSize',fs)
% legend('10% of rest')