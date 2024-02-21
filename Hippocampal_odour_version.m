% Applying a simplified version of the Mongillo et al. 2008 synaptic theory of working memory model 
% to odour presentation in CA3 to CA1 hippocampal layers

% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_stimulation = 50;
delay_time = 500;
start_time = 200;

% Get non-programmable paramaters
p = get_params_hipp(0.85);

% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);

% get cells to be activated by each odour in layer CA3
[mems, first_input, second_input] = get_odours_hipp(p, degree_overlap, pattern_order);

% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_stimulation)];
input.reactivation = [(start_time+length_stimulation+delay_time) (start_time+length_stimulation+length_stimulation+delay_time)];

% Generate memory to run simulation
M = get_memory_hipp(p);

% simulate dynamics across trials and record spike count in reacitvation
% period for each CA1 output neuron
% 
n_trials = 5;
spikes_x_trials = zeros(n_trials, p.out);
% 
i = 1;

% for i = 1:n_trials
    M = simulate_dynapics_hipp(p, C, J, input, M, mems);

    spikes = M.spikelog(p.in+1:p.full, :);
    spikes = spikes(:, input.reactivation(1):input.reactivation(2));
    spikes_out = sum(spikes, 2);
    % spikes_x_trials(i, 1:p.out) = spikes';
    spikes = M.spikelog(1:p.in, :);
    spikes = spikes(:, input.reactivation(1):input.reactivation(2));
    spikes_in = sum(spikes, 2);

    
    subplot(2, 1, 1)
    heatmap(spikes_in)
    subplot(2, 1, 2)
    heatmap(spikes_out)
% end


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

%plot  current logs
subplot(ns, 1, 4)
av_e = mean(M.Iext_log(1:p.in, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('external to input layer','FontSize',fs)

subplot(ns, 1, 5)
av_i = mean(M.Irec_log(1:p.out, :), 1); 
plot(1:p.SimLength, av_i)
ylabel('input to output layer','FontSize',fs)

subplot(ns, 1, 6)
av_e = mean(M.Iext_log(p.in+1:p.full, :), 1); 
plot(1:p.SimLength, av_e)
ylabel('external input to output layer','FontSize',fs)

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
legend('mem1', 'mem2')


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
ylabel('spike log','FontSize',fs)
legend('10% of rest')

%
subplot(ns, 1, 8)
spikeMat = M.spikelog(p.in+1:p.full, :);

% out_cells = 1:100;
% 
% sum(C([mems{1}], :))>3
% spikeMat = spikeMat(randperm(size(spikeMat,1)),:);
% spikeMat = spikeMat(1:p.full/10, :);
for trialCount = 1:size(spikeMat,1)
    hold all;
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        disp(trialCount)
        disp(spikePos)
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'Color', 'm');
        end
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('spike log','FontSize',fs)
legend('10% of rest')