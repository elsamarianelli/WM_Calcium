%% General simulation code
%  Code to implement a model of working memory using "silent" short-term
%  facilitation, as described in:
%  Mongillo et al. (2008) Science 319: 1543-1546
%
%  Version 2
%  Elsa Marianelli, contactable at zcbtetm@ucl.ac.uk

%% Issues
%  1) memory cells firing for whole duration of reactivaiton signal unlike
%  in paper figure
%  2) not sure time from initial or reactition length makes a
%  difference...response is too immediate


%% Define simulation parameters (p) and get input times (in)
[p, in] = get_params(3, ...     % to multiply number of neurones by
                    15000, ...   % simulation length (ms)
                    1.15, ...   % selective stimulation contrast factor
                    1.05,...    % reactivating signal contrast factor
                    0.85, ...   % to multiply excitatory signal 
                    .84,   ...    % SD of external current sigma
                    800, ...   % time between 1st and 2nd pattern 
                    150,  ...   % initial stimulation length 
                    60);        % reactivation length

%% Assign memory, randomly assign neurons to memories, generate synaptic 
%  connectivity matrices
M       = get_memory(p);        % Generate memory for network (M)

%  Assign neurons to memories (mems), generate synaptic connectivity matrix 
%  (C) and synaptic strength matrix (J)
[C, J, mems, first_input, second_input] = connectivity_matrix(p, 1/8, 'AA'); 

%% Simulate dynamics
fire_rate_means = [];
Vm_max_reacts = [];
for delay_time = 500:500:15000
    [M_new, fire_rate_mean,fire_rate_mean_others, Vm_max_react] = simulate_dynamics(p, in, M, C, J, mems, first_input, second_input, delay_time);
    fire_rate_means = [fire_rate_means, fire_rate_mean];
    fire_rate_mean_others_all = [fire_rate_mean_others_all, fire_rate_mean_others];
    Vm_max_reacts = [Vm_max_reacts, Vm_max_react];
    disp(delay_time, fire_rate_means, fire_rate_mean_others_all)
end
% % 
M = M_new;
% % 
% plot(1:1000:13000, fire_rate_means)
% % plot(1:1000:13000, Vm_max_reacts)
% xlabel('time after initial simulation network is pinged')
% % ylabel('peak Vm reached by memory  during reactivation')
% ylabel('mean firing rate in memory  during reactivation')

% n_runs = 50;
% data = zeros(p.N, n_runs.*2);   
% % 1/x = the degree to which the neural populations overlap '--' can put in desired pattern order, AA for repeated patterns
% for i = 1:n_runs  
%     [M_new, slice] = simulate_dynamics(p, in, M, C, J, mems, first_input, second_input);
%     data(:, i) =slice;  
%     disp(i)
% end
% [C, J, mems, first_input, second_input] = connectivity_matrix(M, p, 1/8, 'AC'); 
% for i = n_runs+1:n_runs*2
%     [M_new, slice] = simulate_dynamics(p, in, M, C, J, mems, first_input, second_input);
%     data(:, i) =slice;  
%     disp(i)
% end

%% Plot the output, if required
% plotting vm in both patterns

fs = 10;
ns = 6;
subplot(ns, 1, 1)
for m=1:2
    x = M.V_log(mems{m}, :);
    x(x==0)=NaN;
    x = nanmean(x, 1);
    x(isnan (x)) = 0;
    plot(1:p.SimLength, x)
    hold on;
    x_points = [in.simulation(1), in.simulation(1), in.simulation(2), in.simulation(2)];
    x2_points =[in.reactivation(1), in.reactivation(1), in.reactivation(2), in.reactivation(2)];
    y_points = [0, max(x), max(x), 0];
    color = [0, 0, 1];
    hold on;
    a = fill(x_points, y_points, color,'LineStyle','none');
    a.FaceAlpha = 0.1; 
    hold on;
    a = fill(x2_points, y_points, color,'LineStyle','none');
    a.FaceAlpha = 0.1; 
    hold on;
end
ylabel('Vm in memory','FontSize',fs)

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
%plot recurrent current log adn external input 
subplot(ns, 1, 4)
av_i = mean(M.Irec_log, 1); 
plot(1:p.SimLength, av_i)
ylabel('recurent','FontSize',fs)
subplot(ns, 1, 5)
av_e = mean(M.Iext_log, 1); 
plot(1:p.SimLength, av_e)
ylabel('external','FontSize',fs)

%plotting spike raster
color_ops = { 'b', 'r'};

subplot(ns, 1, 6)
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
                plot([spikePos(spikeCount) spikePos(spikeCount)], ...
                [trialCount-0.4 trialCount+0.4], color_ops{m});
            end
        end
    end
    hold on;
end

spikeMat = M.spikelog; spikeMat(mems{1},:) = [];spikeMat(mems{2},:) = [];
spikeMat = spikeMat(randperm(size(spikeMat,1)),:);
spikeMat = spikeMat(1:p.N/10, :);
for trialCount = 1:size(spikeMat,1)
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'Color', [.5 .5 .5]);
        end
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('spike log','FontSize',fs)