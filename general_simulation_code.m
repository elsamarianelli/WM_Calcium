%% General simulation code

% Code to replicate results in Figure 2 of Mongillo et al's. 2008 
% Synaptic Theory of Working Memory Paper, which models working memory as a 
% property of calcium-mediated short term synaptic facilitation, in a 
% recurrent network of integrate-and-fire neurones.

% Version 2 
% (Elsa Marianelli, contactable at zcbtetm@ucl.ac.uk)

%% current problems 
% 1) memory cells firing for whole duration of reactivaiton signal unlike
% in paper figure
% 2) anything past 0.85 multiplication of excitatory signal essentially
% breaks model --> why does recurrent system have to be balanced out at the
% beggining 
% 3) set up to have 3 memories overlapping
% 4) not sure time from initial or reactition length makes a difference...
%% set up
% Define simulation parameters (p) and get input times (in)... 
[p, in]  = get_params(3, ...   %to multiply number of neurones by
                    2000, ...  %simulation length (ms)
                    1.05, ...  %selective stimulation contrast factor
                    1.015,...  %reactivating signal contrast factor
                    0.85, ...  %to multiply excitatory signal 
                    .9,   ...  % SD of external current sigma
                    100,  ...  %initial stimulation length 
                    10);      %reactivation length

% Generate memory for network (M)
[M] = get_memory(p);

% Assign Neurons to memories (mems)
% generate synaptic connectivity matrix (C) and synaptic strength matrix (J)...
[C, J, mems, first_input, second_input] = connectivity_matrix(M, p, 0.25, 'AC'); %0.25 is the degree of overlap with each 
                                                                                 % initial memory pattern(to overlap odour C with odours A and B)
%% Simulate the dynamics

%  should think about having two options - one where we log everything (for debugging, which uses a lot more memory) and one where we only keep the
%  crucial outputs (but maybe we can deal with this later and just log everything for now)
%initialise variables
SE          = J;
u           = p.U.*(ones(p.Ne, p.Ne));
x           = ones(p.Ne, p.Ne); 
ref_p       = zeros(p.N,1);
delays      = randi([1 5], p.N, p.N).*C;
delay_window      = zeros(p.N, p.N, 5);

for t       = 1 : p.SimLength
    % get updated synaptic efficacy for excitatory-to-excitatory connections (which display short term plasticity)
    M.U_log(:,:, t) = u;
    M.X_log(:,:, t) = x;  
    
    % get external inputs all the time (change strength to 1 for paper set up)
    M.Iext(1:p.Ne) = normrnd(p.mu_e,p.sigma, p.Ne, 1).*p.ex_fact;  
    M.Iext(p.Ne+1:end) = normrnd(p.mu_i, p.sigma, p.Ni, 1).*p.ex_fact; %(3)

    % present memory with a contrast factor applied if in stimulation
    % period or reactivation period
    activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.simulation});
    M.Iext(horzcat(mems{activeMems})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.Ne).*p.SCF;                 
    
    % for when the reloading the same memory 
    % activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.reactivation});
    % M.Iext(horzcat(mems{activeMems})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.Ne).*p.RCF;                 
    % M.Iext_log(:, t) = M.Iext;
    % for when reloading a different memory 
    activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.reactivation});
    if activeMems == 1
        M.Iext(horzcat(mems{2})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.Ne).*p.SCF; %EM: why horzcat? why sum activae mems if its a logical   
    else
    end
                                    % x by 2 to generate index for memory
                                    % pattern used for odour B? not really
                                    % sure how this has been re written so
                                    % need to check
    M.Iext_log(:, t) = M.Iext;

    % The membrane potential in the next time step is dictated by its value in this timestep, and the externally applied current
    M.V_log(:, t) = M.Irec;
    tau_const   = 1./[p.tau_e*ones(p.Ne,1) ; p.tau_i*ones(p.Ni,1)];
    current     = -M.V + M.Irec + M.Iext;                              %(2)
    M.V         = M.V  + (tau_const .*current);

    % Check if the Neuron has exceeded threshold - i.e. fired a spike and is not in refractory period
    % EM: think this order makes more sense for updating RP? check with dan
    fired       = M.V > p.SPE & ref_p==0;
    % update refractory periods for this timestep
    ref_p(ref_p~=0) = ref_p(ref_p~=0) - 1; 
    % update refractory periods for cells which fire - V needs to be set to subrethsold reset potential Vr when refractory period ends
    ref_p(fired) = p.ARP; 
    % reset potentials for cells which fired
    M.V(fired) = [p.Vr_e*ones(sum(fired(1:p.Ne)),1) ; p.Vr_i*ones(sum(fired(p.Ne+1:end)),1)];

    % updating synapse log into future based on transmission delays
    pre     = nonzeros(repmat(fired',p.N,1) .* repmat((1:p.N),p.N,1) .*C);
    post    = nonzeros(repmat(fired',p.N,1) .* repmat((1:p.N)',1,p.N) .*C);
    arrives = nonzeros(repmat(fired',p.N,1) .* delays .*C); %+t
    %update next 5 time step spike logs, rolling window through spike log
    %to save memory EM: ?
    spike_at_synapse = delay_window(:, :, 1);
    delay_window(:, :, 1) = [];
    delay_window = cat(3, delay_window, zeros(p.N, p.N));
    delay_window(post',pre', arrives') = 1;
    
    % reccurent currents calculated as the sum over inputting neurones
    % which fired by the efficacy at that synapse                     (4)
    % spike_at_synapse = M.synapse_log(:, :, t);
    M.Irec    =  SE.*spike_at_synapse;
    M.Irec    =  sum(M.Irec, 2);  
    M.Irec_log(:, t) = M.Irec;

    J_fired = double(fired)'; 
    M.J_fired_log(:, t) = J_fired';
    M.spikelog(:, t) = fired; 
    clear fired;
    % update u and x parameters for E to E connections and update synaptic
    % efficacy matrix
    f_or_s = spike_at_synapse(1:p.Ne, 1:p.Ne);
    u = u + (((p.U -u)/p.tau_facil) + (p.U.*(1-u).*f_or_s));          %(5)
    x = x + (((1-x)/p.tau_decay) - (u.*x.*f_or_s));                   %(6)
    SE_hat = J(1:p.Ne, 1:p.Ne).*u.*x;
    % synaptic efficacy for E to E plastic synapses
    SE(1:p.Ne, 1:p.Ne) = SE_hat;                                      %(7)
    disp(t)
end 

%% plotting output
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

subplot(ns,1,2)
av_u = mean(M.U_log, 1); av_u = squeeze(av_u);
idx = mems{1}; 
av_u_memory = av_u(idx, :); 
av_u_memory = mean(av_u_memory, 1);

av_x = mean(M.X_log, 1); av_x = squeeze(av_x);
idx = mems{1}; 
av_x_memory = av_x(idx, :); 
av_x_memory = mean(av_x_memory, 1);

plot(1:p.SimLength,av_u_memory,'b')
ylabel(first_input,'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
%u and x for second memory 
subplot(ns,1,3)
av_u = mean(M.U_log, 1); av_u = squeeze(av_u);
idx = mems{2}; 
av_u_memory = av_u(idx, :); 
av_u_memory = mean(av_u_memory, 1);

av_x = mean(M.X_log, 1); av_x = squeeze(av_x);
idx = mems{2}; 
av_x_memory = av_x(idx, :); 
av_x_memory = mean(av_x_memory, 1);

plot(1:p.SimLength,av_u_memory,'b')
ylabel(second_input,'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')

%plot recurrent current log 
subplot(ns, 1, 4)
av_i = mean(M.Irec_log, 1); 
plot(1:p.SimLength, av_i)
ylabel('recurent','FontSize',fs)
subplot(ns, 1, 5)
av_e = mean(M.Iext_log, 1); 
plot(1:p.SimLength, av_e)
ylabel('external','FontSize',fs)

%plotting raster
color_ops = { 'k', 'r'};

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
            [trialCount-0.4 trialCount+0.4], 'g');
        end
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('spike log','FontSize',fs)
