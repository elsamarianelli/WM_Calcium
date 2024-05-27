function M      = simulate_dynamics_hipp(p, C, J, input, M, mems)
%% Code to simulate activity of a CA1 (output) population driven by CA3 input and calcium
%  facilitation in CA3 --> CA1 connections, for a single trial in which two
%  input patterns are presented to CA3 at specified times. Calcium and
%  synaptic dynamics are both implemented using Euler integration.
%  Transmission delays are currently assumed to be a single time step (1ms)
%  for simplicity.

% Identify synapses of interest, to log data for           
pre_overlap = intersect(mems{1}, mems{2});
[i,j]       = find(C(pre_overlap,:)>0);
mem_ovlp    = sub2ind(size(C),pre_overlap(i)',j); clear i j
mem1_unq    = mems{1}(~ismember(mems{1},pre_overlap));
[i,j]       = find(C(mem1_unq,:)>0);
mem1_unq    = sub2ind(size(C),mem1_unq(i)',j); clear i j
mem2_unq    = mems{2}(~ismember(mems{2},pre_overlap)); clear pre_overlap
[i,j]       = find(C(mem2_unq,:)>0);
mem2_unq    = sub2ind(size(C),mem2_unq(i)',j); clear i j

%  Loop through time
for t       = 1 : p.SimLength

    %% [1] Deal with external and feed-forward synaptic inputs

    % First, generate random external input
    M.Iext              = normrnd(p.mu_e,p.sigma, p.full, 1);

    % Second, input mem 1 at specified times
    activeMems          = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{input.simulation ; input.reactivation});
    M.Iext(horzcat(mems{activeMems})) = M.Iext(horzcat(mems{activeMems})).*p.CF; clear activeMems

    % Third, get synaptic currents (from CA3 to CA1)
    M.Isyn              = [zeros(p.in,1) ; sum((J.*M.x.*M.u.*M.delay_win(:,:,1)),1)'];

    % Finally, update short-term plasticity parameters for E to E connections
    M.u                 = M.u + (((p.U - M.u)/p.tau_facil) + (p.U .* (1-M.u) .* M.delay_win(:,:,1) .*C));
    M.x                 = M.x + (((1 - M.x)/p.tau_decay) - (M.u .*M.x .* M.delay_win(:,:,1) .* C));

    %% [2] Update membrane potential and check for spiking

    % Update membrane potential
    M.V                 = M.V + ((-M.V + M.Iext + M.Isyn) ./(p.tau_e*ones(p.full,1)));

    % Check for spiking
    fired               = M.V > p.SPE & M.ref_p == 0;

    % Update refractory periods for neurons that fired
    M.ref_p(M.ref_p~=0) = M.ref_p(M.ref_p~=0) - 1;
    M.ref_p(fired)      = p.ARP;

    % Reset membrane potential for cells which fired
    M.V(fired)          = p.Vr_e*ones(sum(fired),1);


    %% [3] Send spiking outputs along the delay lines
    M.delay_win(:,:,1)  = 0;                              % Eliminate spikes that have just been transmitted
    M.delay_win         = circshift(M.delay_win,-1,3);    % Move all spikes one timestep closer

    % Update synapse log based on transmission delays
    pre                 = nonzeros(repmat(fired(1:p.in),1,p.out) .* repmat((1:p.in)',1,p.out) .* C);
    post                = nonzeros(repmat(fired(1:p.in),1,p.out) .* repmat((1:p.out),p.in,1) .* C);
    arrives             = nonzeros(repmat(fired(1:p.in),1,p.out) .* M.delays .* C);
    ind                 = sub2ind(size(M.delay_win), pre, post, arrives);
    M.delay_win(ind)    = 1; clear pre post arrives ind


    %% [4] Log what needs logging

    % log synaptic parameters for independent and overlapping
    % components of each memory
    M.U_mem1_log(t)     = mean(M.u(mem1_unq));
    M.X_mem1_log(t)     = mean(M.x(mem1_unq));
    M.U_mem2_log(t)     = mean(M.u(mem2_unq));
    M.X_mem2_log(t)     = mean(M.x(mem2_unq));
    M.U_ovlp_log(t)     = mean(M.u(mem_ovlp));
    M.X_ovlp_log(t)     = mean(M.x(mem_ovlp));

    % log external input
    M.spikelog(:, t)    = fired; clear fired

end
clear mem1_unq mem2_unq mem_ovlp t

end