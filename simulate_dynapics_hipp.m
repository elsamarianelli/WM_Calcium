function M = simulate_dynapics_hipp(p, C, J, input, M, mems)
    
    %  Initialise variables
    u           = p.U.*(ones(p.in, p.out));
    x           = ones(p.in, p.out);
    ref_p       = zeros(p.full,1);
    fired       = zeros(p.full, 1);
    %  Loop through time
    for t       = 1 : p.SimLength
    
        % [1] Deal with external and recurrent synaptic inputs
    
        % First, generate random external input
        M.Iext(1:p.in) = normrnd(p.mu_e,p.sigma, p.in, 1);

        % Second, input mem 1 at specified times
        activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{input.simulation});
        M.Iext(horzcat(mems{activeMems})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.in).*p.CF; clear activeMems;
        
        % and mem2
        activeMems      = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{input.reactivation});
        if activeMems   == 1
            M.Iext(horzcat(mems{2})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.in).*p.CF;
        end
        clear activeMems
    
        % % Third, get synaptic currents (from CA3 to CA1)
        M.Irec              = sum((J.*x.*u),1); % (4)
    
        % Finally, update short-term plasticity parameters for E to E connections
        u(1:p.in,1:p.out) = u(1:p.in,1:p.out)  + (((p.U -u(1:p.in,1:p.out))/p.tau_facil) + (p.U.*(1-u(1:p.in,1:p.out)).*fired(1:p.in).*C));
        x(1:p.in,1:p.out) = x(1:p.in,1:p.out) + (((1-x(1:p.in,1:p.out) )/p.tau_decay) - (u(1:p.in,1:p.out) .*x(1:p.in,1:p.out) .*fired(1:p.in).*C));

        % [2] Update membrane potential and check for spiking
    
        % Update membrane potential
        M.V(1:p.in)                = M.V(1:p.in)         + ((-M.V(1:p.in)         + M.Iext) ./(p.tau_e*ones(p.in,1)));
        M.V(p.in+1:p.full)         = M.V(p.in+1:p.full)  + ((-M.V(p.in+1:p.full)  + M.Irec') ./(p.tau_e*ones(p.in,1)));

        % Check for spiking
        fired               = M.V > p.SPE & ref_p == 0;
    
        % Update refractory periods for neurons that fired
        ref_p(ref_p~=0)     = ref_p(ref_p~=0) - 1;
        ref_p(fired)        = p.ARP;
    
        % Reset membrane potential for cells which fired
        M.V(fired)          = p.Vr_e*ones(sum(fired(1:p.full)),1);

      
        % % [4] Log what needs logging

        % log synaptic parameters for each memory
        mem_syn_inds = sub2ind(size(u),repmat(mems{1},1,length(mems{1})),sort(repmat(mems{1},1,length(mems{1}))));
        mem_syn_inds(C(mem_syn_inds)==1);
        M.U_mem1_log(t) = mean(u(mem_syn_inds(C(mem_syn_inds)==1)));
        M.X_mem1_log(t) = mean(x(mem_syn_inds(C(mem_syn_inds)==1))); 

        mem_syn_inds = sub2ind(size(u),repmat(mems{2},1,length(mems{2})),sort(repmat(mems{2},1,length(mems{2}))));
        mem_syn_inds(C(mem_syn_inds)==1);
        M.U_mem2_log(t) = mean(u(mem_syn_inds(C(mem_syn_inds)==1)));
        M.X_mem2_log(t) = mean(x(mem_syn_inds(C(mem_syn_inds)==1)));
 
        % log external input
        M.Iext_log(:, t) = M.Iext;
        % log recurrent input
        M.Irec_log(:, t) = M.Irec;
        M.spikelog(:, t) = fired; 

    end
end