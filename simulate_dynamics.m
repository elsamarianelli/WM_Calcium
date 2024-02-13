function [M, spike_slice] = simulate_dynamics(p, in, M, C, J, mems, first_input, second_input)                                                                                 

    %initialise variables     
    SE          = J;
    u           = p.U.*(ones(p.Ne, p.Ne));
    x           = ones(p.Ne, p.Ne); 
    ref_p       = zeros(p.N,1);
    delays      = randi([1 5], p.N, p.N).*C;
    delay_window      = zeros(p.N, p.N, 5);
    
    for t       = 1 : p.SimLength

        % [1] INPUTS

        % get external inputs all the time (change strength to 1 for paper set up)
        M.Iext(1:p.Ne) = normrnd(p.mu_e,p.sigma, p.Ne, 1).*p.ex_fact;  
        M.Iext(p.Ne+1:end) = normrnd(p.mu_i, p.sigma, p.Ni, 1).*p.ex_fact; %(3)
   
        % get initial input when presenting first memory
        activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.simulation});
        M.Iext(horzcat(mems{activeMems})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.Ne).*p.SCF.*p.ex_fact;              

        % pinging whole network with weak excitation to reactivate memory
        if strcmp(first_input, second_input)
            activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.reactivation});
            if activeMems == 1
                M.Iext(1:p.Ne) = normrnd(p.mu_e,p.sigma, p.Ne, 1).*p.ex_fact.*p.RCF;  
                M.Iext(p.Ne+1:end) = normrnd(p.mu_i, p.sigma, p.Ni, 1).*p.ex_fact.*p.RCF;               
                M.Iext_log(:, t) = M.Iext;
            end
        else
        % or, if different pattern is being presented, activate that pattern with the original contrast factor
            activeMems = cellfun(@(x) any(x(:,1)<=t & x(:,2)>=t),{in.reactivation});
            if activeMems == 1
                M.Iext(horzcat(mems{2})) = normrnd(p.mu_e,p.sigma,1,sum(activeMems)*p.f*p.Ne).*p.RCF.*p.ex_fact;
            end
        end
        
        % log external input
        M.Iext_log(:, t) = M.Iext;
    
        % reccurent currents calculated as the sum over inputting neurones which fired by the efficacy at that synapse                     
        spike_at_synapse = delay_window(:, :, 1);
        M.Irec    =  SE.*spike_at_synapse;                                 %(4)
        M.Irec    =  sum(M.Irec, 2);  
        % log recurrent input
        M.Irec_log(:, t) = M.Irec;
        M.V_log(:, t) = M.Irec;


        % [2] SPIKING 
        
        % The membrane potential in the next time step is dictated by its value in this timestep, and the externally applied current
        tau_const   = 1./[p.tau_e*ones(p.Ne,1) ; p.tau_i*ones(p.Ni,1)];
        current     = -M.V + M.Irec + M.Iext;                              %(2)
        M.V         = M.V  + (tau_const .*current);
        
        % Check if the Neuron has exceeded threshold - i.e. fired a spike and is not in refractory period, log 
        fired       = M.V > p.SPE & ref_p==0;
        J_fired = double(fired)'; 
        M.J_fired_log(:, t) = J_fired';
        M.spikelog(:, t) = fired;

        % update refractory periods 
        ref_p(ref_p~=0) = ref_p(ref_p~=0) - 1; 
        ref_p(fired) = p.ARP; 

        % reset membraine potentials for cells which fired
        M.V(fired) = [p.Vr_e*ones(sum(fired(1:p.Ne)),1) ; p.Vr_i*ones(sum(fired(p.Ne+1:end)),1)];
        

        % [3] OUTPUT - send the spikes forward in time, update your short term plasticity variable

        % updating synapse log into future based on transmission delays
        pre     = nonzeros(repmat(fired',p.N,1) .* repmat((1:p.N),p.N,1) .*C);
        post    = nonzeros(repmat(fired',p.N,1) .* repmat((1:p.N)',1,p.N) .*C);
        arrives = nonzeros(repmat(fired',p.N,1) .* delays .*C); %+t
       

        delay_window = circshift(delay_window,-1,3); 
        delay_window(:, :, 5) = 0;
        ind = sub2ind([size(delay_window)], post, pre, arrives);
        delay_window(ind) = 1;
        % delay_window(post, pre, arrives) = 1;
% 
        % update u and x parameters for E to E connections and use to update synaptic efficacy
        u = u + (((p.U -u)/p.tau_facil) + (p.U.*(1-u).*spike_at_synapse(1:p.Ne, 1:p.Ne)));          %(5)
        x = x + (((1-x)/p.tau_decay) - (u.*x.*spike_at_synapse(1:p.Ne, 1:p.Ne)));                   %(6)
        SE_hat = J(1:p.Ne, 1:p.Ne).*u.*x; SE(1:p.Ne, 1:p.Ne) = SE_hat;                              %(7)
        
        % log synaptic paramaters for each memory
        av_u_1 = (mean(u, 1)); M.U_mem1_log(t) = mean(av_u_1(mems{1}));  
        av_u_2 = (mean(u, 1)); M.U_mem2_log(t) = mean(av_u_2(mems{2}));  
        av_x_1 = (mean(x, 1)); M.X_mem1_log(t) = mean(av_x_1(mems{1}));  
        av_x_2 = (mean(x, 1)); M.X_mem2_log(t) = mean(av_x_2(mems{2}));  
        
        disp(t)


    end 
    
    spike_slice = M.V_log(:, in.reactivation(2));
end