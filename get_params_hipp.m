function p = get_params_hipp(factor)

    p.in      = 200;                          % Number of CA3 (input) cells
    p.out     = 200;                          % Number of CA1 (ouput) cells
    p.full    = p.in+p.out;                   % total number of neurones
    p.f       = 0.2;                          % Coding level
    p.c       = 0.2;                         % Probability of synaptic contact
    
    p.j_p     = 2;                            % potentiated level of E-->E synapses
    p.j_b     = 0;                            % baseline level of E-->E synapses 
    p.SPE     = 20;                           % spike emmision threshold
    p.Vr_e    = 17;                           % reset potentials
    p.tau_e   = 14;                           % membrane time constant
    p.ARP     = 2;                            % absolute refractory period
    
    p.mu_e =  23.10.*factor;                  %mean external current
    p.sigma = 0.87*factor;                    %standard deviation of external current
    p.SimLength = 1500;                       %length of simulation (ms)
    p.tau_decay = 200;                        %recovery of synaptic resources 
    p.tau_facil = 800;                       %recovery time of utilisation factor 
    p.U = 0.2;                                %Baseline utilisation factor
    p.X = 1;                                  %Baseline synaptic materials factor
    p.CF = 1.1;                               % stimulation contrast factor
    p.RF = 1.1;                             %second odour stim factor (check)
    
end
