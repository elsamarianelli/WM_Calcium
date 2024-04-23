function p      = get_params_hipp(p)
%% generates a data structure with all parameters needed for the model
%  (most are taken from Mongillo et al. Science 2008)

    p.in        = 200;              % Number of CA3 (input) cells
    p.out       = 100;              % Number of CA1 (ouput) cells
    p.full      = p.in+p.out;       % Total number of neurones
    p.f         = 0.2;              % Coding level of input layer
    p.f_o       = 0.2;              % Coding level of ouput layer
    p.c         = 0.2;              % Probability of synaptic contact
    
    p.max_delay = 1;                % Maximum CA3 -> CA1 axonal / synaptic delay
    p.j_p       = 1;                % potentiated strength of E-->E synapses    
    p.j_b       = 0;                % background strength of E-->E synapses 
    p.SPE       = 20;               % spike emission threshold
    p.Vr_e      = 16;               % reset potentials
    p.tau_e     = 15;               % membrane time constant
    p.ARP       = 2;                % absolute refractory period
    
    p.mu_e      = 23.10*p.scaleF;   % Mean external current
    p.sigma     = 0.87*p.scaleF;    % Standard deviation of external current
    p.tau_decay = 200;              % Recovery of synaptic resources 
    p.tau_facil = 1500;             % Recovery time of utilisation factor 
    p.U         = 0.2;              % Baseline utilisation factor
    p.X         = 1;                % Baseline synaptic materials factor
    p.CF        = 1.05;             % Stimulation contrast factor

end
