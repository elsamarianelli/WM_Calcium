function [p, in] = get_params(factor, simlength, SCF, RCF, ex_fact, sigma, delay_time, sim_time, react_time)
%function to get all parameters and input tiems for memories as given in
%Mongillo et al 2007 SOM paper
    %% parameters
    p.Ne      = 80.*factor;                   % Number of excitatory neurons
    p.Ni      = 20.*factor;                   % Number of inhibitory neurons
    p.N       = p.Ne + p.Ni;                  % Total number of neurons
    p.f       = 0.1;                          % Coding level
    p.p_m     = 2;                            % Number of memories (non-overlapping)
    p.c       = 0.2;                          % Probability of synaptic contact
    p.J_ie    = 0.135;                        % synaptic efficacy E-->I
    p.J_ei    = 0.25;                         % synaptic efficacy I-->E
    p.J_ii    = 0.20;                         % synaptic efficacy I-->I
    p.J_b     = 0.10;                         % baseline level of E-->E synapses
    p.J_p     = 0.45;                         % potenitated level of E-->E synapses
    p.SPE     = 20;                           % spike emmision threshold
    p.Vr_e    = 16; p.Vr_i = 13;              % reset potentials
    p.tau_e   = 15; p.tau_i = 10;             % membrane time constant
    p.ARP     = 2;                            % absolute refractory period
    p.mu_e =  23.10; p.mu_i = 21;             %mean external current
    p.sigma = sigma;                          %standard deviation of external current
    p.SimLength = simlength;                  %length of simulation (ms)
    p.tau_decay = 200;                        %recovery of synaptic resources 
    p.tau_facil = 1500;                       %recovery time of utilisation factor 
    p.U = 0.2;                                %Baseline utilisation factor
    p.ex_fact = ex_fact;                      %factor by which excitatory input is multiplied by 
    p.SCF = SCF;                              %selective stimulation contrast factor
    p.RCF = RCF;                              %reactivating signal contrast factor
    %% input details 
    % times that memory is 'on', ms
    in.simulation = [200 (200+sim_time)];
    in.reactivation = [200+delay_time (200+delay_time+react_time)];
end
