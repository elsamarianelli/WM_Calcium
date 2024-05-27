function M      = get_memory_hipp(p)

%% Assign memory
M.V             = p.Vr_e * ones(p.full,1);                    % Membrane voltage
M.Iext          = zeros(p.full,1);                            % External input to each neuron
M.Isyn          = zeros(p.full,1);                            % Synaptic input to each neuron
M.u             = p.U.*(ones(p.in, p.out));                   % Short term facilitation
M.x             = p.X.*ones(p.in, p.out);                     % Synaptic resources
M.ref_p         = zeros(p.full,1);                            % Refractory period
M.delays        = randi([1 p.max_delay], p.in, p.out);        % Axonal / synaptic connection delays
M.delay_win     = zeros(p.in, p.out, p.max_delay);

% Initialise logging variables
M.spikelog      = zeros(p.full, p.SimLength);
M.U_mem1_log    = zeros(1, p.SimLength);
M.X_mem1_log    = zeros(1, p.SimLength);
M.U_mem2_log    = zeros(1, p.SimLength);
M.X_mem2_log    = zeros(1, p.SimLength);
M.U_ovlp_log    = zeros(1, p.SimLength);
M.X_ovlp_log    = zeros(1, p.SimLength);


end