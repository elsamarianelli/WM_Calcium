%% Assign memory
function M = get_memory(p)

M.V           = [p.Vr_e * ones(p.Ne,1) ; p.Vr_i * ones(p.Ni,1)];  % Membrane voltage
M.Iext        = zeros(p.N,1);                   % Ip.Nput to each Neuron
M.Irec        = zeros(p.N,1);                   % Input to each Neuron
M.C           = zeros(p.N,p.N);                   % Conp.Nectivity matrix
M.J           = zeros(p.N,p.N);                   % Synaptic weight matrix

M.U_log       = zeros(p.Ne,p.Ne, p.SimLength);
M.X_log       = zeros(p.Ne,p.Ne, p.SimLength);
M.spikelog    = zeros(p.N, p.SimLength);
M.V_log       = zeros(p.N, p.SimLength);
M.Irec_log    = zeros(p.N, p.SimLength);
M.Iext_log    = zeros(p.N, p.SimLength);
M.J_fired_log = zeros(p.N, p.SimLength);
M.synapse_log = zeros(p.N, p.N, p.SimLength+5); 

M.U_mem1_log  = zeros(1, p.SimLength);
M.U_mem2_log  = zeros(1, p.SimLength);
M.X_mem1_log  = zeros(1, p.SimLength);
M.X_mem2_log  = zeros(1, p.SimLength);
end