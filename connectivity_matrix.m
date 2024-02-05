function [C, J, mems] = connectivity_matrix(M, p)

J = M.J; C = M.C;
cells       = 1 : p.Ne;                       % List of all Neurons
mems        = cell(p.p_m,1);                  % Array of Neurons in each memory
for i       = 1 : p.p_m

    % Assign a random selection of remaining unassigned Neurons to that
    % memory
    cells   = cells(randperm(length(cells)));
    mems{i} = sort(cells(1:p.f*p.Ne));
    cells(ismember(cells,mems{i})) = [];
    
    % Generate synaptic connections from that memory to all other Neurons
    for j   = 1 : p.N
        inputs                  = mems{i}(~ismember(mems{i},j));
        inputs                  = inputs(randperm(length(inputs)));
        C(inputs(1:p.c*length(mems{i})),j) = 1; clear inputs
    end
    clear j

    % Set the weight of 'potentiated', 'baseline', and exc->inh synaptic 
    % connections from Neurons in that memory to others in the network
    J(mems{i},1:p.Ne)             = p.J_b * C(mems{i},1:p.Ne);
    J(mems{i},mems{i})          = p.J_p * C(mems{i},mems{i});
    J(mems{i},p.Ne+1:p.N)           = p.J_ie * C(mems{i},p.Ne+1:p.N);
    
end
clear i

%  Generate synaptic connections from unassigned Neurons to each cell, and
%  from inhibitory Neurons to each cell
i_cells     = p.Ne+1 : p.N;
for i       = 1 : p.N
    inputs  = cells(~ismember(cells,i));
    inputs  = inputs(randperm(length(inputs)));
    C(inputs(1:p.c*(1-p.f*p.p_m)*p.Ne),i) = 1; clear inputs
    inputs  = i_cells(~ismember(i_cells,i));
    inputs  = inputs(randperm(length(inputs)));
    C(inputs(1:p.c*p.Ni),i)         = 1; clear inputs
end
clear i_cells i

%  Set synaptic weights from unassigned Neurons to all other cells
J(cells,1:p.Ne)                     = p.J_b * C(cells,1:p.Ne);
J(cells,p.Ne+1:p.N)                 = p.J_ie * C(cells,p.Ne+1:p.N); clear cells 

%  Set synaptic weights for each of the inhibitory connections
J(p.Ne+1:p.N,1:p.Ne)                    = p.J_ei * C(p.Ne+1:p.N,1:p.Ne);
J(p.Ne+1:p.N,p.Ne+1:p.N)                = p.J_ii * C(p.Ne+1:p.N,p.Ne+1:p.N);

end