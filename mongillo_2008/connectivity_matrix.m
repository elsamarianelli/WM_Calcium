function [C, J, mems, first_input, second_input] = connectivity_matrix(p, degree_overlap, pattern_order)

J           = zeros(p.N,p.N);                   % Connectivity matrix
C           = zeros(p.N,p.N);                   % Synaptic weight matrix
cells       = 1 : p.Ne;                         % List of all neurons
mems        = cell(p.p_m,1);                    % Array of neurons in each memory

% generating seperate patterns for memory input 1 and 2 for 3 different
% odours, order can be choses in function

%pattern A
num_overlap_cells = (p.f*p.Ne).*degree_overlap;
cells       = cells(randperm(length(cells)));
pattern_A   = sort(cells(1:p.f*p.Ne));
cells(ismember(cells,pattern_A)) = [];
% pattern B
AB_overlap  = pattern_A(1:num_overlap_cells);
cells       = cells(randperm(length(cells)));
non_overlap = sort(cells(1:(p.f*p.Ne)-num_overlap_cells));
pattern_B   = [AB_overlap, non_overlap];
cells(ismember(cells,pattern_B)) = [];
% pattern C 
AC_overlap  = pattern_A(end-num_overlap_cells+1:end);
BC_overlap  = pattern_B(end-num_overlap_cells+1:end);
cells       = cells(randperm(length(cells)));
non_overlap = sort(cells(1:(p.f*p.Ne)*(1-degree_overlap*2)));
pattern_C   = [AC_overlap, BC_overlap, non_overlap];

mems{1}     = eval(strcat('pattern', '_', pattern_order(1)));
mems{2}     = eval(strcat('pattern', '_', pattern_order(2)));
first_input = pattern_order(1); second_input = pattern_order(2);

for i       = 1 : p.p_m
    
    % Generate synaptic connections from that memory to all other Neurons
    for j   = 1 : p.N
        inputs                  = mems{i}(~ismember(mems{i},j));
        inputs                  = inputs(randperm(length(inputs)));
        C(inputs(1:round(p.c*length(mems{i}))),j) = 1; clear inputs
    end
    clear j

    % Set the weight of 'potentiated', 'baseline', and exc->inh synaptic 
    % connections from Neurons in that memory to others in the network
    J(mems{i},1:p.Ne)           = p.J_b * C(mems{i},1:p.Ne);
    J(mems{i},mems{i})          = p.J_p * C(mems{i},mems{i});
    J(mems{i},p.Ne+1:p.N)       = p.J_ie * C(mems{i},p.Ne+1:p.N);
    
end
clear i

%  Generate synaptic connections from unassigned Neurons to each cell, and
%  from inhibitory Neurons to each cell
i_cells     = p.Ne+1 : p.N;
for i       = 1 : p.N
    inputs  = cells(~ismember(cells,i));
    inputs  = inputs(randperm(length(inputs)));
    C(inputs(1:round(p.c*(1-p.f*p.p_m)*p.Ne)),i) = 1; clear inputs
    inputs  = i_cells(~ismember(i_cells,i));
    inputs  = inputs(randperm(length(inputs)));
    C(inputs(1:p.c*p.Ni),i)     = 1; clear inputs
end
clear i_cells i

%  Set synaptic weights from unassigned Neurons to all other cells
J(cells,1:p.Ne)                 = p.J_b * C(cells,1:p.Ne);
J(cells,p.Ne+1:p.N)             = p.J_ie * C(cells,p.Ne+1:p.N); clear cells 

%  Set synaptic weights for each of the inhibitory connections
J(p.Ne+1:p.N,1:p.Ne)            = p.J_ei * C(p.Ne+1:p.N,1:p.Ne);
J(p.Ne+1:p.N,p.Ne+1:p.N)        = p.J_ii * C(p.Ne+1:p.N,p.Ne+1:p.N);

end