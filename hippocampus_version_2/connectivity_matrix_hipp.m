function [C, J] = connectivity_matrix_hipp(p, overlap_control, mems_all)
% generates random connectivity matrix (C) for CA3 --> CA1 feedforward 
% connections in network based on preset coding level c, 
% and Synaptic weight matrix J. 
% turn ON to control CA1 populations which overlap and then generate
% connectivity BACK to CA3 randomly, 
% turn OFF to have original CA3 to CA1 random connectivity where overlap is
% set in CA1

if overlap_control == "OFF"
% CA3 overlap controlled

    C = zeros(p.in,p.out);            % Connectivity matrix 
    J = ones(p.in, p.out).*p.j_b;

    for i = 1:p.in
        rand = randperm(p.out);   
        C(i, rand(1:p.out*p.c)) = 1;
    end
    
    J(C==1) = p.j_p;                  % Synaptic weight matrix

elseif overlap_control == "ON"
% CA1 overlap controlled
    CA3_populations = mems_all;
    % get CA1 odour populations which CA3 odour populations will be allowed
    % to connect to 
    CA1_populations = get_odours_hipp(p, p.degree_overlap_CA1, "ON");

    C = zeros(p.in,p.out);            % Connectivity matrix
    J = ones(p.in, p.out).*p.j_b;
    
    % Connectivity matrix where the populations of CA1 that CA3 populations
    % are allowed to project into are constrained 
    for population = 1:length(CA3_populations)
        % get corresponding subpopulations for each odour
        subpop_CA3 = CA3_populations{population};
        subpop_CA1 = CA1_populations{population};
        for i = 1:length(subpop_CA3)
            % randomly project to subsection of CA1 cell
            perm = randperm(length(subpop_CA1));
            subpop_CA1_shuff = subpop_CA1(perm);
            C(subpop_CA3(i), subpop_CA1_shuff(1:(length(subpop_CA1))*p.c*4)) = 1;
        end
    end

    J(C==1) = p.j_p;                  % Synaptic weight matrix

end
