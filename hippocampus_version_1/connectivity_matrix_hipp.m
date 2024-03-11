function [C, J] = connectivity_matrix_hipp(p)

    %% connectivity matrix
    C = zeros(p.in,p.out);                       % Connectivity matrix
    J = ones(p.in, p.out).*p.j_b;

    for i = 1:p.in
        rand = randperm(p.out);   
        C(i, rand(1:p.out*p.c)) = 1;
    end
    
    J(C==1) = p.j_p;                              %synaptic weight matrix
end
