function [mems] = get_odours_hipp(p, degree_overlap, pattern_order)

    %% memory input
    cells       = 1 : p.in;                   % List of all neurons
    mems        = cell(2,1);                  % Array of neurons in each memory
    
    % generating seperate patterns for memory input 1 and 2 for 3 different
    % odours, order can be choses in function
    
    %pattern A
    num_overlap_cells = (p.f*p.in).*degree_overlap;
    cells       = cells(randperm(length(cells)));
    pattern_A   = sort(cells(1:p.f*p.in));
    cells(ismember(cells,pattern_A)) = [];
    % pattern B
    AB_overlap  = pattern_A(1:num_overlap_cells);
    cells       = cells(randperm(length(cells)));
    non_overlap = sort(cells(1:(p.f*p.in)-num_overlap_cells));
    pattern_B   = [AB_overlap, non_overlap];
    cells(ismember(cells,pattern_B)) = [];
    % pattern C 
    AC_overlap  = pattern_A(end-num_overlap_cells+1:end);
    BC_overlap  = pattern_B(end-num_overlap_cells+1:end);
    cells       = cells(randperm(length(cells)));
    non_overlap = sort(cells(1:int16((p.f*p.in)*(1-degree_overlap*2))));
    pattern_C   = [AC_overlap, BC_overlap, non_overlap];
    
    mems{1}     = eval(strcat('pattern', '_', pattern_order(1)));
    mems{2}     = eval(strcat('pattern', '_', pattern_order(2)));


end
