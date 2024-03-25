function [mems] = get_odours_hipp(p)

    %% memory input
    cells       = 1 : p.in;                   % List of all neurons
    mems        = cell(2,1);                  % Array of neurons in each memory
    
    % generating seperate patterns for memory input 1 and 2 for 3 different
    % odours, order can be choses in function
    
    %pattern A
    num_overlap_cells   = round((p.f*p.in).*p.degree_overlap);
    cells               = cells(randperm(length(cells)));
    pattern_A           = sort(cells(1:round(p.f*p.in)));
    cells(ismember(cells,pattern_A)) = [];
    
    % pattern B
    AB_overlap          = pattern_A(1:num_overlap_cells);
    cells               = cells(randperm(length(cells)));
    non_overlap         = cells(1:(round(p.f*p.in))-num_overlap_cells);
    pattern_B           = sort([AB_overlap, non_overlap]);    
    
    mems{1}             = eval(strcat('pattern', '_', p.pattern_order(1)));
    mems{2}             = eval(strcat('pattern', '_', p.pattern_order(2)));


end
