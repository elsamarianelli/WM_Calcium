function [mems] = get_odours_hipp(cells,coding_level,degree_overlap)
% Function to generate 3 lists of cells which will be activated by each odour
% with variable degree of overlap

%% memory input
nCells      = length(cells);
nOverlap    = round((coding_level*nCells).*degree_overlap);
mems        = cell(3,1);                  % Array of neurons in each memory

% Pattern A
cells       = cells(randperm(nCells));
pattern_A   = cells(1:coding_level*nCells);
cells(ismember(cells,pattern_A)) = [];

% Pattern B
AB_overlap  = pattern_A(1:nOverlap);
cells       = cells(randperm(length(cells)));
non_overlap = cells(1:(coding_level*nCells)-nOverlap);
pattern_B   = [AB_overlap, non_overlap];
cells(ismember(cells,pattern_B)) = [];

% Pattern C
AC_overlap  = pattern_A(end-nOverlap+1:end);
BC_overlap  = pattern_B(end-nOverlap+1:end);
cells       = cells(randperm(length(cells)));
non_overlap = cells(1:round((coding_level*nCells)-nOverlap*2));
pattern_C   = [AC_overlap, BC_overlap, non_overlap];

mems{1}     = pattern_A; clear pattern_A
mems{2}     = pattern_B; clear pattern_B
mems{3}     = pattern_C; clear pattern_C

end
