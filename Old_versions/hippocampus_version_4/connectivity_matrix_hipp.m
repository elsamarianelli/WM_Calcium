function [C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles)
% Generates random connectivity matrix (C) and synaptic weight matrix J for
% CA3 --> CA1 feedforward connections in network based on preset coding
% level c. Connectivity can be constrained either such that:
% (i) every CA1 cell receives exactly the same number of inputs, but CA3
% cells can project to a variable number of outputs (fixIn  = true); or
% (ii) every CA3 cell sends exactly the same number of outputs, but CA1    
% cells can receive a variable number of inputs (fixIn = false)

% Assign some memory
fixIn       = true;
C           = zeros(p.in, p.out);            % Connectivity matrix
J           = ones(p.in, p.out).*p.j_b;      % Synaptic weight matrix     

% Set the random connectivity
if fixIn

    % Check that the coding level is such that there are sufficient input
    % neurons to generate projections from
    if p.c > p.f
        error('This level of synaptic connectivity cannot be supported by the coding level in the output layer')
    end

    % Generate an array of inputs that do not form part of any ensemble
    allIn   = 1:p.in;
    allIn([ca3_ensembles{:}]) = [];

    % Loop through each output cell
    for c   = 1 : p.out

        % Identify which odours it responds to, and therefore which CA3
        % populations it might have projections from
        inMems      = cellfun(@(x) ismember(c+p.in,x),ca1_ensembles);
        inMems      = unique(horzcat(ca3_ensembles{inMems}));

        % If it does not respond to any odours, select randomly from allOut
        if isempty(inMems)
            thisIn  = allIn(randperm(length(allIn)));
        else
            thisIn  = inMems(randperm(length(inMems)));
        end
        C(thisIn(1:round(p.in*p.c)), c) = 1; clear inMems thisIn

    end
    clear c

else

    % Check that the coding level is such that there are sufficient output
    % neurons to project to
    if p.c > p.f_o
        error('This level of synaptic connectivity cannot be supported by the coding level in the output layer')
    end

    % Generate an array of outputs that do not form part of any ensemble
    allOut  = p.in + (1:p.out);
    allOut([ca1_ensembles{:}]-p.in) = [];

    % Loop through each input cell
    for c   = 1 : p.in

        % Identify which odours it responds to, and therefore which CA1
        % populations it might project to
        inMems      = cellfun(@(x) ismember(c,x),ca3_ensembles);
        inMems      = unique(horzcat(ca1_ensembles{inMems}));

        % If it does not respond to any odours, select randomly from allOut
        if isempty(inMems)
            thisOut = allOut(randperm(length(allOut)));
        else
            thisOut = inMems(randperm(length(inMems)));
        end
        C(c, thisOut(1:round(p.out*p.c))-p.in) = 1; clear inMems thisOut

    end
    clear c
        
end

% Set the strength of potentiated connections above baseline
J(C==1)     = p.j_p;