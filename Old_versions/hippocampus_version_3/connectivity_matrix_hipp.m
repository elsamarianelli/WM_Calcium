function [C, J] = connectivity_matrix_hipp(p, ca3_ensembles, ca1_ensembles)
% generates random connectivity matrix (C) for CA3 --> CA1 feedforward
% connections in network based on preset coding level c, and synaptic
% weight matrix J. Connectivity can be constrained either such that:
% (i) every CA1 cell receives exactly the same number of inputs, but CA3
% cells can project a variable number of outputs (fixIn  = true); or
% (ii) every CA3 cell sends exactly the same number of outputs, but CA1    
% cells can receive a variable number of inputs (fixIn = false)

% Assign some memory
fixIn       = true;
% C           = zeros(p.out,p.in);            % Connectivity matrix
% J           = ones(p.out,p.in).*p.j_b;      % Synaptic weight matrix
C           = zeros(p.in, p.out);            % Connectivity matrix          %% E.M. I think that the in and out are the wrong way round in the connectivity matrix
J           = ones(p.in, p.out).*p.j_b;      % Synaptic weight matrix     

% Set the random connectivity
if fixIn

    % Check that the coding level is such that there are sufficient output
    % neurons to project to
    if p.c > p.f_o
        error('This level of synaptic connectivity cannot be supported by the coding level in the output layer')
    end

    % Generate an array of outputs that do not form part of any ensemble
    allOut  = p.in + (1:p.out);
    allOut(horzcat(cell2mat(ca1_ensembles))-p.in) = [];     %% E.M. before there was ca1_ensembles{:} but this only indexed the 3rd cell

    % Loop through each input cell
    for c   = 1 : p.in

        % Identify which odours it responds to, and therefore which CA1
        % populations it might project to
        inMems      = cellfun(@(x) ismember(c,x),ca3_ensembles);
        inMems      = unique(horzcat(ca1_ensembles{inMems}));

        % If it does not respond to any odours, select randomly from allOut
        if isempty(inMems)
            rand    = allOut(randperm(length(allOut)));
        else
            rand    = inMems(randperm(length(inMems)));
        end
        C(c, rand(1:round(p.out*p.c))-p.in) = 1; clear inMems rand  %% E.M. c comes first, input first? - should change p.in

    end
    clear c

else    %% E.M. need to change this

    % Check that the coding level is such that there are sufficient input
    % neurons to generate projections from
    if p.c > p.f
        error('This level of synaptic connectivity cannot be supported by the coding level in the output layer')
    end

    % Generate an array of inputs that do not form part of any ensemble
    allIn   = 1:p.in;
    allIn(horzcat(cell2mat(ca1_ensembles))-p.in) = [];

    % Loop through each output cell
    for c   = 1 : p.out

        % Identify which odours it responds to, and therefore which CA3
        % populations it might have projections from
        inMems      = cellfun(@(x) ismember(c+p.in,x),ca1_ensembles);
        inMems      = unique(horzcat(ca3_ensembles{inMems}));

        % If it does not respond to any odours, select randomly from allOut
        if isempty(inMems)
            rand    = allIn(randperm(length(allIn)));
        else
            rand    = inMems(randperm(length(inMems)));
        end
        C(rand(1:round(p.in*p.c)), c) = 1; clear inMems rand    % E.M. flipped this aswell 

    end
    clear c

end

% Set the strength of potentiated connections above baseline
J(C==1)     = p.j_p;