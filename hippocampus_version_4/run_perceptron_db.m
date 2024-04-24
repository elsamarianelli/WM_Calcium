function[output,error, w] = run_perceptron_db(data)
%% Function to iteratively train a perceptron on simulated CA1 output

%  Extract some parameters, set some parameters, assign some memory
iterations  = 500;                                  % Number of training blocks
bias        = 0.5;                                  % Threshold
alpha       = 0.01;                                  % Learning rate
n_trials   	= size(data,1);                         % Number of trials
n_ca1       = size(data,2)-1;                       % Number of CA1 inputs
w           = zeros(n_ca1,1);                       % Connection weights
output      = nan(1,n_trials*iterations);           % Network output (lick / no lick)
error       = nan(1,n_trials*iterations);           % Error

%  Run the dynamics
for i       = 1 : iterations
    for t	= 1 : n_trials
        o	= data(t,1:n_ca1)*w;                    % Compute output
        o   = double(o>=bias);                      % Threshold
        d   = data(t,n_ca1+1) - o;                  % Error term
        w   = w + alpha * d .* data(t,1:n_ca1)';	% Update weights
        output((i-1)*n_trials+t)	= o; clear o    % Record output
        error((i-1)*n_trials+t)     = d; clear d    % Record error
    end
    clear t
    disp(['iteration', i])
end
clear i