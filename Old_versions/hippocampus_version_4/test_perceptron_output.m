function[performance_test] = test_perceptron_output(data, w)
%% Function to test weights in trained perceptron on newly generated test data

%  Extract some parameters, set some parameters, assign some memory
bias        = 0.5;                                  % Threshold
n_trials   	= size(data,1);                         % Number of trials
n_ca1       = size(data,2)-1;                       % Number of CA1 inputs
output      = nan(1,n_trials);           % Network output (lick / no lick)
error       = nan(1,n_trials);           % Error
%  Run the dynamics

for t	= 1 : n_trials
    o	= data(t,1:n_ca1)*w;                    % Compute output
    o   = double(o>=bias);                      % Threshold
    d   = data(t,n_ca1+1) - o;   
    % Error term
    output(t)	= o; clear o    % Record output
    error(t)   = d; clear d    % Record error
end
clear t
clear i

performance_test = length(find(error == 0))/length(error);
