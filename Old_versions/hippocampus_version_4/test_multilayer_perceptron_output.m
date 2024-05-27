function[performance_test] = test_multilayer_perceptron_output(data, w1, w2)
%% Function to test weights in trained perceptron on newly generated test data
%  Extract some parameters, set some parameters, assign some memory
n_trials   	= size(data,1);              % Number of trials
n_ca1       = size(data,2)-1;            % Number of CA1 inputs
output      = nan(1,n_trials);           % Network output (lick / no lick)
error       = nan(1,n_trials);           % Error
%  Run the dynamics

% Sigmoid activation function
sigmoid = @(x) 1 ./ (1 + exp(-x));

for t	= 1 : n_trials

    % input to hidden
    hidden_input = data(t, 1:n_ca1) * w1;
    hidden_output = sigmoid(hidden_input);

    % hidden to output
    final_input = hidden_output * w2;
    o = double(final_input >= 0);            % Binary threshold (reward/no reward)

    % Error calculation
    d = data(t, n_ca1+1) - o;
    output(t)	= o; clear o    
    error(t)   = d; clear d   

end
clear t
clear i

performance_test = length(find(error == 0))/length(error);
