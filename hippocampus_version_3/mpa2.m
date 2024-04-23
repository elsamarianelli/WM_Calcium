% Sample data generation for demonstration
Y = [1 1 1 0 0 0]; % where 1 = reward, 0 = no reward
X = [0 0 1 0 0 1
              0 1 0 0 1 0 
              1 0 0 1 0 0 
              0 1 0 1 0 0 
              0 0 1 0 1 0 
              1 0 0 0 0 1];
data = [X Y'];
function [output, error, w1, w2] = run_perceptron_db(data)
%% Function to iteratively train a perceptron with one hidden layer on simulated CA1 output

% Extract some parameters, set some parameters, assign some memory
iterations  = 1000;                                 % Number of training blocks
alpha       = 0.01;                                 % Learning rate
n_trials    = size(data,1);                        % Number of trials
n_ca1       = size(data,2)-1;                      % Number of CA1 inputs
n_hidden    = 6;                                   % Number of neurons in the hidden layer

% Initialize weights
w1          = rand(n_ca1, n_hidden) - 0.05;  % Connection weights from input to hidden layer
w2          = rand(n_hidden, 1) - 0.05;      % Connection weights from hidden to output layer
output      = nan(1, n_trials * iterations);       % Network output (lick / no lick)
error       = nan(1, n_trials * iterations);       % Error

% Sigmoid activation function
sigmoid = @(x) 1 ./ (1 + exp(-x));

% Run the dynamics
for i = 1:iterations
    for t = 1:n_trials
        % Forward pass: from input to hidden
        hidden_input = data(t, 1:n_ca1) * w1;      % Linear combination
        hidden_output = sigmoid(hidden_input);     % Non-linear activation

        % Forward pass: from hidden to output
        final_input = hidden_output * w2;          % Linear combination
        o = double(final_input >= 0);              % Threshold (binary output)

        % Compute the error
        d = data(t, n_ca1+1) - o;                  % Error term

        % Backward pass: update weights using the error
        % Update w2 (hidden to output)
        delta_w2 = alpha * d * hidden_output';
        w2 = w2 + delta_w2;

        % Update w1 (input to hidden)
        % Compute derivative for sigmoid at hidden outputs
        sigmoid_derivative = hidden_output .* (1 - hidden_output);
        delta_w1 = alpha * d * w2' .* sigmoid_derivative .* data(t, 1:n_ca1);
        w1 = w1 + delta_w1';

        % Record output and error
        output((i - 1) * n_trials + t) = o;
        error((i - 1) * n_trials + t) = d;
    end
end

end

