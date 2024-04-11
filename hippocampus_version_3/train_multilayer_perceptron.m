function[error, W1, W2] = train_multilayer_perceptron(input, target_output)

%  Set some parameters
iterations  = 100;             % Number of training blocks
bias        = 0.5;             % Threshold
alpha       = 0.1;             % Learning rate
n_trials    = size(input,1);
% Define architecture 
input_layer_size = 6;
hidden_layer_size = 6;
output_layer_size = 1;

% Initialize weights
W1 = rand(hidden_layer_size, input_layer_size) - 0.5;
W2 = rand(output_layer_size, hidden_layer_size) - 0.5;

% Memory
error  = nan(1,n_trials*iterations);   % Error log
output = nan(1,n_trials*iterations);   % Network output log

for i = 1:iterations
    for t = 1:n_trials

        % [1] Forward pass

        % Input to hidden layer
        hidden_input =  (W1 * input(t, :)')' ;                   
        hidden_output = (1.0./(1.0 + exp(-hidden_input)));     % Sigmoid activation function
        
        % Hidden to output layer
        output_layer_input = W2 * hidden_output';
        output_output = 1.0./(1.0 + exp(-output_layer_input)); % Sigmoid activation function for output
   
        % Output is binary (either reward or no reward (1/0))
        output_output   = double(output_output>=bias);
        output((i-1)*n_trials+t) = output_output;               % update output log 

        % [2] Backward pass

        % Derivative of the sigmoid function for output and hidden layer
        output_deriv = output_output * (1 - output_output);
        hidden_deriv = hidden_output .* (1 - hidden_output);
    
        % Output layer error
        output_delta = (output_output - target_output(t)) * output_deriv;
        error((i-1)*n_trials+t) = output_delta;                 % update error log 

        % Hidden layer error
        hidden_delta = (W2' * output_delta')' .* hidden_deriv;
      
        % Update weights
        W2 = W2 - alpha * output_delta' * hidden_output;       % Adjusted weight update for W2
        W1 = W1 - alpha * (hidden_delta' * input(t, :))';      % Adjusted weight update for W1

    end
end

%% 
% 
% pints = 0:0.1:9.9;
% 
% performance = arrayfun(@(x) x - (x - 2)^2, pints);
% 
% figure;
% plot(pints, performance, 'LineWidth', 2);
% title('Pool Performance vs. Alcohol Intake (Pints)');
% 
% xlabel('Pints');
% 
% performance_labels = {'Leave the Pub', 'Inept', 'Embarrasing', 'Below Mediocre', 'Average', 'Passable', 'Kinda Alright', 'Decent?'};
% performance_positions = linspace(min(performance), max(performance), length(performance_labels));
% 
% yticks(performance_positions);
% yticklabels(performance_labels);
% 
% grid on;
