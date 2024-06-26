function [output, error, w1, w2] = run_multilayer_perceptron(data)
    % Function to iteratively train a perceptron with one hidden layer

    % Parameters
    iterations = 1000;                               % Number of training iterations
    alpha = 0.1;                                    % Learning rate
    n_trials = size(data,1);                        % Number of trials
    n_ca1 = size(data,2) - 1;                       % Number of CA1 inputs
    n_hidden = 6;                               % Number of neurons in the hidden layer
    bias = 0.5;                                     % Bias term for weights

    % Initialize weights
    w1 = rand(n_ca1, n_hidden) - bias;              % Weights from input to hidden layer
    w2 = rand(n_hidden, 1) - bias;                  % Weights from hidden to output layer
    output = nan(1, n_trials * iterations);         % Network output
    error = nan(1, n_trials * iterations);          % Error tracking

    % Sigmoid activation function
    sigmoid = @(x) 1 ./ (1 + exp(-x));

    for i = 1:iterations
        for t = 1:n_trials

            % [1] Forward pass

            % input to hidden
            hidden_input = data(t, 1:n_ca1) * w1;
            hidden_output = sigmoid(hidden_input);

            % hidden to output
            final_input = hidden_output * w2;
            o = double(final_input >= 0);            % Binary threshold (reward/no reward)

            % Error calculation
            d = data(t, n_ca1+1) - o;

            %[2] Backward pass

            % Update output layer weights
            delta_w2 = alpha * d * hidden_output';
            w2 = w2 + delta_w2;

            % Update hidden layer weights
            sigmoid_derivative = hidden_output .* (1 - hidden_output);

            % Calculate the gradient for each weight in w1
            gradient_w1 = (data(t, 1:n_ca1)' * (d * (w2' .* sigmoid_derivative)));
            
            % Update weights w1
            w1 = w1 + (alpha * gradient_w1);

            % Store output and error
            output((i - 1) * n_trials + t) = o;
            error((i - 1) * n_trials + t) = d;

        end
    end
end

