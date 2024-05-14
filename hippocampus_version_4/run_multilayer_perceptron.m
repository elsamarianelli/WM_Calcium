function [output, error, w1, w2] = run_multilayer_perceptron(data)
    % Function to iteratively train a perceptron with one hidden layer

    % Parameters
    iterations = 1000;                               % Number of training iterations
    alpha = 0.1;                                   % Learning rate
    n_trials = size(data,1);                        % Number of trials
    n_ca1 = size(data,2) - 1;                       % Number of CA1 inputs
    n_hidden = 6;                                   % Number of neurons in the hidden layer
    bias = 0.5;                                     % Bias term for weights
    momentum = 0.9;                                 % momentum term for w1 updates (to try to prevent getting stuck in local minima)
    delta_w1 = 0;                                   % initial change in w1                                                                                                               
    % Initialize weights
    w1 = rand(n_ca1, n_hidden) - bias;              % Weights from input to hidden layer
    %w2 = rand(n_hidden, 1) - bias;                  % Weights from hidden to output layer

    % what if i already know what my second layer weights should look like
    % --> makes the performance better
    no_reward_weights = -1 + 0.2 * randn(1, 3);
    reward_weights = 1 + 0.2 * randn(1, 3);
    w2 = [no_reward_weights, reward_weights]';

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

            % Calculate the gradient for each weight in w1
            sigmoid_derivative = hidden_output .* (1 - hidden_output);
            gradient_w1 = (data(t, 1:n_ca1)' * (d * (w2' .* sigmoid_derivative)));

            % Update weights w1
            momentum_w1 = delta_w1 * momentum;
            delta_w1 = (alpha * gradient_w1);
            w1 = w1 + delta_w1 + momentum_w1;

            % Store output and error for each trial
            output((i - 1) * n_trials + t) = o;
            error((i - 1) * n_trials + t) = d;

        end
        disp(['iteration ', num2str(i)])
    end
end

