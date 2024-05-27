function performance_accuracy = run_perceptron(data, n_trials, p)
%% code to run perceptron with an online learning algorithm 
%  MATLAB function training was being done using batch learning (updates
%  weights with cumulative errors after each epoch) 
%  This is to do online learning (updating weights after each trial
%  presentation)
%  and just to have a better idea of what's ouput in general.

    % partition data
    train_data = data(1:(n_trials*0.8), :);
    test_data = data((n_trials*0.8)+1:end, :);
    
    % seperate into input and expected ouput 
    x = train_data(1:(n_trials*0.8), 1:p.out); 
    y = train_data(1:(n_trials*0.8), end);
    % shuffle_y = y(randperm(length(y)));
    % y = shuffle_y;
    
    %% train weights and bias on train data
    w = ones(size(x, 2)+1, 1)';
    % alpha sets step size (speed of learning)
    alpha = 0.001;
    n = size(x,1);
    error_log = [];
    
    % training perceptron weights
    err = 1;
    bias = 1; % Bias value (to shift decision boundary)
    
    while err > 0
        % initialize accumulated error
        err = 0;
        % for every set of inputs
        for i = 1:n
    
            % compute the actual output
            o = dot(w, [x(i,:), bias]) >= 0; % Adding bias term
            
            % compute error
            e = y(i) - o;
    
            % update weights if the actual output does not equal the expected output (if abs(e)>0)
            w = w + alpha * e * [x(i,:), bias]; % Updating weights with bias
            
            % accumulate error
            err = err + abs(e);

        end
        disp('run')
        disp(err)
        error_log = [error_log, err];
    end
    
    % plot learning over runs 
    accuracy = (error_log/n);
    % run = 1:size(accuracy, 2);
    plot(accuracy)
    hold on
    
    % legend('delay of 400', 'delay of 2000')
    % xlabel('run')
    % ylabel('error')
    
    %% assessing performance on test data
    x_test = test_data(1:(n_trials*0.2), 1:p.out); 
    y_test = test_data(1:(n_trials*0.2), end);
    
    % Initialize success counter
    success_count = 0;
    
    % Test the perceptron on the test data
    for i = 1:size(x_test, 1)
        
        % Compute actual output
        o = dot(w, [x_test(i,:), bias]) >= 0; 
        
        % Compare to labeled output
        if o == y_test(i)
            success_count = success_count + 1;
        end

    end
    
    % Calculate accuracy
    performance_accuracy = success_count / size(x_test, 1);

    % disp('Perceptron Accuracy:')
    % disp(accuracy);
end