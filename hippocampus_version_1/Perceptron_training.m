%% code to run perceptron with an online learning algorithm 
%  MATLAB function training was being done using batch learning (aka updates
%  weights with cumulative errors after each epoch) 
%  This is to do online learning (aka updating weights after each trial
%  presentation)
%  and just to have a better idea of what's ouput in general.

%% get input
% programmable paramaters 
degree_overlap = 0.2;
pattern_order = 'AB';
length_first = 30;
lenght_second = 30;
delay_time = 500;
start_time = 200;
% Get non-programmable paramaters
p = get_params_hipp(0.85);
% Get connectivity matrix and synaptic efficacy matrix
[C, J] = connectivity_matrix_hipp(p);
% Times that memory is 'on', ms
input.simulation = [start_time (start_time+length_first)];
input.reactivation = [(start_time+lenght_second+delay_time) (start_time+lenght_second+lenght_second+delay_time)];
% generate memory
M = get_memory_hipp(p);

%% generate training data 
n_trials = 6.*100;
data_overlap = get_train_data(C, J, input, n_trials, 0.2, p);
data_no_overlap = get_train_data(C, J, input, n_trials, 0.0, p);
data = data_no_overlap;

% partition data
train_data = data(1:(n_trials*0.8), :);
test_data = data((n_trials*0.8)+1:end, :);

% seperate into input and expected ouput 
x = train_data(1:(n_trials*0.8), 1:p.out); 
y = train_data(1:(n_trials*0.8), end);
% shuffle_y = y(randperm(length(y)));
% y = shuffle_y;

w = ones(size(x, 2)+1, 1)';
% alpha sets the speed of learning
alpha = 0.001;
n = size(x,1);
error_log = [];

% training perceptron weights
err = 1;
bias = 1; % Bias value

while err > 0
    % initialize accumulated error
    err = 0;
    % for every set of inputs
    for i = 1:n

        % compute the actual output
        o = dot(w, [x(i,:), bias]) >= 0; % Adding bias term
        
        % compute error
        e = y(i) - o;

        % update weights if the actual output does not equal the expected
        % output (if abs(e)>0)
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
run = 1:size(accuracy, 2);
figure
plot(accuracy)
hold on
% testing ouput on test data
x_test = test_data(1:(n_trials*0.2), 1:p.out); 
y_test = test_data(1:(n_trials*0.2), end);



% Assuming x_test is your test data matrix and y_test is the corresponding labels

% Initialize success counter
success_count = 0;

% Test the perceptron on the test data
for i = 1:size(x_test, 1)
    % Compute actual output
    o = dot(w, [x_test(i,:), bias]) >= 0; % Adding bias term
    
    % Compare to labeled output
    if o == y_test(i)
        success_count = success_count + 1;
    end
end

% Calculate accuracy
accuracy = success_count / size(x_test, 1);
disp('Perceptron Accuracy:')
disp(accuracy);

% % running perceptron on test data to see how much of the time it gives
% % correct output
% success_log = [];
% for i = 1:size(x_test, 1)
%     % compute actual ouput
%     o = 1;
%     o(dot(w,[x_test(i, :), bias])<0) = 0;
%     % compare to labelled ouput an log result
%     check = abs(o-y_test(i));
%     success_log = [success_log, check];
% end
% 
% rate_success = sum(success_log)./size(success_log, 2);
