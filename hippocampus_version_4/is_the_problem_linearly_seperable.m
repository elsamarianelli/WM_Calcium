%% can the perceptron (single layer) solve the problem in the first place
% if I assume that the activity of certain subpopulations gives me
% information about which of the odours were presented in which order (1st
% is smaller subpopulation of overlap and second is broader population
% which is always active) then I can try to treat the problem as a linear
% seperability problem to see if it's solvable using a single layer
% perceptron as I have done 

% Plot the pairs
marker_size = 100;
colour_reward =  [0, 128, 128] / 255;
colour_no_reward = [.7 .7 .7];
figure;
scatter(2, 2, marker_size, 'filled', 'MarkerFaceColor', colour_reward, 'DisplayName', 'CB (reward)'); hold on;
scatter(1, 1, marker_size, 'filled', 'MarkerFaceColor', colour_reward, 'DisplayName', 'BA (reward)');
scatter(3, 3, marker_size, 'filled', 'MarkerFaceColor', colour_reward, 'DisplayName', 'AC (reward)');
scatter(2, 1, marker_size, 'filled', 'MarkerFaceColor', colour_no_reward, 'DisplayName', 'AB (no reward)');
scatter(1, 3, marker_size, 'filled', 'MarkerFaceColor', colour_no_reward, 'DisplayName', 'CA (no reward)');
scatter(3, 2, marker_size, 'filled', 'MarkerFaceColor', colour_no_reward, 'DisplayName', 'BC (no reward)');

xlabel('Active odour populations');
ylabel('Active overlapping populations');
labels_x = {'Ano', 'Bno', 'Cno'};
labels_y = {'1', '2', '3'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels_x, 'YTick', 1:3, 'YTickLabel', labels_y);
ylim([0 4]); xlim([0 4])
title('subpopulation activity');
legend('Location', 'best');
grid on;

%% try to solve perceptron based on simplified data 
% 0 = not firing, 1 = fired in first odour but not since then, 
% 2 = fired in second odour but wasnt potentiated, 3 = potentiated and
% fired in second odour

column_labels = {'A', 'B', 'C', 'e1', 'e2', 'e3'};

% event_data = [1 2 0 3 2 1
%         0 1 2 1 3 2
%         2 0 1 2 1 3
%         2 1 0 3 1 2
%         0 2 1 2 3 1
%         1 0 2 1 2 3];
% 
target_output = [1 1 1 0 0 0]; % where 1 = reward, 0 = no reward
% input_data = [event_data, target_output'];

% data simpler 

input = [0 0 1 0 0 1
              0 1 0 0 1 0 
              1 0 0 1 0 0 
              0 1 0 1 0 0 
              0 0 1 0 1 0 
              1 0 0 0 0 1];
input_data = [input, target_output'];

% run perceptorn just with this simplified input for many trials to see if
% it can learn weights
%[performance_accuracy, error, w]    = run_perceptron_db(input_data);
[error, W1, W2] = train_multilayer_perceptron(input, target_output);
%%  Debug plot (requires fastsmooth function)
figure;
plot(fastsmooth((error),10)), set(gca,'FontSize',18), axis square
xlabel('Trial Number','FontSize',24), ylabel('Moving Average Error','FontSize',24)

[performance_test] = test_perceptron_output(input_data, w);
box off;