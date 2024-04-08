%% are my odours linearly sperable?
% if I assume that the activity of certain subpopulations gives me
% information about which of the odours were presented in which order (1st
% is smaller subpopulation of overlap and second is broader population
% which is always active) then I can try to treat the problem as a linear
% seperability problem to see if it's solvable using a single layer
% perceptron as I have done 


% Plot the pairs
marker_size = 100;

figure;
scatter(2, 2, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'CB (reward)'); hold on;
scatter(1, 1, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'BA (reward)');
scatter(3, 3, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'AC (reward)');
scatter(2, 1, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'AB (no reward)');
scatter(1, 3, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'CA (no reward)');
scatter(3, 2, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'BC (no reward)');

xlabel('Active odour populations');
ylabel('Active overlapping populations');
labels_x = {'A', 'B', 'C'};
labels_y = {'1', '2', '3'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels_x, 'YTick', 1:3, 'YTickLabel', labels_y);
ylim([0 4]); xlim([0 4])
title('subpopulation activity');
legend('Location', 'best');
% grid on;

%% or alternatively based on activation strength

% columns are A, B, C, AB, BC, CA
a = [1 2 0 3 2 1
    0 1 2 1 3 2
    2 0 1 2 1 3
    2 1 0 3 1 2
    0 2 1 2 3 1
    1 0 2 1 2 3];

x_reward = 
% Plot the pairs
marker_size = 100;

figure;
scatter(2, 2, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'CB (reward)'); hold on;
scatter(1, 1, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'BA (reward)');
scatter(3, 3, marker_size, 'filled', 'MarkerFaceColor', 'b', 'DisplayName', 'AC (reward)');
scatter(2, 1, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'AB (no reward)');
scatter(1, 3, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'CA (no reward)');
scatter(3, 2, marker_size, 'filled', 'MarkerFaceColor', 'r', 'DisplayName', 'BC (no reward)');

xlabel('Active subpopulations');
ylabel('Active subpopulations');
labels_x = {'A', 'B', 'C'};
labels_y = {'1', '2', '3'};
set(gca, 'XTick', 1:3, 'XTickLabel', labels_x, 'YTick', 1:3, 'YTickLabel', labels_y);
ylim([0 4]); xlim([0 4])
title('subpopulation activity');
legend('Location', 'best');
% grid on;
