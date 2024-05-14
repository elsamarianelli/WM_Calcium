function plotError(error, subplotIndex, color, performance_test, ncols)
% function to plot the error over perceptron iterations
    subplot(ncols, 2, subplotIndex);
    plot(fastsmooth(abs(error), 3000), 'Color', color);
    set(gca, 'FontSize', 12); 
    axis square;
    xlabel('Trial Number', 'FontSize', 10);
    ylabel('Moving Average Error', 'FontSize', 10);
    title(sprintf('Performance = %.2f', performance_test), "FontSize",12); % Assuming performance_test is a scalar
end