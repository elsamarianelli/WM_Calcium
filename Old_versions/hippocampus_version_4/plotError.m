function plotError(error, subplotIndex, color, performance_test, ncols)
% function to plot the error over perceptron iterations
    subplot(ncols, 2, subplotIndex);
    plot(fastsmooth(abs(error), 3000), 'Color', color);
    set(gca, 'FontSize', 8); 
    axis square;
    xlabel('Trial Number', 'FontSize', 8);
    ylabel('Moving Error', 'FontSize', 8);
    title(sprintf('%.2f', performance_test), "FontSize",8); % Assuming performance_test is a scalar
end
