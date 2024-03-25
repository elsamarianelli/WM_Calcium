function [mean_firing_second_odour] = get_mean_firing_second_odour(p, C, J, input, M, mems, lenght_second)
%plotting mean spiking in overlapping cells vs non overlapping cells
%during second odour, filtering for cells which recived increasing number of inputs

figure

for filter = 1:3
    
    ind1 = find((sum(C(mems{1}, :))>filter)); 
    ind2 = find((sum(C(mems{2}, :))>filter)); 
    coi = intersect(ind1,ind2);
    full = 1:p.out; full(coi) = [];
    
    n_trials = 10;
    
    overlapping_log = zeros(n_trials, 1);
    non_overlapping_log = zeros(n_trials, 1);
    
    for i = 1:n_trials

        M = simulate_dynapics_hipp(p, C, J, input, M, mems);
    
        overlapping = M.spikelog(p.in + coi, input.reactivation(1):input.reactivation(2));
        non_overlapping = M.spikelog(p.in+full,  input.reactivation(1):input.reactivation(2));
    
        mean_overlapping = sum( overlapping, "all")/(size(overlapping, 1).*lenght_second.*0.001);
        mean_non_overlapping = sum(non_overlapping, "all")/(size(non_overlapping, 1).*lenght_second.*0.001);
    
        overlapping_log(i) = mean_overlapping; non_overlapping_log(i) = mean_non_overlapping;
        disp(i)

    end
    
    
    % plot non overlapping cells
    if filter == 1 
        non_over = mean(non_overlapping_log);
        scatter((1*ones(n_trials)), non_overlapping_log,'MarkerEdgeColor','r', 'MarkerFaceColor','r')
        hold on;
        plot([1-0.1 1.+0.1], [non_over non_over], 'Color', 'r', 'LineWidth',2);
        hold on;
    end
    over = mean(overlapping_log); 
    scatter(filter+(1*ones(n_trials)),overlapping_log, 'MarkerEdgeColor','b', 'MarkerFaceColor','b');
    hold on;
    plot([filter+1-0.1 filter+1+0.1], [over over],'Color', 'b', 'LineWidth',2);

end
xticks([1 2 3 4 5 6])
xticklabels({'non-overlap', '>1 input', '>2 input', '>3 input', '>4 input', '>5 input'})

xlabel('CA1 cells receiving non overlapping input and overlapping input from CA3')
ylabel('mean firing rate during second odour presentation (Hz)')

L1 = plot(nan, nan, 'color', 'r');
L2 = plot(nan, nan, 'color', 'b');
legend([L1, L2 ], {'non-overlapping' 'overlapping'}, 'Location', 'southeast')

mean_firing_second_odour = gcf;
end