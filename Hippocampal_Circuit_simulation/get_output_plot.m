function[fig] = get_output_plot(M, p, mems, C)

figure
fs = 10;
ns = 4; 
n = 0;
pattern_order = p.pattern_order;
% plotting synaptic parameters u an x for synapses potentiated by first
% odour which are reactivated during second odour presentation
ax1 = subplot(ns,1,n+1);
av_u_memory = M.U_mem2_log;
av_x_memory = M.X_mem2_log;
plot(1:p.SimLength,av_u_memory,'b')
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','best')
ylim([0 1])
ylabel(pattern_order(2),'FontSize',fs)
box off

%u and x for non overlapping 2nd odor cells
ax2 = subplot(ns,1,n+2);
av_u_memory = M.U_ovlp_log;
av_x_memory = M.X_ovlp_log;
plot(1:p.SimLength,av_u_memory,'b')
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','best')
ylim([0 1])
ylabel(strcat(pattern_order(1), ' overlap ', pattern_order(2)))
box off

%plot  current logs
% ax3 = subplot(ns, 1, n+3);
% av_e = mean(M.Iext_log(1:p.in, :), 1); 
% plot(1:p.SimLength, av_e)
% ylabel('ext-->CA3','FontSize',fs)
% 
% ax4 = subplot(ns, 1, n+4);
% av_i = mean(M.Irec_log(1:p.out, :), 1); 
% plot(1:p.SimLength, av_i)
% ylabel('CA3-->CA1','FontSize',fs)
% 
% ax5 = subplot(ns, 1, n+5);
% av_e = mean(M.Iext_log(p.in+1:p.full, :), 1); 
% plot(1:p.SimLength, av_e)
% ylabel('ext-->CA1','FontSize',fs)

%plotting spike raster for CA3 cells odour 1 and odour 2
color_ops = { 'b', 'r'};
ops = {0, 2};
ax6 = subplot(ns, 1, n+3);
for m = 1:2
    tVec = 1:p.SimLength;
    spikeMat = M.spikelog(mems{m}, :);
    hold all; 
    for trialCount = 1:size(spikeMat,1)
        if sum(spikeMat(trialCount, :)) == 0
            continue
        else
            spikePos = tVec(find(spikeMat(trialCount, :)));
            for spikeCount = 1:length(spikePos)
                plot([spikePos(spikeCount)+ops{m} spikePos(spikeCount)+ops{m}], ...
                [trialCount-0.4 trialCount+0.4], color_ops{m});
                
            end
        end
    end
    hold on;
end

spikeMat = M.spikelog; spikeMat(mems{1},:) = []; spikeMat(mems{2},:) = [];
spikeMat = spikeMat(randperm(size(spikeMat,1)),:);
spikeMat = spikeMat(1:p.full/10, :);
for trialCount = 1:size(spikeMat,1)
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'Color', 'm');
        end
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 p.SimLength])
ylabel('CA3','FontSize',fs)
box off

% 
% L1 = plot(nan, nan, 'color', 'b');
% L2 = plot(nan, nan, 'color', 'r');
% L3 = plot(nan, nan, 'color', 'm');
% legend([L1, L2, L3], {'odour 1' 'odour 2', '10% over non-odour cells'}, 'Location', 'best')

% plotting spike raster for CA1 ouput cell firing
ind1 = find((sum(C(mems{2}, :))>p.in.*0.01)); 
ind2 = find((sum(C(mems{1}, :))>p.in.*0.01)); 
coi = intersect(ind1,ind2);
full = 1:p.out; full(coi) = [];

% overlapping with more than 2 inputs
ax7 = subplot(ns, 1, n+4);
spikeMat = M.spikelog(p.in+coi, :);
for trialCount = 1:size(spikeMat,1)
    hold all;
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount)+1 spikePos(spikeCount)+1], ...
            [trialCount-1 trialCount+1], 'Color', 'm');
        end
    end
end
hold on;
% others (non overlapping and those receiving <2 inputs)
spikeMat = M.spikelog(p.in+full, :);
for trialCount = 1:size(spikeMat,1)
    hold all;
    if sum(spikeMat(trialCount, :)) == 0
        continue
    else
        spikePos = tVec(find(spikeMat(trialCount, :)));
        for spikeCount = 1:length(spikePos)
            plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-1 trialCount+1], 'Color', 'b');
        end
    end
end
% ylim([0 p.in]);
xlim([0 p.SimLength])
ylabel('CA1','FontSize',fs)
box off


% L1 = plot(nan, nan, 'color', 'm');
% L2 = plot(nan, nan, 'color', 'b');
% legend([L1, L2], {'overlapping ' 'non-overlapping'}, 'Location', 'best')
box off
fig = gcf;
shg
end