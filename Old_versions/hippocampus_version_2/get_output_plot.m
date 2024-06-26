function [fig] = get_output_plot(M,pattern_order, p, mems, C, ca3_ensembles, ca1_ensembles)

stim                    = cell(2,1); 
first                   = double(upper(p.pattern_order(1))) - 64; 
stim{1}                 = ca3_ensembles{first}; %clear first
second                  = double(upper(p.pattern_order(2))) - 64; 
stim{2}                 = ca3_ensembles{second}; %clear second

figure
% plotting vm in both patterns
fs = 10;
ns = 4; 
n = 0;

% plotting synaptic parameters u an x for synapses potentiated by first
% odour which are reactivated during second odour presentation
ax1 = subplot(ns,1,n+1);
av_u_memory = M.U_mem2_log;
av_x_memory = M.X_mem2_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(pattern_order(1),'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('mem_2')
%u and x for non overlapping 2nd odor cells
ax2 = subplot(ns,1,n+2);
av_u_memory = M.U_ovlp_log;
av_x_memory = M.X_ovlp_log;
plot(1:p.SimLength,av_u_memory,'b')
ylabel(pattern_order(2),'FontSize',fs)
hold on
plot(1:p.SimLength,av_x_memory,'r')
legend('u', 'x' ,'Location','southeast')
ylim([0 1])
ylabel('overlapping')

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

L1 = plot(nan, nan, 'color', 'b');
L2 = plot(nan, nan, 'color', 'r');
L3 = plot(nan, nan, 'color', 'm');
legend([L1, L2, L3], {'odour 1' 'odour 2', '10% over non-odour cells'}, 'Location', 'best')

% plotting spike raster for CA1 ouput cell firing
% get CA1 cells which are able to receive input from both A and B
CA3_overlap = intersect(ca3_ensembles{first}, ca3_ensembles{second});
ind = find((sum(C(CA3_overlap, :))>4)); 
CA1_overlap_second =[(intersect(ca1_ensembles{second}, ind+p.in))];
CA1_second = setdiff(ca1_ensembles{second}, CA1_overlap_second);

% overlapping with more than 2 inputs
ax7 = subplot(ns, 1, n+4);
spikeMat = M.spikelog(CA1_overlap_second, :);
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
% others (non overlapping )
spikeMat = M.spikelog(CA1_second, :);
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


L1 = plot(nan, nan, 'color', 'm');
L2 = plot(nan, nan, 'color', 'b');
legend([L1, L2], {'>x overlapping odour 2 ' 'other odour 2'}, 'Location', 'southwest')

fig = gcf;
shg
end