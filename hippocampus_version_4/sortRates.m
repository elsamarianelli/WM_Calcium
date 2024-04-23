function[out]   = sortRates(p,M,input,ca3_ensembles,ca1_ensembles)
%% Script to compute mean firing rates of each sub-population
out.labels      = {'Background','A','B','C'};
ca3ids          = 1 : p.in;
ca3mem          = horzcat(ca3_ensembles{:});
ca3ids(ca3mem)  = []; clear ca3mem
ca1ids          = p.in + (1 : p.out);
ca1mem          = horzcat(ca1_ensembles{:});
ca1ids(ca1mem-p.in)  = []; clear ca1mem

ca1rates        = cellfun(@(x) M.spikelog(x,:),ca1_ensembles,'UniformOutput',false);
ca3rates        = cellfun(@(x) M.spikelog(x,:),ca3_ensembles,'UniformOutput',false);

out.ca3Rates	= cellfun(@(x) mean(x,1),ca3rates,'UniformOutput',false);
out.ca3stim1    = mean(M.spikelog(:,input.simulation(1):input.simulation(2)),2) ./ (diff(input.simulation)/1000);
out.ca3stim2    = mean(M.spikelog(:,input.reactivation(1):input.reactivation(2)),2) ./ (diff(input.reactivation)/1000);
out.ca3Rates    = [mean(M.spikelog(ca3ids,:),1) ; vertcat(out.ca3Rates{:})]; clear ca3ids ca3rates
out.ca3stim1Mn  = sum(out.ca3Rates(:,input.simulation(1):input.simulation(2)),2) ./ (diff(input.simulation)/1000);
out.ca3delayMn  = sum(out.ca3Rates(:,input.simulation(2)+1:input.reactivation(1)-1),2) ./ ((input.reactivation(1)-input.simulation(2)-2)/1000);
out.ca3stim2Mn  = sum(out.ca3Rates(:,input.reactivation(1):input.reactivation(2)),2) ./ (diff(input.reactivation)/1000);

out.ca1Rates    = cellfun(@(x) mean(x,1),ca1rates,'UniformOutput',false);
out.ca1stim1    = mean(M.spikelog(:,input.simulation(1):input.simulation(2)),2) ./ (diff(input.simulation)/1000);
out.ca1stim2    = mean(M.spikelog(:,input.reactivation(1):input.reactivation(2)),2)  ./ (diff(input.reactivation)/1000);
out.ca1Rates	= [mean(M.spikelog(ca1ids,:),1) ; vertcat(out.ca1Rates{:})]; clear ca1ids ca1rates
out.ca1stim1Mn  = sum(out.ca1Rates(:,input.simulation(1):input.simulation(2)),2) ./ (diff(input.simulation)/1000);
out.ca1delayMn  = sum(out.ca3Rates(:,input.simulation(2)+1:input.reactivation(1)-1),2) ./ ((input.reactivation(1)-input.simulation(2)-2)/1000);
out.ca1stim2Mn	= sum(out.ca1Rates(:,input.reactivation(1):input.reactivation(2)),2) ./ (diff(input.reactivation)/1000);