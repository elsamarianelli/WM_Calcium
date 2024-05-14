function [data, data_test] = loadData(folderName)
    load(fullfile(folderName, 'spikeCounts.mat'));  % Loads 'spikeCounts'
    load(fullfile(folderName, 'spikeCounts_test.mat'));  % Loads 'spikeCounts_test'
    data = spikeCounts;
    data_test = spikeCounts_test; 
end