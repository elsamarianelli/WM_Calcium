function[] = save_data_varying_delay_times(C, J, p, ca3_ensembles, delay_time, n_trials, main_folder)
% function to create a folder of test and train data for perceptron for a
% variety of delay times
    %% Set parameters for the simulation
    p.delay_time            = delay_time;            % Delay between odour presentations (ms)
    p.SimLength             = p.start_time+p.length_first+p.delay_time+p.length_second+200;

    % set extra paramaters for single trial
    input.simulation        = [p.start_time p.start_time+p.length_first];
    input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];

    % get times to use to train perceptorn 
    time3 = input.reactivation(1);
    time4 = time3+100; 
    
    % %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
    [spikeCounts,~]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
    [spikeCounts_test, ~]	= get_train_data_db(C, J, 10, p, ca3_ensembles, time3, time4);
    
    %% Save data and settings info to current working directory
    folderName = num2str(delay_times);
    
    % Create the folder if it doesn't already exist
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    % Save the matrix in a .mat file
    save(fullfile(main_folder, folderName, 'spikeCounts.mat'), 'spikeCounts');
    save(fullfile(main_folder, folderName, 'spikeCounts_test.mat'), 'spikeCounts_test');
    % Save the structure in a .mat file
    save(fullfile(main_folder, folderName, 'myStruct.mat'), 'p');

end