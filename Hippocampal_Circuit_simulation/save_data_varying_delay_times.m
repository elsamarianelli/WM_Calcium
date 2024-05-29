function[] = save_data_varying_delay_times(C, J, p, ca3_ensembles, n_trials, variable, main_folder)
% function to create a folder of test and train data for perceptron for a
% variety of delay times
    %% set parameters for the simulation
    p.SimLength             = p.start_time+p.length_first+p.delay_time+p.length_second+200;

    % set extra paramaters for single trial
    input.simulation        = [p.start_time p.start_time+p.length_first];
    input.reactivation      = [p.start_time+p.length_first+p.delay_time p.start_time+p.length_first+p.delay_time+p.length_second];

    % get times to use to train perceptorn 
    time3 = input.reactivation(1);
    time4 = time3+100; 
    
    % %%  Simulate hippocampal dynamics  over many trials labelling with reward/no reward
    [spikeCounts,~]	= get_train_data_db(C, J, n_trials, p, ca3_ensembles, time3, time4);
    [spikeCounts_test, ~]	= get_train_data_db(C, J, 1, p, ca3_ensembles, time3, time4);
    
    %% Save data and settings info to current working directory
    % folderName = ([num2str(p.tau_facil) '_' num2str(p.tau_decay)]);
    folderName = num2str(variable);
    
    % Full path for the new folder
    fullFolderPath = fullfile(main_folder, folderName);
    
    % Check if the directory exists; if not, create it
    if ~exist(fullFolderPath, 'dir')
        mkdir(fullFolderPath);
    end
    
    % Save in .mat file
    save(fullfile(fullFolderPath, 'spikeCounts.mat'), 'spikeCounts');
    save(fullfile(fullFolderPath, 'spikeCounts_test.mat'), 'spikeCounts_test');
    save(fullfile(fullFolderPath, 'myStruct.mat'), 'p');
    save(fullfile(fullFolderPath, 'Con_matrix.mat'), 'C');
    save(fullfile(fullFolderPath, 'Ca3_ensembles.mat'), 'ca3_ensembles');

end