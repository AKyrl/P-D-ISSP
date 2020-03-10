% preprocess_64_eeg.m
% clear all

%% config
eeg_file_folder = '../64_channel_Smarting_EEG_data/';

%% main process
% some static config
Number_of_train_subjects = 12;
Number_of_mat_per_train_subjects = 4;
Number_of_test_subjects = 12;
Number_of_mat_per_test_subjects = 4;



% add third_party
%addpath('../third_party/biosig4c++-2.0.2_mex');
% addpath('../third_party/eeglab2019_1');
% addpath('../third_party/eeglab2019_1/functions/adminfunc');
% addpath('../third_party/eeglab2019_1/functions/popfunc');
% biosigpathfirst  

% load eeg data
for subject_idx = 1:Number_of_train_subjects
    current_dir = dir(strcat(eeg_file_folder,'Subject',num2str(subject_idx),'/'));
    for file_idx = 1:Number_of_gdf_per_person
        EEG{subject_idx, file_idx} = pop_biosig(strcat(eeg_file_folder,'Subject',num2str(subject_idx),'/',current_dir(file_idx+2).name));
    end
end

% eegread
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_gdf_per_person
        AADEEG{subject_idx, data_idx} = eegread(EEG{subject_idx, data_idx}, 33027, 33025);
    end
end

% downsample again
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_gdf_per_person
        AADEEG_downsampled{subject_idx, data_idx} = zeros(ceil(size(AADEEG{subject_idx, data_idx}).*[1,20/EEG{subject_idx, data_idx}.srate]));
        for channel_idx = 1:size(AADEEG{subject_idx, data_idx}, 1)
            AADEEG_downsampled{subject_idx, data_idx}(channel_idx, :) = resample(double(AADEEG{subject_idx, data_idx}(channel_idx, :)), 20, EEG{subject_idx, data_idx}.srate);
        end
    end
end
fs = 20;

% BPF
[b,a] = butter(5,[1 9]/(fs/2),'bandpass');
for subject_idx = 1:Number_of_train_subjects
    for data_idx = 1:Number_of_gdf_per_person
        for channel_idx = 1:size(AADEEG{subject_idx, data_idx}, 1)
            AADEEG_downsampled{subject_idx, data_idx}(channel_idx, :) = filter(b,a,AADEEG_downsampled{subject_idx, data_idx}(channel_idx, :));
        end
    end
end

