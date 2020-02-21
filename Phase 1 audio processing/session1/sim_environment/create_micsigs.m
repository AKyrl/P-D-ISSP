% create_micsigs.m
% clc; clear all; close all
%% config
length_limit = 10; % (s)

%% define filenames of waveforms

speech_filename{1}      = 'audio_files/speech1.wav';
speech_filename{2}      = 'audio_files/speech2.wav';

whitenoise_filename{1}  = 'audio_files/whitenoise_signal_1.wav';
whitenoise_filename{2}  = 'audio_files/whitenoise_signal_2.wav';

dry_filename{1}         = 'audio_files/part1_track1_dry.wav';
dry_filename{2}         = 'audio_files/part1_track2_dry.wav';

bubble_filename         = 'audio_files/Babble_noise1.wav';

% noise_filename = dry_filename;
noise_filename = whitenoise_filename;

%% load waveforms

[source{1}, source_Fs{1}]   = audioread(speech_filename{1});
[source{2}, source_Fs{2}]   = audioread(speech_filename{2});


[noise{1}, noise_Fs{1}]   = audioread(noise_filename{1});
[noise{2}, noise_Fs{2}]   = audioread(noise_filename{2});

%% load RIRs
load Computed_RIRs.mat

%% preprocess for waveform (resample, truncation)
original_source{1}  = resample(source{1}, fs_RIR, source_Fs{1});
% original_source{2}  = resample(source{2}, fs_RIR, source_Fs{2});

original_noise{1}   = resample(noise{1}, fs_RIR, noise_Fs{1});
original_noise{2}   = resample(noise{2}, fs_RIR, noise_Fs{2});

original_source{1}  = original_source{1}(1:length_limit*fs_RIR);
% original_source{2}  = original_source{2}(1:length_limit*fs_RIR);

original_noise{1}  = original_noise{1}(1:length_limit*fs_RIR);
original_noise{2}  = original_noise{2}(1:length_limit*fs_RIR);

noisy_source{1}     = original_source{1} + original_noise{1};
% noisy_source{2}     = original_source{2} + original_noise{2};

%% filter data
Mic = zeros(    length(original_source)+length(noisy_source), ...
                size(RIR_sources,2), ...
                length_limit*fs_RIR);

for i=1:length(original_source)
    for j=1:size(RIR_sources,2)
        Mic(i,j,:)= fftfilt(RIR_sources(:,j,i), original_source{i});
    end
end
% for i=1:length(noisy_source)
%     for j=1:size(RIR_sources,2)
%         Mic{j,i} = fftfilt(RIR_sources(:,j,i), noisy_source{i});
%     end
% end


% for i=1:length(original_noise)
%     for j=1:size(RIR_noise,2)
%         Mic{j,i+length(original_source)} = fftfilt(RIR_noise(:,j,i), original_noise{i});
%     end
% end

%% visualize seperate result

% for j=1:size(Mic, 1)
%     figure('name', ['Mic',num2str(j)])
%     hold on
%     for i=1:size(Mic, 2)
%         plot(Mic{j,i})
%     end
% end

%% result fustion and visualization

% figure('name', 'Mic_1')
% hold on
% for j=1:size(Mic, 1)
%     Mic_1{j} = Mic{j,1};
%     for i=2:size(Mic, 2)
%         Mic_1{j} = Mic{j,i} + Mic_1{j};
%     end
%     plot(Mic_1{j})
% end

%% play result with bubble_noise
% Mic_2 = fftfilt(RIR_sources(:,1,1), original_source{1}(1:length_limit));
% Mic_2 = Mic_2 + fftfilt(RIR_noise(:,1,1), bubble(1:length_limit));
% soundsc(Mic_2,fs_RIR)


%% old code
% %% mySA - preprocess
% 
% %% mySA - RIRs estimation
% xy_mic = zeros(S.L_mic,2);
% xy_audio = zeros(S.L_audio,2);
% 
% 
% for k = 1:S.nmic
%     xy_mic (k,:) = [get(S.hpc_mic(k),'Xdata') get(S.hpc_mic(k),'Ydata')];
% end
% 
% for k = 1:S.L_audio
%     xy_audio (k,:) = [get(S.hpc_audio(k),'Xdata') get(S.hpc_audio(k),'Ydata')];
% end
% 
% if S.L_noise > 0
%     xy_noise = zeros(S.L_noise,2);
%     for k = 1:S.L_noise
%         xy_noise (k,:) = [get(S.hpc_noise(k),'Xdata') get(S.hpc_noise(k),'Ydata')];
%     end
% else
%     xy_noise = [];
% end
% 
% %% mySA - final process (save)
% create_rirs(xy_mic,xy_audio,xy_noise,S.rdim*[1 1],S.reverb,S.fs,S.lRIR);
% ed = msgbox('RIRs created and stored in Computed_RIRs.mat!');
% uiwait(ed);
