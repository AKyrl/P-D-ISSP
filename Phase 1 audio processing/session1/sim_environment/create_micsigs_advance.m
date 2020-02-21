% create_micsigs.m
% clc; clear all; close all
%% config
length_limit = 10; % (s)
upsampling_FS = 441000;

%% config of impulse response
impulse_positions = ["s30";"s60";"s90";"s-30";"s-60";"s-90"];
impulse_path = '../Head_mounted_IRs/Head_mounted_IRs/';
L1_impulse_filename = '/HMIR_L1.wav';
L2_impulse_filename = '/HMIR_L2.wav';
R1_impulse_filename = '/HMIR_R1.wav';
R2_impulse_filename = '/HMIR_R2.wav';

%% define filenames of waveforms and load waveforms

% speech_filename{1}      = 'audio_files/speech1.wav';
% speech_filename{2}      = 'audio_files/speech2.wav';
% 
% whitenoise_filename{1}  = 'audio_files/whitenoise_signal_1.wav';
% whitenoise_filename{2}  = 'audio_files/whitenoise_signal_2.wav';

dry_filename{1}         = 'audio_files/part1_track1_dry.wav';
% dry_filename{2}         = 'audio_files/part1_track2_dry.wav';
% 
% bubble_filename         = 'audio_files/Babble_noise1.wav';

% noise_filename = dry_filename;
% noise_filename = whitenoise_filename;
% 
% [source{1}, source_Fs{1}]   = audioread(speech_filename{1});
% [source{2}, source_Fs{2}]   = audioread(speech_filename{2});
% 
% 
% [noise{1}, noise_Fs{1}]   = audioread(noise_filename{1});
% [noise{2}, noise_Fs{2}]   = audioread(noise_filename{2});

[target_signal, target_Fs]   = audioread(dry_filename{1});

%% load impulse response
for i=1:length(impulse_positions)
    [impulse_L1{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), L1_impulse_filename));
    impulse_L1{i}  = resample(impulse_L1{i}, upsampling_FS, Fs);
    [impulse_L2{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), L2_impulse_filename));
    impulse_L2{i}  = resample(impulse_L2{i}, upsampling_FS, Fs);
    [impulse_R1{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), R1_impulse_filename));
    impulse_R1{i}  = resample(impulse_R1{i}, upsampling_FS, Fs);
    [impulse_R2{i}, Fs]   = audioread(strcat(impulse_path, impulse_positions(i), R2_impulse_filename));
    impulse_R2{i}  = resample(impulse_R2{i}, upsampling_FS, Fs);
end


%% preprocess for waveform
target_signal  = resample(target_signal, upsampling_FS, target_Fs);
target_signal  = target_signal(1:length_limit* upsampling_FS);

% original_source{1}  = resample(source{1}, fs_RIR, source_Fs{1});
% original_source{2}  = resample(source{2}, fs_RIR, source_Fs{2});
% 
% original_noise{1}   = resample(noise{1}, fs_RIR, noise_Fs{1});
% original_noise{2}   = resample(noise{2}, fs_RIR, noise_Fs{2});
% 
% original_source{1}  = original_source{1}(1:length_limit);
% original_source{2}  = original_source{2}(1:length_limit);
% 
% original_noise{1}  = original_noise{1}(1:length_limit);
% original_noise{2}  = original_noise{2}(1:length_limit);
% 
% noisy_source{1}     = original_source{1} + original_noise{1};
% noisy_source{2}     = original_source{2} + original_noise{2};

%% filter data

for i=1:length(impulse_positions)
    Mic{1,i} = fftfilt(impulse_L1{i}, target_signal);
    Mic{2,i} = fftfilt(impulse_L2{i}, target_signal);
    Mic{3,i} = fftfilt(impulse_R1{i}, target_signal);
    Mic{4,i} = fftfilt(impulse_R2{i}, target_signal);
end

%% combine signal from each speaker

for i=1:4
    sound{i} = Mic{i,1};
end

for i=2:length(impulse_positions)
    for j=1:4
        sound{j} = sound{j} + Mic{j,i};
    end
end


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
