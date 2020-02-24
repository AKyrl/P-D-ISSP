% create_micsigs.m
clear all; 
%% config
length_limit = 10; % (s)

%% define filenames of waveforms

speech_filename{1}      = 'audio_files/speech1.wav';
speech_filename{2}      = 'audio_files/speech2.wav';

whitenoise_filename{1}  = 'audio_files/whitenoise_signal_1.wav';
whitenoise_filename{2}  = 'audio_files/whitenoise_signal_2.wav';

dry_filename{1}         = 'audio_files/part1_track1_dry.wav';
dry_filename{2}         = 'audio_files/part1_track2_dry.wav';

bubble_filename{1}         = 'audio_files/Babble_noise1.wav';

% noise_filename = dry_filename;
noise_filename = bubble_filename;

%% load waveforms

[source{1}, source_Fs{1}]   = audioread(speech_filename{1});
[source{2}, source_Fs{2}]   = audioread(speech_filename{2});


[noise{1}, noise_Fs{1}]   = audioread(noise_filename{1});
% [noise{2}, noise_Fs{2}]   = audioread(noise_filename{2});

%% load RIRs
load Computed_RIRs.mat

%% preprocess for waveform (resample, truncation)
original_source{1}  = resample(source{1}, fs_RIR, source_Fs{1});
original_source{2}  = resample(source{2}, fs_RIR, source_Fs{2});

original_noise{1}   = resample(noise{1}, fs_RIR, noise_Fs{1});
% original_noise{2}   = resample(noise{2}, fs_RIR, noise_Fs{2});

VAD=abs(original_source{1}(:,1))>std(original_source{1}(:,1))*1e-3;
original_source{1} =  original_source{1}(VAD==1,1);
original_source{1}  = original_source{1}(1:length_limit*fs_RIR);

VAD=abs(original_source{2}(:,1))>std(original_source{2}(:,1))*1e-3;
original_source{2} =  original_source{2}(VAD==1,1);
original_source{2}  = original_source{2}(1:length_limit*fs_RIR);
%soundsc(original_source{1},fs_RIR)

original_noise{1}  = original_noise{1}(1:length_limit*fs_RIR);
% original_noise{2}  = original_noise{2}(1:length_limit*fs_RIR);

noisy_source{1}     = original_source{1} + original_noise{1};
% noisy_source{2}     = original_source{2} + original_noise{2};

%% filter data
Mic_noise = zeros(      length(original_noise), ...
                        size(RIR_sources,2), ...
                        length_limit*fs_RIR);

Mic = zeros(    length(original_source), ...
                size(RIR_sources,2), ...
                length_limit*fs_RIR);

for i=1:length(original_noise)
    for mic_idx=1:size(RIR_noise,2)
        Mic_noise(i,mic_idx,:)= fftfilt(RIR_noise(:,mic_idx,i), original_noise{i});
    end
end

for i=1:length(original_source)
    for mic_idx=1:size(RIR_sources,2)
        Mic(i,mic_idx,:)= fftfilt(RIR_sources(:,mic_idx), original_source{i});
    end
end





%% for week 3 (adding noise and filter out weak signal)
for channel_idx=1:size(Mic(:,:,:),1)
    for mic_idx=1:size(RIR_sources,2)
        microphone_power(channel_idx,mic_idx)  = var(Mic(channel_idx,mic_idx,:),0,'all');
        noise = wgn(1,size(Mic(channel_idx,mic_idx,:),3),10*log10(0.1*microphone_power(channel_idx,1)));
        noise = reshape(noise,[1,1,length(noise)]);
        noise = noise + Mic_noise(1,mic_idx,:);
        noise_power(channel_idx,mic_idx) = var(noise,0,'all');
        speech_noise(channel_idx,mic_idx,:) = Mic(channel_idx,mic_idx,:) + noise;
    end
end

SNR = 10*log10(microphone_power./noise_power);

sound = speech_noise(1,1,:);
sound = reshape(sound,[1,length(sound)]);
% soundsc(sound,fs_RIR)



% for i=1:length(noisy_source)
%     for j=1:size(RIR_sources,2)
%         Mic{j,i} = fftfilt(RIR_sources(:,j,i), noisy_source{i});
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
