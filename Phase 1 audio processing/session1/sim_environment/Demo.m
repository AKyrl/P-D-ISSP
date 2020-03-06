% Demo.m
clear all

%% config
length_limit    = 200;
change_time     = 5;
mic_distance    = 5;
DOA_list        = [-90; -90; -60; -60; -90];
STFT_L          = 1024;
STFT_overlap    = 50;
number_of_source_channel = 1;
LMS_step_size   = 0.1;

%% main process

[mic_signals fs]    = load_audio_files(length_limit);
[received_signal]   = signal_fusion(mic_signals, change_time, fs);

DAS_out     = zeros(1,size(received_signal,2));
speech      = zeros(size(received_signal,1),size(received_signal,2));
GSC_out     = zeros(1,size(received_signal,2));
GSC_removed = zeros(1,size(received_signal,2));
DOA_overtime= zeros(number_of_source_channel,size(received_signal,2));

DOA_idx = 1;

for time_index_start=1:change_time*fs:size(received_signal,2)
    time_index_end = time_index_start+fs*change_time-1;
    
    [DOA_est] = DOA_estimation( received_signal(:,time_index_start:time_index_end), ...
        mic_distance, STFT_L, STFT_overlap, fs, number_of_source_channel);
    
    DOA_est_sorted = sort(DOA_est,'descend');
    if length(DOA_est_sorted)<number_of_source_channel
        DOA_est_sorted(length(DOA_est_sorted)+1:number_of_source_channel)=0;
    end    M1_filename_header          = 'HML_M1';
    M2_filename_header          = 'HML_M2';
    M3_filename_header          = 'HMR_M1';
    M4_filename_header          = 'HMR_M2';

    DOA_overtime(:,time_index_start:time_index_end) = DOA_est_sorted'*ones(1,time_index_end-time_index_start+1);
    
%     DOA_est = DOA_list(DOA_idx);
%     DOA_idx = DOA_idx+1;
%     if DOA_idx>5
%         DOA_idx = 1;
%     end
   
    
    [DAS_out(time_index_start:time_index_end) speech(:,time_index_start:time_index_end) ]= ...
        DAS_BF(received_signal(:,time_index_start:time_index_end), mic_distance, DOA_est(1), fs);
    [GSC_out(:,time_index_start:time_index_end) GSC_removed(:,time_index_start:time_index_end) ]= ...
        GSC(DAS_out(time_index_start:time_index_end), ...
        speech(:,time_index_start:time_index_end), ...
        STFT_L, LMS_step_size);
end

%% visualize
figure
p(1)=subplot(3,1,1);
plot(squeeze(mic_signals(1,1,:)));
title('original');
p(2)=subplot(3,1,2);
plot(DAS_out);
title('DAS');
p(3)=subplot(3,1,3);
plot(GSC_out);
title('GSC');

for i=1:4
    SNR(i) = SNR_cal(squeeze(mic_signals(1,i,:)));
end
SNR_DAS = SNR_cal(DAS_out);
SNR_GSC = SNR_cal(GSC_out);

pause;


%% load signals
function [mic_signals fs] = load_audio_files(length_limit)

    impulse_path                = '../Head_mounted_real/';
    positions_filename_header   = ["-90_part2_track1";"90_part2_track2";"-60_part2_track1";"60_part2_track2"];
    M1_filename_header          = 'LMA_M1';
    M2_filename_header          = 'LMA_M2';
    M3_filename_header          = 'LMA_M3';
    M4_filename_header          = 'LMA_M4';
%     M1_filename_header          = 'HML_M1';
%     M2_filename_header          = 'HML_M2';
%     M3_filename_header          = 'HMR_M1';
%     M4_filename_header          = 'HMR_M2';

    % load fs
    [audio, fs]  = audioread(strcat(impulse_path, M1_filename_header,'_',positions_filename_header(1),'.wav'));

    mic_signals = zeros(    length(positions_filename_header), ...
                    4, ...
                    length_limit*fs);

    for i=1:length(positions_filename_header)
        [audio, fs]  = audioread(strcat(impulse_path, M1_filename_header,'_',positions_filename_header(i),'.wav'));
        mic_signals(i,1,:) = audio(1:length_limit*fs);
        [audio, fs]  = audioread(strcat(impulse_path, M2_filename_header,'_',positions_filename_header(i),'.wav'));
        mic_signals(i,2,:) = audio(1:length_limit*fs);
        [audio, fs]  = audioread(strcat(impulse_path, M3_filename_header,'_',positions_filename_header(i),'.wav'));
        mic_signals(i,3,:) = audio(1:length_limit*fs);
        [audio, fs]  = audioread(strcat(impulse_path, M4_filename_header,'_',positions_filename_header(i),'.wav'));
        mic_signals(i,4,:) = audio(1:length_limit*fs);
    end
end

%% generate signals
function [received_signal] = signal_fusion(mic_signals, change_time, fs)
    received_signal = zeros(4, size(mic_signals,3));
    switch_list = [1,2; 1,4; 3,4; 3,2; 1,4];
    time_index  = 1;
    switch_list_idx = 1;
    
    while((time_index+fs*change_time-1)<=size(mic_signals,3))
        received_signal(:,time_index:time_index+fs*change_time-1) = ...
            squeeze(mic_signals(switch_list(switch_list_idx,1),:,time_index:time_index+fs*change_time-1)) + ...
            squeeze(mic_signals(switch_list(switch_list_idx,2),:,time_index:time_index+fs*change_time-1));
        time_index = time_index+fs*change_time;
        switch_list_idx = switch_list_idx+1;
        if(switch_list_idx>size(switch_list, 1))
            switch_list_idx = 1;
        end
    end
end

%% DOA estimation
function [DOA_est] = DOA_estimation(target_signals, mic_distance_, STFT_L, STFT_overlap, fs, number_of_source_channel)

    %% STFT
    target_signals_stft = zeros( size(target_signals,1), STFT_L, ...
        round(size(target_signals,2)/(round((100-STFT_overlap)*STFT_L/100)))-2);
    for i=1:size(target_signals,1)
        target_signals_stft(i,:,:) = stft(target_signals(i,:), ...
            'Window',hanning(STFT_L),...
            'OverlapLength',round(STFT_overlap*STFT_L/100));
    end


    %% evalute the pseudospectrum of each frequency bins
    angles                              = [0:0.5:180];
    Music_pseudospectrum_of_each_bins   = zeros(STFT_L/2-1, length(angles));

    target_signals_stft_power = zeros(size(target_signals_stft,1),STFT_L/2-1);

    for frequency_bin_idx=1+STFT_L/2:STFT_L-1

        w = (frequency_bin_idx-STFT_L/2)*(fs/STFT_L);
        recorded_signal_at_w  = ...
            zeros(size(target_signals_stft,1), size(target_signals_stft,3));
        for i=1:size(target_signals_stft,1)
            recorded_signal_at_w(i,:) = ...
                target_signals_stft(i,frequency_bin_idx,:);
            target_signals_stft_power(i,frequency_bin_idx-STFT_L/2) = ...
                mean(abs(recorded_signal_at_w(i,:)).^2);
        end

        correlation_matrix = recorded_signal_at_w*recorded_signal_at_w';
        [EigenVector EigenValue ] = eig(correlation_matrix);
        EigenVector = EigenVector(:,1:size(EigenVector,2)-number_of_source_channel);

        Mic_position = [0:mic_distance_:mic_distance_*(size(target_signals_stft,1)-1)];

        for angle_idx=1:length(angles)
            manifold_vector = exp((cosd(angles(angle_idx))*Mic_position./340./100).*1i.*w'*2*pi);
            manifold_vector = reshape(manifold_vector, [size(target_signals_stft,1) 1]);
            Music_pseudospectrum_of_each_bins(frequency_bin_idx-STFT_L/2, angle_idx) = ...
                1/(manifold_vector'* EigenVector* EigenVector'* manifold_vector);
        end

    end

    % Music_pseudospectrum = geomean(abs(Music_pseudospectrum_of_each_bins));
    % Music_pseudospectrum = mean(abs(Music_pseudospectrum_of_each_bins));
    % weighted mean
    target_signals_stft_power = sum(target_signals_stft_power);
    target_signals_stft_power(find(target_signals_stft_power==0)) = 1e-10;
    Music_pseudospectrum = sum(transpose(target_signals_stft_power).*abs(Music_pseudospectrum_of_each_bins()))./sum(target_signals_stft_power);

    Music_pseudospectrum_sorted = findpeaks(Music_pseudospectrum);
    if sum(find(Music_pseudospectrum_sorted==max(Music_pseudospectrum)))==0
        Music_pseudospectrum_sorted = [max(Music_pseudospectrum) Music_pseudospectrum_sorted];
    end
    Music_pseudospectrum_sorted = sort(Music_pseudospectrum_sorted,'descend');

    for DOA_idx=1:min(number_of_source_channel,length(Music_pseudospectrum_sorted))
        DOA_est(DOA_idx) = angles(find(Music_pseudospectrum==Music_pseudospectrum_sorted(DOA_idx),1));
    end
    
end

%% DAS
function [DAS_out speech] = DAS_BF(received_signal, mic_distance_, DOA_est, fs)
    delay = abs(cosd(DOA_est)*mic_distance_/100/340);

    samples = ceil(delay*fs);

    speech = zeros(size(received_signal,1),size(received_signal,2));
    speech_DAS = zeros(1,size(received_signal,2));

    if(DOA_est<90)
        for mic_idx=1:size(received_signal,1)
            speech_temp = received_signal(mic_idx,:);
            speech_temp = reshape(speech_temp,size(speech_DAS));
            speech(mic_idx, 1:end-(mic_idx-1)*samples) = speech_temp((mic_idx-1)*samples+1:end);
            speech_DAS(1:end-(mic_idx-1)*samples) = speech(mic_idx,1:end-(mic_idx-1)*samples)+ speech_DAS(1:end-(mic_idx-1)*samples);
        end
    elseif(DOA_est>=90)
         for mic_idx=1:size(received_signal,1)
            count = size(received_signal,1)- mic_idx+1;
            speech_temp = received_signal(mic_85idx,:);;
            speech_temp = reshape(speech_temp,size(speech_DAS));
            speech(mic_idx, 1:end-(count-1)*samples) = speech_temp((count-1)*samples+1:end);
            speech_DAS(1:end-(count-1)*samples) = speech(mic_idx,1:end-(count-1)*samples)+ speech_DAS(1:end-(count-1)*samples);
        end
    else
        print 'you are doomed';
    end

    speech_DAS = speech_DAS/size(received_signal,1);

    DAS_out = 2*speech_DAS;

end

%% GSC
function [GSC_out GSC_removed] = GSC(DAS_out, speech, STFT_L, LMS_step_size)


    % generate essential data
    Griffiths_Jim_matrix = zeros(size(speech,1)-1, size(speech,1));
    Griffiths_Jim_matrix(:,1) = 1;
    for i=1:size(speech,1)-1
        Griffiths_Jim_matrix(i,i+1)=-1;
    end

    % get noise reference
    noise_reference = Griffiths_Jim_matrix*speech;

    % get VHD reference
    % VAD=abs(DAS_out)>std(DAS_out)*1e-3;
    % VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
    VAD = VAD_cal(DAS_out);

    % adaptive filter over the time
    adaptive_filter_weight = ones(size(noise_reference,1),STFT_L)./STFT_L;
    % adaptive_filter_weight = zeros(size(noise_reference,1),STFT_L);

    % reserve variable for denoise result
    GSC_out = DAS_out;
    GSC_removed = zeros(size(DAS_out));

    silence_threshold = 0;
    Window_step = 1;
    error_index = 1;
   
    % d as value
    for time_index=STFT_L+1:size(DAS_out,2)-STFT_L
        noise_reference_part  = noise_reference(:,time_index-STFT_L+1:time_index);
        pre_denoise_result =  adaptive_filter_weight.*noise_reference_part;
        if sum(VAD(time_index-silence_threshold:time_index))==0 %% update when long silence?
            error{error_index} = DAS_out(time_index) - sum(pre_denoise_result,'all');
            weight_update_coeff = noise_reference_part./(norm(noise_reference_part,'fro').^2).*LMS_step_size.*error{error_index};
            adaptive_filter_weight = adaptive_filter_weight+weight_update_coeff;
            error_index = error_index+1;
        end
        denoise_result =  adaptive_filter_weight.*noise_reference_part;
        GSC_out(time_index) = DAS_out(time_index)-sum(denoise_result,'all');
        GSC_removed(time_index) = GSC_removed(time_index) + sum(denoise_result,'all');
    end

end

%% SNR
function [VAD] = VAD_cal(sig, detect_threshold, mean_step, mean_threshold)

    if nargin < 2
        detect_threshold = 1.5*1e-2;
    end
    if nargin < 3
        mean_step = 1000;
    end
    if nargin < 4
        mean_threshold = 0.97;
    end
    % VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
    VAD=abs(sig)>std(sig)*detect_threshold;
    
    for i=1:mean_step:length(VAD)-mean_step-1
        VAD(i:i+mean_step-1) = mean(VAD(i:i+mean_step-1))>mean_threshold;
    end
end

function [SNR] = SNR_cal(sig, VAD)
    sig(isnan(sig))=0;
    
    if nargin < 2
        VAD = VAD_cal(sig);
    end
    
    noise_power = var(sig(VAD==0));
    speech_power = var(sig(VAD==1))-noise_power;
    SNR = 10*log10(speech_power./noise_power);
end