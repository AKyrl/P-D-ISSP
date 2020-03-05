% GSC.m

run config.m

%% define signal
% run create_micsigs.m;
run load_micsigs_linear_array.m;


%% Estimate DOA using MUSIC_narrrowband
str2num(positions_filename_header(1))
run MUSIC_wideband_linear_array.m;
DOA_est = 90-str2num(positions_filename_header(1))

%% Estimate DAS using DAS_BF
run DAS_BF_linear_array.m;

%% Main program

% generate essential data
Griffiths_Jim_matrix = zeros(size(RIR_sources,2)-1, size(RIR_sources,2));
Griffiths_Jim_matrix(:,1) = 1;
for i=1:size(RIR_sources,2)-1
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
% % d as vector
% for time_index=STFT_L+1:Window_step:size(DAS_out,2)-STFT_L
%     noise_reference_part  = noise_reference(:,time_index-STFT_L+1:time_index);
% %     if sum(VAD(time_index-silence_threshold:time_index))==0 %% update when long silence?
% %     if round(mean(VAD(time_index:time_index+STFT_L-1)))==0 
%     if VAD(time_index)==0 
%         pre_denoise_result =  adaptive_filter_weight.*noise_reference_part;
%         error{error_index} = DAS_out(time_index:time_index+STFT_L-1) - sum(pre_denoise_result);
%         weight_update_coeff = noise_reference_part./(norm(noise_reference_part,'fro').^2).*LMS_step_size.*error{error_index};
%         adaptive_filter_weight = adaptive_filter_weight+weight_update_coeff;
%         error_index = error_index+1;
%     end
%     denoise_result =  adaptive_filter_weight.*noise_reference_part;
%     GSC_out(time_index:time_index+STFT_L-1) = DAS_out(time_index:time_index+STFT_L-1)-sum(denoise_result);
%     GSC_removed(time_index:time_index+STFT_L-1) = GSC_removed(time_index:time_index+STFT_L-1) + sum(denoise_result);
% end

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
    GSC_out(time_index) = DAS_out(time_index)-sum(pre_denoise_result,'all');
    GSC_removed(time_index) = GSC_removed(time_index) + sum(pre_denoise_result,'all');
end


% adaptive filter from dsp project...

%         % denoise-1
%         modulated_signal_fft = fft(modulated_signal(data_package_idx+trainblock_package_per_group_N+package_idx_offset_trainblock+package_idx_offset_data, OFDM_prefix_N+1:end));
%         % LMS
%         draft_denoise_signal                = conj(last_lms_weight) .* modulated_signal_fft;
%         draft_denoise_signal                = reshape(draft_denoise_signal,[length(draft_denoise_signal) 1]);
%         predict_binary_data                 = qam_demod(draft_denoise_signal(on_off_bit_range), QAM_M);
%         ideal_OFDM_frame                    = zeros(1,length(draft_denoise_signal));
%         ideal_OFDM_frame(on_off_bit_range)  = qam_mod(predict_binary_data(1:size(on_off_bit_range)*log2(QAM_M)), QAM_M);
%         ideal_OFDM_frame(-on_off_bit_range+OFDM_frame_N+2) ...
%                                             = conj(ideal_OFDM_frame(on_off_bit_range));
%                                         
%         if max(modulated_signal_fft)~=0 || min(modulated_signal_fft)~=0
%             lms_update_vector = ...
%                 LMS_step_size ./ (LMS_precision_error + conj(modulated_signal_fft) .* modulated_signal_fft).* ...
%                 modulated_signal_fft .* conj(ideal_OFDM_frame - conj(last_lms_weight) .* modulated_signal_fft);
%         else
%             lms_update_vector = zeros(size(last_lms_weight));
%         end
%         
%         channel_est(:,data_package_idx+package_idx_offset_data) = ...
%             last_lms_weight + lms_update_vector;
%             
%         last_lms_weight = ...
%             reshape(channel_est(:,data_package_idx+package_idx_offset_data), size(last_lms_weight));
%         
%                 % denoise-2
%         denoised_modulated_signal_fft = conj(last_lms_weight) .* modulated_signal_fft;

%% visualize
% figure
% hold on;
% plot(reshape(Mic(1,1,:),size(speech_DAS)));
% plot(speech_DAS);
% plot(GSC_out);
% legend('received','DAS','GSC')
% VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
figure
p(1)=subplot(3,1,1);
plot(reshape(Mic(1,1,:),size(speech_DAS)));
title('original');
p(2)=subplot(3,1,2);
plot(speech_DAS);
title('DAS');
p(3)=subplot(3,1,3);
plot(GSC_out);
title('GSC');
% linkaxes(p,'xy')
% ylim([-0.01,0.01])

for i=1:4
    SNR(i) = SNR_cal(squeeze(Mic(1,i,:)));
end
SNR_DAS = SNR_cal(DAS_out);
SNR_GSC = SNR_cal(GSC_out);


%% SNR
function [VAD] = VAD_cal(sig)
    % VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
    VAD=abs(sig)>std(sig)*1e-2;
    mean_step = 1000;
    for i=1:mean_step:length(VAD)-mean_step-1
        VAD(i:i+mean_step-1) = mean(VAD(i:i+mean_step-1))>0.95;
    end
end

function [SNR] = SNR_cal(sig)
    VAD = VAD_cal(sig);
    noise_power = var(sig(VAD==0));
    speech_power = var(sig(VAD==1))-noise_power;
    SNR = 10*log10(speech_power./noise_power);
    % soundsc(GSC_out,fs_RIR)
end