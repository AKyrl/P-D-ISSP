% GSC.m

run config.m

%% define signal
% run create_micsigs.m;
positions_filename_header   = ["-60"];
type_filename_header = "part2_track1";
run load_micsigs_linear_array.m;
Mic2=Mic;
positions_filename_header   = ["60"]; 
type_filename_header = "part2_track2";
run load_micsigs_linear_array.m;
Mic1=Mic;
Mic=Mic1+Mic2;



%% Estimate DOA using MUSIC_narrrowband
%str2num(positions_filename_header(1))
run MUSIC_wideband_linear_array.m;
DOA_est = str2num(positions_filename_header(1))

%% Estimate DAS using DAS_BF
run DAS_BF_linear_array.m;

%% Main program
VAD_ideal = VAD_cal_ideal(squeeze(Mic1(1,1,:)));

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
    GSC_out(time_index) = DAS_out(time_index)-sum(pre_denoise_result,'all');
    GSC_removed(time_index) = GSC_removed(time_index) + sum(pre_denoise_result,'all');
end

%% run again with VAD ideal
% adaptive filter over the time
adaptive_filter_weight = ones(size(noise_reference,1),STFT_L)./STFT_L;
% adaptive_filter_weight = zeros(size(noise_reference,1),STFT_L);

% reserve variable for denoise result
GSC_out_ideal = DAS_out;
GSC_removed_ideal = zeros(size(DAS_out));

silence_threshold = 0;
Window_step = 1;
error_index = 1;
% d as value
for time_index=STFT_L+1:size(DAS_out,2)-STFT_L
    noise_reference_part  = noise_reference(:,time_index-STFT_L+1:time_index);
    pre_denoise_result =  adaptive_filter_weight.*noise_reference_part;
    if sum(VAD_ideal(time_index-silence_threshold:time_index))==0 %% update when long silence?
        error_ideal{error_index} = DAS_out(time_index) - sum(pre_denoise_result,'all');
        weight_update_coeff = noise_reference_part./(norm(noise_reference_part,'fro').^2).*LMS_step_size.*error{error_index};
        adaptive_filter_weight = adaptive_filter_weight+weight_update_coeff;
        error_index = error_index+1;
    end
    denoise_result =  adaptive_filter_weight.*noise_reference_part;
    GSC_out_ideal(time_index) = DAS_out(time_index)-sum(denoise_result,'all');
    GSC_removed_ideal(time_index) = GSC_removed(time_index) + sum(denoise_result,'all');
end


%% visualize
% figure
% hold on;
% plot(reshape(Mic(1,1,:),size(speech_DAS)));
% plot(speech_DAS);
% plot(GSC_out);
% legend('received','DAS','GSC')
% VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
figure
p(1)=subplot(4,1,1);
plot(reshape(Mic(1,1,:),size(speech_DAS)));
title('original');
p(2)=subplot(4,1,2);
plot(speech_DAS);
title('DAS');
p(3)=subplot(4,1,3);
plot(GSC_out);
title('GSC');
p(4)=subplot(4,1,4);
plot(GSC_out_ideal);
title('GSC_ideal');
% linkaxes(p,'xy')
% ylim([-0.01,0.01])

for i=1:4
    SNR(i) = SNR_cal(squeeze(Mic(1,i,:)));
end
SNR_DAS = SNR_cal(DAS_out);
SNR_GSC = SNR_cal(GSC_out);
SNR_GSCi = SNR_cal(GSC_out_ideal);

for i=1:4
    SNR_ideal(i) = SNR_cal(squeeze(Mic(1,i,:)),VAD_ideal);
end
SNR_DAS_ideal = SNR_cal(DAS_out,VAD_ideal);
SNR_GSC_ideal = SNR_cal(GSC_out,VAD_ideal);
SNR_GSCi_ideal = SNR_cal(GSC_out_ideal,VAD_ideal);

%% SNR
function [VAD] = VAD_cal_ideal(sig)
    VAD=abs(sig)>std(sig)*1e-3;
    mean_step = 20;
    for i=1:mean_step:length(VAD)-mean_step-1
        VAD(i:i+mean_step-1) = mean(VAD(i:i+mean_step-1))>0.85;
    end
end

function [VAD] = VAD_cal(sig)
    % VAD=abs(squeeze(Mic(1,1,:)))>std(squeeze(Mic(1,1,:)))*1e-3;
    VAD=abs(sig)>std(sig)*1e-2;
    mean_step = 1000;
    for i=1:mean_step:length(VAD)-mean_step-1
        VAD(i:i+mean_step-1) = mean(VAD(i:i+mean_step-1))>0.95;
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