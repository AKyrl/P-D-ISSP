
%% Start of Beemforming

Q = size(Mic,2);

delay = abs(cosd(DOA_est)*mic_distance/100/340);

samples = ceil(delay*fs_RIR);

speech = zeros(size(Mic,2),size(Mic,3));
speech_DAS = zeros(1,size(Mic,3));
% noise_DAS = zeros(1,size(noise,3));

if(DOA_est<90)
    for mic_idx=1:size(Mic,2)
        speech_temp = Mic(1,mic_idx,:);
        speech_temp = reshape(speech_temp,size(speech_DAS));
        speech(mic_idx, 1:end-(mic_idx-1)*samples) = speech_temp((mic_idx-1)*samples+1:end);
        speech_DAS(1:end-(mic_idx-1)*samples) = speech(mic_idx,1:end-(mic_idx-1)*samples)+ speech_DAS(1:end-(mic_idx-1)*samples);
    end
elseif(DOA_est>90)
     for mic_idx=1:size(Mic,2)
         count = size(RIR_sources,2) - mic_idx+1;
        speech_temp = Mic(1,mic_idx,:);
        speech_temp = reshape(speech_temp,size(speech_DAS));
        speech(mic_idx, 1:end-(count-1)*samples) = speech_temp((count-1)*samples+1:end);
        speech_DAS(1:end-(count-1)*samples) = speech(mic_idx,1:end-(count-1)*samples)+ speech_DAS(1:end-(count-1)*samples);
    end
else
    print 'you are doomed';
    
end

speech_DAS = speech_DAS/size(RIR_sources,2);
% noise_DAS = noise_DAS/size(RIR_sources,2);

DAS_out = 2*speech_DAS;
% figure
% plot(reshape(Mic(1,1,:),size(speech_DAS)));
% hold on;
% plot(speech_DAS);

% % soundsc(DAS_out,fs_RIR)
% VAD=abs(speech_DAS)>std(speech_DAS)*1e-3;
% speech_power_DAS = var(speech_DAS(VAD==1));
% % noise_power_DAS= var(noise_DAS(VAD==1));
speech_power_DAS = var(speech_DAS);
noise_power_DAS= var(noise_DAS);


for channel_idx=1:size(Mic(:,:,:),1)
    for mic_idx=1:size(Mic,2)
        speech_in_DAS(channel_idx,mic_idx) = var(Mic(channel_idx,mic_idx,:),0,'all');
    end
end

% SNR_DAS = 10*log10(speech_power_DAS./noise_power_DAS);


