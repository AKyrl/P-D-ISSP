%% Setup
clear all;
run create_micsigs.m;

%% Estimate DOA using MUSIC_narrrowband
run MUSIC_narrowband.m;

%% Start of Beemforming

Q = size(RIR_sources,3);

delay = abs(cosd(DOA_est)*mic_distance/100/340);

samples = ceil(delay*fs_RIR);

speech_DAS = zeros(size(Mic,3),1);
noise_DAS = zeros(size(noise,3),1);

if(DOA_est<90)
    for i=1:size(RIR_sources,2)
        speech = Mic(1,i,:);
        speech = reshape(speech,size(speech_DAS));
        speech_DAS(1:end-(i-1)*samples) = speech((i-1)*samples+1:end)+ speech_DAS(1:end-(i-1)*samples);
        
        noise = reshape(noise,size(noise_DAS));
        noise_DAS(1:end-(i-1)*samples) = noise((i-1)*samples+1:end)+ noise_DAS(1:end-(i-1)*samples);
    end
elseif(DOA_est>90)
     for i=1:size(RIR_sources,2)
         count = size(RIR_sources,2) - i+1;
        speech = Mic(1,i,:);
        speech = reshape(speech,size(speech_DAS));
        speech_DAS(1:end-(count-1)*samples) = speech((count-1)*samples+1:end)+ speech_DAS(1:end-(count-1)*samples);
        
        noise = reshape(noise,size(noise_DAS));
        noise_DAS(1:end-(count-1)*samples) = noise((count-1)*samples+1:end)+ noise_DAS(1:end-(count-1)*samples);
        
    end
else
    print 'you are doomed';
    
end

speech_DAS = speech_DAS/size(RIR_sources,2);
noise_DAS = noise_DAS/size(RIR_sources,2);

DAS_out = 2*(speech_DAS+noise_DAS);
figure
plot(reshape(Mic(1,1,:),size(speech_DAS)));
hold on;
plot(speech_DAS);

soundsc(DAS_out,fs_RIR)
speech_power_DAS = var(speech_DAS);
noise_power_DAS= var(noise_DAS);

SNR_DAS = 10*log10(speech_power_DAS./noise_power_DAS);


