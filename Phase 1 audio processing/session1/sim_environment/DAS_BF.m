%% Setup
clear all;
run create_micsigs.m;

%% Estimate DOA using MUSIC_narrrowband
run MUSIC_narrowband.m;

%% Start of Beemforming

Q = size(RIR_sources,3);

delay = abs(cosd(DOA_est)*mic_distance/100);

speech_DAS = zeros(size(Mic,3),1);
noise_DAS = zeros(size(noise,3),1);

if(DOA_est<90)
    for i=1:size(RIR_sources,2)
        speech = Mic(1,i,:);
        speech = reshape(speech,size(speech_DAS));
        speech_DAS = i*delay*speech+ speech_DAS;
        
        noise = reshape(noise,size(noise_DAS));
        noise_DAS = i*delay*noise + noise_DAS;
    end
elseif(DOA_est>90)
     for i=size(RIR_sources,2):-1:1
        speech = Mic(1,i,:);
        speech = reshape(speech,size(speech_DAS));
        speech_DAS = i*delay*speech+ speech_DAS;
        
        noise = reshape(noise,size(noise_DAS));
        noise_DAS = i*delay*noise + noise_DAS;
    end
else
    print 'you are doomed';
    
end

speech_DAS = speech_DAS/size(RIR_sources,2);
noise_DAS = noise_DAS/size(RIR_sources,2);

DAS_out = 10*(speech_DAS+noise_DAS);
figure
plot(sound);
hold on;
plot(DAS_out);

soundsc(DAS_out,fs_RIR)
speech_power_DAS = var(speech_DAS);
noise_power_DAS= var(noise_DAS);

SNR_DAS = 10*log10(speech_power_DAS./noise_power_DAS);


