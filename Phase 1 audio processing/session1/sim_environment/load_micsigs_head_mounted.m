% load_micsigs_head_mounted.m

fs = 44100;
data_length = 10;

impulse_path        = '../Head_mounted_real/';
% positions_filename_header   = ["90","-60","-90"];
positions_filename_header   = ["-90"];
% type_filename_header = "noise";
% type_filename_header = "part1_track1";
type_filename_header = "part2_track1";
L1_filename_header = 'HML_M1';
L2_filename_header = 'HML_M2';
R1_filename_header = 'HMR_M1';
R2_filename_header = 'HMR_M2';



Mic = zeros(    length(positions_filename_header), ...
                4, ...
                data_length*fs);

for i=1:length(positions_filename_header)
    [audio, Fs]  = audioread(strcat(impulse_path, L1_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,1,:) = audio(1:data_length*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, L2_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,2,:) = audio(1:data_length*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, R1_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,3,:) = audio(1:data_length*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, R2_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,4,:) = audio(1:data_length*fs);
end
