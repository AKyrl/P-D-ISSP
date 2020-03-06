% load_micsigs_linear_array.m

impulse_path        = '../Head_mounted_real/';
% positions_filename_header   = ["90","-60","-90"];
% positions_filename_header   = ["-60"];
% % type_filename_header = "noise";
% % type_filename_header = "part1_track1";
% type_filename_header = "part2_track1";
% type_filename_header = "part2_track2";
M1_filename_header = 'LMA_M1';
M2_filename_header = 'LMA_M2';
M3_filename_header = 'LMA_M3';
M4_filename_header = 'LMA_M4';



Mic = zeros(    length(positions_filename_header), ...
                4, ...
                length_limit*fs);

for i=1:length(positions_filename_header)
    [audio, Fs]  = audioread(strcat(impulse_path, M1_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,1,:) = audio(1:length_limit*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, M2_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,2,:) = audio(1:length_limit*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, M3_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,3,:) = audio(1:length_limit*fs);
    [audio, Fs]  = audioread(strcat(impulse_path, M4_filename_header,'_',positions_filename_header(i),'_',type_filename_header,'.wav'));
    Mic(i,4,:) = audio(1:length_limit*fs);
end
