clear all
close all
clc

addpath /home/liuzzil2/fieldtrip-20190812/
ft_defaults


%% Co-register MRI from fiducial positions
subn = ['24071' ; '24172'; '24138'; '24103'; '23490';
    '24213'; '24201' ; '23911'; '24208'; '24199';
    '22695'; '22694'; '24175'; '24216'; '23732'; '23951'];

for sn = 1 %[1,3,4,6,7,8,9,11,14,15,16] %Subjects showing enough variation in mood
    sub = subn(sn,:);
    freqband = [40 90];
% Checekd on sub1: max(svd(C)) = 1e-25 for all freqs, min(svd(C)) becomes 
% much smaller at higher freqs! 1e-31 for theta, 1e-40 for beta, 1e-44 for gamma
    mu = 0.001; %mu = 1e-14 = 4*min(svd(C_beta)) 
    
%     mmi_mood_beamformer(sub,freqband,mu); 
    mmi_buttonpress_beamformer(sub,freqband,0.15,mu); 
    
end