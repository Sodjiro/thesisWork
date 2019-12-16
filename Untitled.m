clear all

dataDir = './Detected_motions/A1_15_12/2nd_jump_sigma3_alpha_90/';
resultsDir = './Detected_motions/A1_15_12/2nd_jump_sigma3_alpha_90/';

% pyrType = 'octave';
% pyrType = 'halfOctave';
% pyrType = 'smoothHalfOctave';
pyrType = 'quarterOctave';

fileName='VID_20191215_083309.mp4_2nd_95_102.avi';
samplingRate=60;
alpha=90;
sigma=3;
temporalFilter = @FIRWindowBP; 
%%

% magfiles(:,1)=(0:1:samplingRate/2-1)';
% magfiles(:,2)=(1:1:samplingRate/2)';
magfiles(:,1)=[2.9 3.0 3.1];
magfiles(:,2)=[3.0 3.1 3.2];

for alpha=60:10:200
    resultsDir = fullfile(resultsDir,num2str(alpha));
    mkdir(resultsDir);
    freqMag(dataDir, resultsDir,pyrType,fileName,samplingRate,alpha,sigma,temporalFilter,magfiles)
end