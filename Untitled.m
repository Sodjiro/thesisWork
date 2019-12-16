clear all

dataDir = './Detected_motions/A1_15_12/after_1st_jump_sigma3_alpha_90/';
resultsDir = './Detected_motions/A1_15_12/after_1st_jump_sigma3_alpha_90/';

% pyrType = 'octave';
% pyrType = 'halfOctave';
% pyrType = 'smoothHalfOctave';
pyrType = 'quarterOctave';

fileName='VID_20191215_082949.mp4_1st_122_127.avi';
samplingRate=60;
alpha=90;
sigma=3;
temporalFilter = @FIRWindowBP; 
%%

% magfiles(:,1)=(0:1:samplingRate/2-1)';
% magfiles(:,2)=(1:1:samplingRate/2)';
magfiles(:,1)=[3.0];
magfiles(:,2)=[3.3];


freqMag(dataDir, resultsDir,pyrType,fileName,samplingRate,alpha,sigma,temporalFilter,magfiles)
