
function [] = freqMag(dataDir, resultsDir,pyrType,fileName,samplingRate,alpha,sigma,temporalFilter,magfiles) 

% Paths for the phase-based method
addpath(fullfile(pwd, 'PhaseBased'));
addpath(fullfile(pwd, 'pyrToolsExt'));
addpath(fullfile(pwd, 'Filters'));

%%
inFile = fullfile(dataDir, fileName);

%%
magfiles(magfiles==0)=0.01;
magfiles(magfiles==samplingRate/2)=samplingRate/2-0.01;


for i=1:length(magfiles)
loCutoff = magfiles(i,1);    % Hz
hiCutoff = magfiles(i,2);    % Hz
phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'temporalFilter', temporalFilter,'scaleVideo', 1);
end

for n=1:length(magfiles)
    nameOfFile = sprintf('%0.2f-%0.2f.avi',magfiles(n,1),magfiles(n,2));
    nameOfFile =fullfile(resultsDir, nameOfFile);
    ymotion=motionframeModif(nameOfFile,threshPer,magfiles(n,2));
    magfiles(n,3)=ymotion; %number of motion detected frames
end
%%
dlmwrite(sprintf('%s%s_1filter.txt',resultsDir,fileName),magfiles);
