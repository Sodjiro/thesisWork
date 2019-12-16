clear;

% Adds directories to MATLAB path

% Paths for the linear method
addpath(fullfile(pwd, 'Linear'));
addpath(fullfile(pwd, 'Util'));
addpath(fullfile(pwd, 'matlabPyrTools'));
addpath(fullfile(pwd, 'matlabPyrTools', 'MEX'));

% Paths for the phase-based method
addpath(fullfile(pwd, 'PhaseBased'));
addpath(fullfile(pwd, 'pyrToolsExt'));
addpath(fullfile(pwd, 'Filters'));

%%

dataDir = './Detected_motions';

resultsDir = 'Detected_motions/';
% mkdir(resultsDir);

% defaultPyrType = 'octave';
% defaultPyrType = 'halfOctave'; % Half octave pyramid is default as discussed in paper
% defaultPyrType = 'smoothHalfOctave';
defaultPyrType = 'quarterOctave';

fileName='VID_20191215_083309.mp4_2nd_95_102.avi';
inFile = fullfile(dataDir, fileName);
samplingRate = 60; % Hz
alpha = 90;    
sigma = 3;         % Pixels
temporalFilter = @FIRWindowBP; 
pyrType = defaultPyrType;

threshPer=0.000;
incr=0.0001;

% a=VideoReader(inFile);

% while 1
% ymotion=motionframe(inFile,threshPer);
%     if ymotion<a.NumberOfFrames/3
%         break
%     end
% threshPer=threshPer+incr;
% ymotion
% end

fprintf('Threshold is taken as %.3f percent of frame area.\n', threshPer*100);
%%
% magfiles(:,1)=(0:1:samplingRate/2-1)';
% magfiles(:,2)=(1:1:samplingRate/2)';

% magfiles(:,1)=(2:1:9)';
% magfiles(:,2)=(3:1:10)';
% magfiles(:,1)=[2.3 2.4 3.6 3.9 5.3 5.7 6.3 8.1];
% magfiles(:,2)=[2.4 2.7 3.7 4.2 5.4 6.2 6.4 8.5];
% magfiles(:,1)=[8.3 11.3 13.1 14.9];
% magfiles(:,2)=[8.7 11.6 13.4 15.3];
magfiles(:,1)=[2.4 3.9 5.7 8.3];
magfiles(:,2)=[2.7 4.2 6.1 8.5];


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
%     ymotion=motionframe(nameOfFile,threshPer);
    magfiles(n,3)=ymotion; %number of motion detected frames
%     magfiles(n,4)=ymotion;
end
%%
dlmwrite(sprintf('%s%s_1filter.txt',resultsDir,fileName),magfiles);

% medianfreq=median(unique(magfiles(:,3))); % Use for threshold
% medianfreq=a.NumberOfFrames/10;
medianfreq=0;
idx = magfiles(:,3)<medianfreq;
magfiles(idx,3) = 0;
dlmwrite(sprintf('%s%s_1filter0.txt',resultsDir,fileName),magfiles);

for n=1:length(magfiles(:,1))
    if magfiles(n,3)>0
        for i=magfiles(n,1):(magfiles(n,2)-magfiles(n,1))/10:magfiles(n,2)-(magfiles(n,2)-magfiles(n,1))/10
        loCutoff = i;    % Hz
        hiCutoff = i+(magfiles(n,2)-magfiles(n,1))/10;    % Hz
        phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'temporalFilter', temporalFilter,'scaleVideo', 1);
        end
    end
end

for n=1:length(magfiles(:,1))
    if magfiles(n,3)>0
        b=1;
        for i=magfiles(n,1):(magfiles(n,2)-magfiles(n,1))/10:magfiles(n,2)-(magfiles(n,2)-magfiles(n,1))/10
        nameOfFile = sprintf('%0.2f-%0.2f.avi',i,i+(magfiles(n,2)-magfiles(n,1))/10);
        nameOfFile =fullfile(resultsDir, nameOfFile);
        ymotion=motionframeModif(nameOfFile,threshPer,i+(magfiles(n,2)-magfiles(n,1))/10);
        magfiles2(b,n)=ymotion;
        b=b+1;
        end
    end
end



magfiles2(magfiles2==0) = NaN;
[M,IND] = max(magfiles2,[],'omitnan');
for i=1:1:length(IND)
     if IND(i)>1 && IND(i)<length(magfiles2(:,1))
         freq(i,1)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i)-1);
         freq(i,2)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i));
         
     elseif IND(i)==length(magfiles2(:,1)) && i~=length(IND)
         if M(i)<M(i+1)
            freq(i+1,1)= magfiles(i+1,1)+(magfiles(i+1,2)-magfiles(i+1,1))/10*(IND(i+1)-1);
            freq(i+1,2)= magfiles(i+1,1)+(magfiles(i+1,2)-magfiles(i+1,1))/10*(IND(i+1));
         else
            freq(i,1)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i)-1);
            freq(i,2)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i));
         end
     
     elseif IND(i)==length(magfiles2(:,1)) && i==length(IND)
         freq(i,1)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i)-1);
         freq(i,2)= magfiles(i,1)+(magfiles(i,2)-magfiles(i,1))/10*(IND(i));
     end
end

freq=round(freq,2);
freq( all(~freq,2), : ) = [];

dlmwrite(sprintf('%s%s_2filter.txt',resultsDir,fileName),magfiles2);
dlmwrite(sprintf('%s%s_freq.txt',resultsDir,fileName),freq);

%% Functions
% 
% function [ymotion] = motionframeModif(nameOfFile, threshPer,freqmag) 
% 
% v = VideoReader(nameOfFile);
% 
% df=floor(v.FrameRate/(2*freqmag));
% if df==0 
%     df=1;
% end
% ymotion=0;
% nmotion=0;
% 
% for r=1:(v.NumberOfFrames-df)
%     
%     i=read(v,r);
%     i=rgb2gray(i);
%     j=read(v,r+df);
%     j=rgb2gray(j);
%     z=imabsdiff(i,j);%IMABSDIFF Absolute difference of two images.
%     z=im2bw(z);      %IM2BW Convert image to binary image by thresholding.
%     [a b ~]=size(z);
%     threshold=a*b*threshPer; % Threshold for motion detection
%     
%     x=sum(z(:));
%      
%          
%     if(x>threshold)
%         ymotion=ymotion+1;
%     else
%         nmotion=nmotion+1;
%         
%     end   
% end
% end
% function [ymotion] = motionframe(nameOfFile, threshPer) 
% 
% v = VideoReader(nameOfFile);
% 
% ymotion=0;
% nmotion=0;
% 
% for r=1:(v.NumberOfFrames-1)
%     
%     i=read(v,r);
%     i=rgb2gray(i);
%     j=read(v,r+1);
%     j=rgb2gray(j);
%     z=imabsdiff(i,j);%IMABSDIFF Absolute difference of two images.
%     z=im2bw(z);      %IM2BW Convert image to binary image by thresholding.
%     [a b ~]=size(z);
%     threshold=a*b*threshPer; % Threshold for motion detection
%     
%     x=sum(z(:));
%      
%          
%     if(x>threshold)
%         ymotion=ymotion+1;
%     else
%         nmotion=nmotion+1;
%         
%     end   
% end
% end
% 
