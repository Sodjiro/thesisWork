function [ymotion] = motionframeModif(nameOfFile, threshPer,freqmag) 

v = VideoReader(nameOfFile);

df=floor(v.FrameRate/(2*freqmag));
if df==0 
    df=1;
end
ymotion=0;
nmotion=0;

for r=1:(v.NumberOfFrames-df)
    
    i=read(v,r);
    i=rgb2gray(i);
    j=read(v,r+df);
    j=rgb2gray(j);
    z=imabsdiff(i,j);%IMABSDIFF Absolute difference of two images.
    z=im2bw(z);      %IM2BW Convert image to binary image by thresholding.
    [a b ~]=size(z);
    threshold=a*b*threshPer; % Threshold for motion detection
    
    x=sum(z(:));
     
         
    if(x>threshold)
        ymotion=ymotion+1;
    else
        nmotion=nmotion+1;
        
    end   
end


% fprintf('Motion detected in %d frames.\n', ymotion);
% fprintf('Motion not detected in %d frames.\n', nmotion);
% fprintf('Total number of frames: %d.\n', v.NumberOfFrames);