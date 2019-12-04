threshPer=0.001;
incr=0.001;



while 1
nameOfFile = strcat('0.avi');
ymotion=motionframe(nameOfFile,threshPer);
    if ymotion==0
        break
    end
threshPer=threshPer+incr;
ymotion
end

fprintf('Threshold is taken as %.3f percent of frame area.\n', threshPer*100);

for n=1:1:16
    nameOfFile = strcat(num2str(n),'.avi');
    ymotion=motionframe(nameOfFile,threshPer);
    detected(n)=ymotion; %number of motion detected frames
end

function [ymotion] = motionframe(nameOfFile, threshPer) 

v = VideoReader(nameOfFile);

ymotion=0;
nmotion=0;

for r=1:(v.NumberOfFrames-1)
    
    i=read(v,r);
    i=rgb2gray(i);
    j=read(v,r+1);
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
end
