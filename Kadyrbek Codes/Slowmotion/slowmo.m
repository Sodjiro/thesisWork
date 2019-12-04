% MATLAB program to convert video into slow motion 
% clc;clear;close all; 
  
 % load the video. 
obj = VideoReader('a4(1)crop_15_20sr50-FIRWindowBP-band21.00-22.00-sr50-alpha20-mp0-sigma5-scale0.67-frames1-250-halfOctave-STAR.avi');   
  
% Write in new variable 
obj2= VideoWriter('slowmo2.avi');     
  
% decrease framerate  
obj2.FrameRate = 5;               
open(obj2); 
  
% for reading frames one by one 
while hasFrame(obj)               
    k = readFrame(obj);  
  
    % write the frames in obj2.          
    obj2.writeVideo(k);           
end
  
close(obj2); 