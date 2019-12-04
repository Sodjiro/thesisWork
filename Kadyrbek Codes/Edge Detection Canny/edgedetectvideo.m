function [] = edgedetectvideo(nameOfFile,T_Low,T_High) 
inputName = strcat(nameOfFile);
outputName = strcat(nameOfFile,'_canny_',num2str(round(T_Low,3)),'_',num2str(round(T_High,3)),'.avi');

a = VideoReader(inputName);
vidObj = VideoWriter(outputName);
vidObj.Quality = 100;
vidObj.FrameRate = a.FrameRate;
open(vidObj);

for img = 1:a.NumberOfFrames
    b = read(a, img);
    b2 = edgedetect(b,T_Low,T_High);
    writeVideo(vidObj,b2)
    
end
close(vidObj);

