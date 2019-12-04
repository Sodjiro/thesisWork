function [] = extractscene(nameOfFile, beginTime, endTime, SceneClassifier) 
inputName = strcat(nameOfFile);
outputName = strcat(nameOfFile,'_',num2str(SceneClassifier),'_',num2str(beginTime),'_',num2str(endTime),'.avi');
a = VideoReader(inputName);
beginFrame = round(beginTime * a.FrameRate);
endFrame = round(endTime * a.FrameRate);
vidObj = VideoWriter(outputName);
vidObj.Quality = 100;
vidObj.FrameRate = a.FrameRate;
open(vidObj);
c = read(a, beginFrame);
[~,RECT] = imcrop(c);% cropping the frame [XMIN YMIN WIDTH HEIGHT]

for img = beginFrame:endFrame
    b = read(a, img);
    b2 = imcrop(b,RECT);% cropping the frame [XMIN YMIN WIDTH HEIGHT]
    writeVideo(vidObj,b2)
    
end
close(vidObj);
