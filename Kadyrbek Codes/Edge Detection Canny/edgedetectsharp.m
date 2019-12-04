function Canny_img=edgedetectsharp(rgbImage,sharpnum)

grayImage = rgb2gray(rgbImage);

f1=-1/256*[1 4 6 4 1, 4 16 24 16 4, 6 24 -476 24 6, 4 16 24 16 4,1 4 6 4 1];

for i=1:sharpnum
grayImage = filter2(f1,grayImage);
end
Canny_img = edge(grayImage, 'Canny');