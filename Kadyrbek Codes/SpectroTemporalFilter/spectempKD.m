


a = VideoReader('data/a4(4)crop.avi');
b = read(a, 1);
b2 = imcrop(b, [390,130,20,0]);
b3=zeros(size(b2));

for img = 1:60
    b = read(a, img);
    b2 = imcrop(b, [390,130,20,0]);
    b3=[b3;b2];
end

% Read in original RGB image.
rgbImage  = b3;
subplot(2, 2, 1);
imshow(rgbImage)
axis('on', 'image');
title('Original Image')
% Convert to gray scale.

grayImage1 = rgb2gray(rgbImage);

% min33 = ordfilt2(grayImage1,1,ones(3,3));
% grayImage = ordfilt2(min33,9,ones(3,3));
% grayImage  = imsharpen(grayImage1,'Radius',1.2,'Amount',2,'Threshold',0.);

f1=-1/256*[1 4 6 4 1, 4 16 24 16 4, 6 24 -476 24 6, 4 16 24 16 4,1 4 6 4 1];
grayImage = filter2(f1,grayImage1);

subplot(2, 2, 2);
imshow(grayImage/(max(max(grayImage))))
title('Sharp Image')
Canny_img = edge(grayImage1, 'Canny');
subplot(2, 2, 3);
imshow(Canny_img, [])
axis('on', 'image');
title('Edge Detected Image of orig. image')
% Get edges
Canny_img = edge(grayImage, 'Canny');
subplot(2, 2, 4);
imshow(Canny_img, [])
axis('on', 'image');
title('Edge Detected Image of sharp image')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0.05, 1, 0.95]);
