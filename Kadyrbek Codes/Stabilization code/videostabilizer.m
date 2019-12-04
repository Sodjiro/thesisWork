
%Step 1. Read Frames from a Movie File
% filename = 'shaky_car.avi';
filename = 'ruler_home.avi';


videoFReader = vision.VideoFileReader(filename, ...
                                      'ImageColorSpace', 'Intensity',...
                                      'VideoOutputDataType', 'double');
videoFWriter=vision.VideoFileWriter('stabilizedKD.avi',...
    'FrameRate',videoFReader.info.VideoFrameRate);
hTM = vision.TemplateMatcher('ROIInputPort', true, ...
                            'BestMatchNeighborhoodOutputPort', true);
                 
hVideoOut = vision.VideoPlayer('Name', 'Video Stabilization');
hVideoOut.Position(1) = round(0.4*hVideoOut.Position(1));
hVideoOut.Position(2) = round(1.5*(hVideoOut.Position(2)));
hVideoOut.Position(3:4) = [650 350];

a = VideoReader(strcat(filename));
c = read(a, 1);
sizec=size(c);
[~,RECT] = imcrop(c);% cropping the frame [XMIN YMIN WIDTH HEIGHT]

% RECT=[83.5100,  196.5100,   13.9800,   15.9800];
% RECT=[43,217,20,20];
pos.template_orig = round(RECT(1:2)); % [x y] upper left corner
pos.template_size = round(RECT(3:4));   % [width height]
pos.search_border = round([sizec(1)*0.1 sizec(2)*0.1]);   % max horizontal and vertical displacement
% pos.search_border = [20 20];   % max horizontal and vertical displacement
pos.template_center = floor((pos.template_size-1)/2);
pos.template_center_pos = (pos.template_orig + pos.template_center - 1);
fileInfo = info(videoFReader);
W = fileInfo.VideoSize(1); % Width in pixels
H = fileInfo.VideoSize(2); % Height in pixels
BorderCols = [1:pos.search_border(1)+4 W-pos.search_border(1)+4:W];
BorderRows = [1:pos.search_border(2)+4 H-pos.search_border(2)+4:H];
sz = fileInfo.VideoSize;
TargetRowIndices = ...
  pos.template_orig(2)-1:pos.template_orig(2)+pos.template_size(2)-2;
TargetColIndices = ...
  pos.template_orig(1)-1:pos.template_orig(1)+pos.template_size(1)-2;
SearchRegion = pos.template_orig - pos.search_border - 1;
Offset = [0 0];
Target = zeros(18,22);
firstTime = true;



while ~isDone(videoFReader)
    input = step(videoFReader);
      % Find location of Target in the input video frame
      if firstTime
        Idx = int32(pos.template_center_pos);
        MotionVector = [0 0];
        firstTime = false;
      else
        IdxPrev = Idx;
        ROI = [SearchRegion, pos.template_size+2*pos.search_border];
        Idx = step(hTM, input, Target, ROI);
        MotionVector = double(Idx-IdxPrev);
      end
      [Offset, SearchRegion] = updatesearch(sz, MotionVector, ...
          SearchRegion, Offset, pos);
      % Translate video frame to offset the camera motion
      Stabilized = imtranslate(input,Offset,'linear');
      
      Target = Stabilized(TargetRowIndices, TargetColIndices);
      % Add black border for display
      Stabilized(:, BorderCols) = 0;
      Stabilized(BorderRows, :) = 0;
      TargetRect = [pos.template_orig-Offset, pos.template_size];
      SearchRegionRect = [SearchRegion, pos.template_size + 2*pos.search_border];
      % Draw rectangles on input to show target and search region
      input = insertShape(input,'Rectangle',[TargetRect;SearchRegionRect],...
                          'Color','white');
                      
      % Display the offset(displacement) values on the input image
      txt = sprintf('(%+05.1f,%+05.1f)',Offset);
      input = insertText(input(:,:,1),[191 215],txt,'FontSize',16,...
                        'TextColor','white','BoxOpacity',0);
      
      % Display video
      step(hVideoOut, [input(:,:,1) Stabilized]);
      
      % Save video
      step(videoFWriter,Stabilized); % Save Stabilized video
%       step(videoFWriter, [input(:,:,1) Stabilized]); % Save collage video(Shaky and Stabilized)
     
end


release(videoFReader);
release(videoFWriter);
release(hVideoOut);



function [Offset, SearchRegion] = updatesearch(sz, MotionVector, SearchRegion, Offset, pos)
% Function to update Search Region for SAD and Offset for Translate

  % check bounds
  A_i = Offset - MotionVector;
  AbsTemplate = pos.template_orig - A_i;
  SearchTopLeft = AbsTemplate - pos.search_border;
  SearchBottomRight = SearchTopLeft + (pos.template_size + 2*pos.search_border);

  inbounds = all([(SearchTopLeft >= [1 1]) (SearchBottomRight <= fliplr(sz))]);

  if inbounds
      Mv_out = MotionVector;
  else
      Mv_out = [0 0];
  end

  Offset = Offset - Mv_out;
  SearchRegion = SearchRegion + Mv_out;

end % function updatesearch
