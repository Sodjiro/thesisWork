v1 = VideoReader('cello.mp4_speed033_1_15.avi');
v2 = VideoReader('cello.mp4_speed033_1_15-FIRWindowBP-band2.50-3.00-sr10-alpha20-mp0-sigma5-scale0.67-frames1-419-halfOctave.avi');
i1 = 0;
i2 = 0;
fig = gcf;
ax1 = subplot(2,2,1, 'Parent', fig);
ax2 = subplot(2,2,2, 'Parent', fig);
while i1 < v1.NumberOfFrames && i2 < v2.NumberOfFrames
        if i1 < v1.NumberOfFrames
            i1 = i1+1;
            if ishandle(ax1)
              image(ax1, v1.read(i1));
            else
              break;    %axes is gone, figure is probably gone too
           end
        end
        if i2 < v2.NumberOfFrames
            i2 = i2+1;
            if ishandle(ax2)
              image(ax2, v2.read(i2));
            else
              break;    %axes is gone, figure is probably gone too
            end
        end
    drawnow
end