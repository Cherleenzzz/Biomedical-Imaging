function   SaveAsVideo(F,name)
% F is the frame captured,e.g.:
% F(size) = struct('cdata',[],'colormap',[]);
% F(i) = getframe(gcf);
% Check qu6_5_1.m line 142 and 150

% name is the video name include mp4, e.g. 'myVideo.mp4'

% create the video writer with 3 fps
  writerObj = VideoWriter(name, 'MPEG-4');
  writerObj.FrameRate = 3;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
end

