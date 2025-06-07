

% if exist([cd,'\images'],'file')
%     rmdir([cd,'\images'],'s')
% end
% mkdir([cd,'\images'])
% 
% 
% for i = 3:105
%    
%     imwrite(F(i).cdata,[cd,'\images\p',num2str(i-2),'.png']);
% end

v = VideoWriter('emdlab');


v.open;

for i = 2:numel(F)
%     im = imread(['p',num2str(i),'.jpg']);
    v.writeVideo(F(i).cdata);
end

v.close
v.FrameRate = 1;
v.saveobj




