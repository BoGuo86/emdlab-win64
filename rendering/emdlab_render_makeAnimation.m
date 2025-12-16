function emdlab_render_makeAnimation(path, NFrames, delayTime)

% set default speed
if nargin<3, delayTime = 0.1; end

filename = 'animation.gif';  % Output GIF file

for k = 1:NFrames
    % Read the image
    img = imread(sprintf('%s\\frame-%d.png', path, k));  % or use full path if needed

    % Convert to indexed image
    [imind, cm] = rgb2ind(img, 256);

    % Write to the GIF file
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delayTime);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

delete(sprintf('%s\\frame-*.png', path))

end